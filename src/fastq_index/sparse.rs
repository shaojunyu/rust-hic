//! Streaming filter for indexed_gzip GZIDX v1. Memory is proportional to
//! checkpoint metadata, not window data; unselected windows never hit disk.
use std::{
    fs::File,
    io::{self, Write},
    path::Path,
};
/// Keep roughly 64 MiB decoded spacing and the last point before each chunk.
pub(crate) const CHECKPOINT_SPACING: u64 = 64 * 1024 * 1024;
pub(crate) const MAX_CHECKPOINT_GAP: u64 = 128 * 1024 * 1024;
pub(crate) fn select_points(positions: &[u64], offsets: &[u64]) -> io::Result<Vec<bool>> {
    if positions.is_empty() || positions.windows(2).any(|w| w[0] > w[1]) {
        return Err(io::Error::other("invalid or unsorted checkpoints"));
    }
    let mut selected = vec![false; positions.len()];
    selected[0] = true;
    selected[positions.len() - 1] = true;
    let mut last = positions[0];
    for (i, &position) in positions.iter().enumerate() {
        if position.saturating_sub(last) >= CHECKPOINT_SPACING {
            selected[i] = true;
            last = position;
        }
    }
    for offset in offsets {
        let n = positions
            .partition_point(|p| p <= offset)
            .checked_sub(1)
            .ok_or_else(|| io::Error::other("no preceding checkpoint"))?;
        selected[n] = true;
    }
    Ok(selected)
}
pub struct SparseWriter {
    file: File,
    offsets: Vec<u64>,
    header: Vec<u8>,
    header_len: Option<usize>,
    selected: Vec<bool>,
    rows: Vec<[u8; 18]>,
    window: usize,
    row: usize,
    within: usize,
}
impl SparseWriter {
    pub fn new(path: &Path, offsets: &[u64]) -> io::Result<Self> {
        Ok(Self {
            file: File::create(path)?,
            offsets: offsets.to_vec(),
            header: Vec::new(),
            header_len: None,
            selected: Vec::new(),
            rows: Vec::new(),
            window: 0,
            row: 0,
            within: 0,
        })
    }
    fn parse_header(&mut self) -> io::Result<()> {
        if &self.header[..5] != b"GZIDX" || self.header[5] != 1 {
            return Err(io::Error::other("unsupported GZIDX"));
        }
        let count = u32::from_le_bytes(self.header[31..35].try_into().unwrap()) as usize;
        self.window = u32::from_le_bytes(self.header[27..31].try_into().unwrap()) as usize;
        if count == 0 || count > 10_000_000 || self.window > 1_048_576 {
            return Err(io::Error::other("invalid GZIDX header"));
        }
        self.rows = self.header[35..]
            .chunks_exact(18)
            .map(|x| x.try_into().unwrap())
            .collect();
        if self.rows.len() != count {
            return Err(io::Error::other("wrong checkpoint count"));
        }
        let positions: Vec<u64> = self
            .rows
            .iter()
            .map(|r| u64::from_le_bytes(r[8..16].try_into().unwrap()))
            .collect();
        self.selected = select_points(&positions, &self.offsets)?;
        let n = self.selected.iter().filter(|x| **x).count() as u32;
        let mut head = self.header[..35].to_vec();
        head[31..35].copy_from_slice(&n.to_le_bytes());
        self.file.write_all(&head)?;
        for (row, keep) in self.rows.iter().zip(&self.selected) {
            if *keep {
                self.file.write_all(row)?;
            }
        }
        self.header.clear();
        self.next_window();
        Ok(())
    }
    fn next_window(&mut self) {
        while self.row < self.rows.len() && self.rows[self.row][17] == 0 {
            self.row += 1;
        }
    }
    pub fn finish(mut self) -> io::Result<()> {
        if self.rows.is_empty() || self.row != self.rows.len() || self.within != 0 {
            return Err(io::Error::other("truncated gzip index stream"));
        }
        self.file.flush()
    }
}
impl Write for SparseWriter {
    fn write(&mut self, mut data: &[u8]) -> io::Result<usize> {
        let total = data.len();
        while !data.is_empty() {
            if self.rows.is_empty() {
                let need = self.header_len.unwrap_or(35);
                let n = (need - self.header.len()).min(data.len());
                self.header.extend_from_slice(&data[..n]);
                data = &data[n..];
                if self.header.len() == 35 && self.header_len.is_none() {
                    let count =
                        u32::from_le_bytes(self.header[31..35].try_into().unwrap()) as usize;
                    if count == 0 || count > 10_000_000 {
                        return Err(io::Error::other("invalid checkpoint count"));
                    }
                    self.header_len = Some(35 + 18 * count);
                }
                if self.header_len == Some(self.header.len()) {
                    self.parse_header()?;
                }
            } else {
                if self.row >= self.rows.len() {
                    return Err(io::Error::other("extra gzip index data"));
                }
                let n = (self.window - self.within).min(data.len());
                if self.selected[self.row] {
                    self.file.write_all(&data[..n])?;
                }
                self.within += n;
                data = &data[n..];
                if self.within == self.window {
                    self.within = 0;
                    self.row += 1;
                    self.next_window();
                }
            }
        }
        Ok(total)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.file.flush()
    }
}
