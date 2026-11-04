//! BGZF pairs and PX2 indexes. Format follows pairix 0.3.9 (MIT license),
//! copyright Soo Lee; see licenses/pairix-MIT.txt.
use anyhow::{bail, Context, Result};
use flate2::{write::DeflateEncoder, Compression};
use rayon::prelude::*;
use std::{
    collections::HashSet,
    fs::{File, OpenOptions},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    sync::mpsc::{sync_channel, SyncSender},
    thread::{self, JoinHandle},
    time::Instant,
};

const BLOCK: usize = 65_280;
struct Bgzf {
    buffer: Vec<u8>,
    uncompressed: u64,
    batch_size: usize,
    sender: Option<SyncSender<Vec<u8>>>,
    worker: Option<JoinHandle<Result<Vec<u64>>>>,
}
impl Bgzf {
    fn new(mut out: BufWriter<File>, threads: usize) -> Result<Self> {
        if !(1..=64).contains(&threads) {
            bail!("compression threads must be 1..64")
        }
        let batch_size = BLOCK * threads * 2;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()?;
        let (tx, rx) = sync_channel::<Vec<u8>>(2);
        let worker = thread::spawn(move || -> Result<Vec<u64>> {
            let mut compressed = 0u64;
            let mut offsets = Vec::new();
            while let Ok(batch) = rx.recv() {
                let blocks = pool.install(|| {
                    batch
                        .par_chunks(BLOCK)
                        .map(Self::encode_block)
                        .collect::<Result<Vec<_>>>()
                })?;
                for block in blocks {
                    if compressed >= (1u64 << 48) {
                        bail!("BGZF file exceeds virtual offset limit")
                    }
                    offsets.push(compressed);
                    out.write_all(&block)?;
                    compressed = compressed
                        .checked_add(block.len() as u64)
                        .context("compressed offset overflow")?;
                }
            }
            offsets.push(compressed);
            out.write_all(&Self::encode_block(&[])?)?;
            out.flush()?;
            Ok(offsets)
        });
        Ok(Self {
            buffer: Vec::with_capacity(batch_size),
            uncompressed: 0,
            batch_size,
            sender: Some(tx),
            worker: Some(worker),
        })
    }
    fn position(&self) -> Result<u64> {
        Ok(self.uncompressed)
    }
    fn encode_block(data: &[u8]) -> Result<Vec<u8>> {
        let mut encoder = DeflateEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(data)?;
        let compressed = encoder.finish()?;
        let size = 18 + compressed.len() + 8;
        if size > 65_536 {
            bail!("BGZF compressed block too large")
        }
        let mut header: [u8; 18] = [
            31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, b'B', b'C', 2, 0, 0, 0,
        ];
        header[16..18].copy_from_slice(&((size - 1) as u16).to_le_bytes());
        let mut result = Vec::with_capacity(size);
        result.extend_from_slice(&header);
        result.extend_from_slice(&compressed);
        result.extend_from_slice(&crc32fast::hash(data).to_le_bytes());
        result.extend_from_slice(&(data.len() as u32).to_le_bytes());
        Ok(result)
    }
    fn submit_blocks(&mut self) -> Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let batch = std::mem::replace(&mut self.buffer, Vec::with_capacity(self.batch_size));
        self.sender
            .as_ref()
            .context("BGZF writer is closed")?
            .send(batch)
            .map_err(|_| anyhow::anyhow!("BGZF compression worker stopped"))?;
        Ok(())
    }
    fn write_bytes(&mut self, mut data: &[u8]) -> Result<()> {
        while !data.is_empty() {
            let n = data.len().min(self.batch_size - self.buffer.len());
            self.buffer.extend_from_slice(&data[..n]);
            self.uncompressed = self
                .uncompressed
                .checked_add(n as u64)
                .context("uncompressed offset overflow")?;
            data = &data[n..];
            if self.buffer.len() == self.batch_size {
                self.submit_blocks()?
            }
        }
        Ok(())
    }
    fn finish(mut self) -> Result<Vec<u64>> {
        self.submit_blocks()?;
        drop(self.sender.take());
        self.worker
            .take()
            .context("BGZF worker missing")?
            .join()
            .map_err(|_| anyhow::anyhow!("BGZF compression worker panicked"))?
    }
}
impl Drop for Bgzf {
    fn drop(&mut self) {
        drop(self.sender.take());
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}
#[derive(Default)]
struct Group {
    name: String,
    chrom1: String,
    chrom2: String,
    bins: Vec<(u32, u64, u64)>,
    linear: Vec<u64>,
}
fn u32write(out: &mut Bgzf, x: u32) -> Result<()> {
    out.write_bytes(&x.to_le_bytes())
}
fn u64write(out: &mut Bgzf, x: u64) -> Result<()> {
    out.write_bytes(&x.to_le_bytes())
}
fn new_file(path: &Path) -> Result<BufWriter<File>> {
    Ok(BufWriter::with_capacity(
        1024 * 1024,
        OpenOptions::new().write(true).create_new(true).open(path)?,
    ))
}
fn write_index(path: &Path, groups: &[Group], lines: u64) -> Result<()> {
    let mut out = Bgzf::new(new_file(path)?, 1)?;
    out.write_bytes(b"PX2.004\x01")?;
    u32write(&mut out, groups.len().try_into()?)?;
    u64write(&mut out, lines)?;
    // Pairix preset: 1-based pairs, columns 2,3,3,4,5,5 and tab delimiter.
    for n in [3, 2, 3, 3, 4, 5, 5] {
        u32write(&mut out, n)?;
    }
    out.write_bytes(&[b'\t', b'|', 0, 0])?;
    u32write(&mut out, b'#' as u32)?;
    u32write(&mut out, 0)?;
    let names = groups.iter().try_fold(0usize, |n, g| {
        n.checked_add(g.name.len() + 1)
            .context("index names overflow")
    })?;
    u32write(&mut out, names.try_into()?)?;
    for group in groups {
        out.write_bytes(group.name.as_bytes())?;
        out.write_bytes(&[0])?;
    }
    for group in groups {
        u32write(&mut out, group.bins.len().try_into()?)?;
        for &(bin, start, end) in &group.bins {
            u32write(&mut out, bin)?;
            u32write(&mut out, 1)?;
            u64write(&mut out, start)?;
            u64write(&mut out, end)?;
        }
        u32write(&mut out, group.linear.len().try_into()?)?;
        let mut last = 0;
        for &offset in &group.linear {
            if offset != 0 {
                last = offset
            };
            u64write(&mut out, last)?;
        }
    }
    out.finish()?;
    Ok(())
}
pub fn pack(input: &Path, output: &Path, threads: usize, log: Option<&Path>) -> Result<()> {
    let started = Instant::now();
    if !(1..=64).contains(&threads) {
        bail!("compression threads must be 1..64")
    }
    if !output.to_string_lossy().ends_with(".pairs.gz") {
        bail!("output must end in .pairs.gz")
    }
    let mut index_name = output.as_os_str().to_os_string();
    index_name.push(".px2");
    let index = PathBuf::from(index_name);
    if output.exists() || index.exists() {
        bail!("output or pairix index already exists")
    }
    let mut log = log
        .map(|path| OpenOptions::new().write(true).create_new(true).open(path))
        .transpose()?;
    let mut reader = BufReader::with_capacity(1024 * 1024, File::open(input)?);
    let mut out = Bgzf::new(new_file(output)?, threads)?;
    let mut line = String::new();
    let mut groups = Vec::<Group>::new();
    let mut seen = HashSet::new();
    let (mut lines, mut records, mut previous_position) = (0u64, 0u64, 0u64);
    let (mut body, mut columns) = (false, None);
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if !line.ends_with('\n') {
            bail!("pairs input has an incomplete final line")
        }
        if line.as_bytes().contains(&0) {
            bail!("NUL bytes are unsupported in pairs files")
        }
        lines += 1;
        let offset = out.position()?;
        if line.starts_with('#') {
            if body {
                bail!("header after pairs body")
            }
            if line.starts_with("#columns:") {
                if columns.is_some() {
                    bail!("duplicate #columns header")
                }
                columns = Some(crate::pairs::format::Format::from_columns(&line)?);
            }
        } else {
            if columns.is_none() {
                bail!("missing #columns header")
            }
            body = true;
            let body_line = line.trim_end_matches(['\r', '\n']);
            let mut fields = body_line.split('\t');
            let _read_id = fields.next().context("missing read ID")?;
            let chrom1 = fields.next().context("missing chromosome 1")?;
            let pos1 = fields.next().context("missing position 1")?;
            let chrom2 = fields.next().context("missing chromosome 2")?;
            let pos2 = fields.next().context("missing position 2")?;
            if fields.count() + 5 != columns.unwrap().len() {
                bail!("pairs row does not match #columns header")
            }
            if [chrom1, chrom2].iter().any(|name| {
                name.is_empty()
                    || *name == "!"
                    || name.contains('|')
                    || name.chars().any(|c| c.is_whitespace() || c == '\0')
            }) {
                bail!("unsupported chromosome name")
            }
            if [pos1, pos2]
                .iter()
                .any(|pos| pos.len() > 1 && pos.starts_with('0'))
            {
                bail!("positions must be decimal without leading zeros")
            }
            let a: u64 = pos1.parse()?;
            let b: u64 = pos2.parse()?;
            if a == 0 || b == 0 || a > (1u64 << 30) || b > (1u64 << 30) {
                bail!("pairix supports 1-based positions in 1..=2^30")
            }
            if groups
                .last()
                .is_none_or(|g| g.chrom1 != chrom1 || g.chrom2 != chrom2)
            {
                let name = format!("{chrom1}|{chrom2}");
                if !seen.insert(name.clone()) {
                    bail!("chromosome pair occurs in multiple blocks")
                }
                if let Some(g) = groups.last_mut() {
                    if let Some(last) = g.bins.last_mut() {
                        last.2 = offset;
                    }
                }
                groups.push(Group {
                    name,
                    chrom1: chrom1.to_owned(),
                    chrom2: chrom2.to_owned(),
                    ..Default::default()
                });
                previous_position = 0;
            }
            if a < previous_position {
                bail!("positions out of order within chromosome pair")
            }
            previous_position = a;
            let g = groups.last_mut().unwrap();
            let tile = ((a - 1) >> 15) as usize;
            let bin = 4681 + tile as u32;
            if g.linear.len() <= tile {
                g.linear.resize(tile + 1, 0);
            }
            if g.linear[tile] == 0 {
                g.linear[tile] = offset;
            }
            if g.bins.last().is_none_or(|last| last.0 != bin) {
                if let Some(last) = g.bins.last_mut() {
                    last.2 = offset;
                }
                g.bins.push((bin, offset, 0));
            }
            records += 1;
            if records % 10_000_000 == 0 {
                let event = format!(
                    "event=pairix_progress records={records} input_bytes={} elapsed_seconds={:.3}",
                    offset,
                    started.elapsed().as_secs_f64()
                );
                eprintln!("{event}");
                if let Some(log) = &mut log {
                    writeln!(log, "{event}")?;
                }
            }
        }
        out.write_bytes(line.as_bytes())?;
    }
    if columns.is_none() {
        bail!("missing #columns header")
    }
    if let Some(g) = groups.last_mut() {
        if let Some(last) = g.bins.last_mut() {
            last.2 = out.position()?;
        }
    }
    let offsets = out.finish()?;
    let virtual_offset = |offset: u64| -> Result<u64> {
        let block: usize = (offset / BLOCK as u64).try_into()?;
        let compressed = *offsets.get(block).context("missing BGZF block offset")?;
        if compressed >= (1u64 << 48) {
            bail!("BGZF file exceeds virtual offset limit")
        }
        Ok((compressed << 16) | (offset % BLOCK as u64))
    };
    for group in &mut groups {
        for (_, start, end) in &mut group.bins {
            *start = virtual_offset(*start)?;
            *end = virtual_offset(*end)?;
        }
        for offset in &mut group.linear {
            if *offset != 0 {
                *offset = virtual_offset(*offset)?;
            }
        }
    }
    write_index(&index, &groups, lines)?;
    let event = format!(
        "event=pairix_done records={records} chromosome_pairs={} elapsed_seconds={:.3} output={}",
        groups.len(),
        started.elapsed().as_secs_f64(),
        output.display()
    );
    eprintln!("{event}");
    if let Some(log) = &mut log {
        writeln!(log, "{event}")?;
    }
    Ok(())
}
