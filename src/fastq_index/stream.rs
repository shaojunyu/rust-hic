//! Paired FASTQ validation and aligner-specific streaming.
use super::{fingerprint, open_indexed_reader, rapidgzip, Manifest};
use anyhow::{bail, Context, Result};
use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    os::fd::AsRawFd,
    os::unix::fs::{FileExt, OpenOptionsExt},
    path::{Path, PathBuf},
    sync::atomic::{AtomicBool, Ordering},
    thread,
    time::Duration,
};

fn open_paired_readers(
    out: &Path,
    chunk: usize,
) -> Result<[io::Take<rapidgzip::RapidgzipReader>; 2]> {
    let m: Manifest = serde_json::from_reader(File::open(out.join("manifest.json"))?)?;
    if m.version != 1
        || m.inputs.len() != 2
        || m.boundaries.len() != m.chunks + 1
        || chunk >= m.chunks
    {
        bail!("invalid FASTQ index or chunk number")
    }
    let mut readers = Vec::with_capacity(2);
    for (mate, input) in m.inputs.iter().enumerate() {
        let fingerprint = fingerprint(
            Path::new(&input.path),
            input.index.clone(),
            input.uncompressed_size,
        )?;
        if (
            fingerprint.size,
            fingerprint.modified_ns,
            &fingerprint.edge_sha256,
        ) != (input.size, input.modified_ns, &input.edge_sha256)
        {
            bail!("input changed: {}", input.path)
        }
        let a = &m.boundaries[chunk];
        let b = &m.boundaries[chunk + 1];
        let (start, end) = if mate == 0 {
            (a.r1_uncompressed_offset, b.r1_uncompressed_offset)
        } else {
            (a.r2_uncompressed_offset, b.r2_uncompressed_offset)
        };
        let length = end
            .checked_sub(start)
            .context("invalid FASTQ chunk boundaries")?;
        let reader = open_indexed_reader(input, out, start, 1)?;
        readers.push(reader.take(length));
    }
    readers
        .try_into()
        .map_err(|_| anyhow::anyhow!("expected two FASTQ mates"))
}

struct PairedRecords<R: BufRead> {
    readers: [R; 2],
    records: [Vec<u8>; 2],
    pairs: u64,
}

impl<R: BufRead> PairedRecords<R> {
    fn new(readers: [R; 2]) -> Self {
        Self {
            readers,
            records: [Vec::new(), Vec::new()],
            pairs: 0,
        }
    }

    fn next(&mut self) -> Result<bool> {
        let first = read_fastq_for_interleaving(&mut self.readers[0], &mut self.records[0])?;
        let second = read_fastq_for_interleaving(&mut self.readers[1], &mut self.records[1])?;
        match (first, second) {
            (false, false) => Ok(false),
            (true, true) if fastq_id(&self.records[0]) == fastq_id(&self.records[1]) => {
                self.pairs = self
                    .pairs
                    .checked_add(1)
                    .context("FASTQ pair count overflow")?;
                Ok(true)
            }
            (true, true) => bail!("FASTQ mate IDs differ at pair {}", self.pairs + 1),
            _ => bail!("FASTQ mate counts differ at pair {}", self.pairs + 1),
        }
    }
}

pub(crate) struct PairedFastqReader {
    paired: PairedRecords<BufReader<io::Take<rapidgzip::RapidgzipReader>>>,
}

pub(crate) struct FastqPair {
    pub name: Vec<u8>,
    pub seq1: Vec<u8>,
    pub qual1: Vec<u8>,
    pub seq2: Vec<u8>,
    pub qual2: Vec<u8>,
}

impl PairedFastqReader {
    pub fn open(index: &Path, chunk: usize) -> Result<Self> {
        Ok(Self {
            paired: PairedRecords::new(
                open_paired_readers(index, chunk)?.map(|r| BufReader::with_capacity(128 * 1024, r)),
            ),
        })
    }

    pub fn next_pair(&mut self) -> Result<Option<FastqPair>> {
        if !self.paired.next()? {
            if self.paired.readers.iter().any(|r| r.get_ref().limit() != 0) {
                bail!("truncated FASTQ chunk")
            }
            return Ok(None);
        }
        fn parts(record: &[u8]) -> Result<(Vec<u8>, Vec<u8>)> {
            let mut lines = record.split(|&b| b == b'\n');
            let _ = lines.next();
            let seq = lines.next().context("missing FASTQ sequence")?;
            let _ = lines.next();
            let qual = lines.next().context("missing FASTQ quality")?;
            Ok((
                seq.strip_suffix(b"\r").unwrap_or(seq).to_vec(),
                qual.strip_suffix(b"\r").unwrap_or(qual).to_vec(),
            ))
        }
        let (seq1, qual1) = parts(&self.paired.records[0])?;
        let (seq2, qual2) = parts(&self.paired.records[1])?;
        Ok(Some(FastqPair {
            name: fastq_id(&self.paired.records[0])[1..].to_vec(),
            seq1,
            qual1,
            seq2,
            qual2,
        }))
    }
}

/// Stream the two mates as interleaved FASTQ for BWA's paired mode.
pub(crate) fn stream_paired<W: Write>(
    out: &Path,
    chunk: usize,
    output: W,
    cancel: &AtomicBool,
) -> Result<u64> {
    let mut paired = PairedRecords::new(
        open_paired_readers(out, chunk)?.map(|reader| BufReader::with_capacity(128 * 1024, reader)),
    );
    let mut writer = BufWriter::with_capacity(128 * 1024, output);
    loop {
        if cancel.load(Ordering::Relaxed) {
            bail!("FASTQ stream cancelled")
        }
        if !paired.next().with_context(|| format!("chunk {chunk}"))? {
            if paired.readers.iter().any(|r| r.get_ref().limit() != 0) {
                bail!("truncated FASTQ chunk")
            }
            writer.flush()?;
            return Ok(paired.pairs);
        }
        writer.write_all(&paired.records[0])?;
        writer.write_all(&paired.records[1])?;
    }
}

/// Feed validated mates to two named pipes. A lagging consumer may require
/// disk-backed spill: Chromap can read a large batch from one FIFO before the
/// other, so a fixed-size in-memory queue would deadlock.
pub(crate) fn stream_split(
    out: &Path,
    chunk: usize,
    paths: [&Path; 2],
    cancel: &AtomicBool,
) -> Result<[u64; 2]> {
    const BUFFER_SIZE: usize = 128 * 1024;
    struct Spill {
        writer: BufWriter<File>,
        read: u64,
        written: u64,
    }
    struct Mate {
        path: PathBuf,
        writer: Option<File>,
        buffer: Vec<u8>,
        start: usize,
        spill: Option<Spill>,
        written: u64,
        done: bool,
    }
    impl Mate {
        fn append(&mut self, record: &[u8], scratch: &Path) -> Result<()> {
            if self.start == self.buffer.len() && self.spill.is_none() {
                self.buffer.clear();
                self.start = 0;
            }
            if self.spill.is_none() && self.buffer.len() - self.start < BUFFER_SIZE {
                self.buffer.extend_from_slice(record);
            } else {
                if self.spill.is_none() {
                    let file = tempfile::tempfile_in(scratch).context("create FASTQ spill file")?;
                    self.spill = Some(Spill {
                        writer: BufWriter::with_capacity(BUFFER_SIZE, file),
                        read: 0,
                        written: 0,
                    });
                }
                let spill = self.spill.as_mut().expect("spill created");
                spill
                    .writer
                    .write_all(record)
                    .context("write FASTQ spill file")?;
                spill.written = spill
                    .written
                    .checked_add(record.len() as u64)
                    .context("FASTQ spill size overflow")?;
            }
            Ok(())
        }
        fn load_spill(&mut self) -> Result<()> {
            if self.start != self.buffer.len() {
                return Ok(());
            }
            let Some(spill) = &mut self.spill else {
                return Ok(());
            };
            spill.writer.flush().context("flush FASTQ spill file")?;
            let n = ((spill.written - spill.read).min(BUFFER_SIZE as u64)) as usize;
            self.buffer.resize(n, 0);
            let mut got = 0;
            while got < n {
                let count = spill
                    .writer
                    .get_ref()
                    .read_at(&mut self.buffer[got..], spill.read + got as u64)
                    .context("read FASTQ spill file")?;
                if count == 0 {
                    bail!("truncated FASTQ spill file")
                }
                got += count;
            }
            self.start = 0;
            spill.read += n as u64;
            if spill.read == spill.written {
                self.spill = None;
            }
            Ok(())
        }
        fn needs_input(&self) -> bool {
            self.start == self.buffer.len() && self.spill.is_none()
        }
    }
    let mut paired = PairedRecords::new(
        open_paired_readers(out, chunk)?
            .map(|reader| BufReader::with_capacity(BUFFER_SIZE, reader)),
    );
    let make = |path: &Path| Mate {
        path: path.to_owned(),
        writer: None,
        buffer: Vec::with_capacity(BUFFER_SIZE),
        start: 0,
        spill: None,
        written: 0,
        done: false,
    };
    let mut mates = [make(paths[0]), make(paths[1])];
    let scratch = paths[0]
        .parent()
        .context("FASTQ pipe has no parent directory")?;
    let mut eof = false;
    loop {
        if cancel.load(Ordering::Relaxed) {
            bail!("FASTQ stream cancelled")
        }
        let mut progressed = false;
        for mate in &mut mates {
            mate.load_spill()?;
        }
        if !eof && mates.iter().any(Mate::needs_input) {
            // Fill both mates together while one side needs more input; the
            // slower side spills to disk once its memory buffer is full.
            let needy = [mates[0].needs_input(), mates[1].needs_input()];
            while (0..2).any(|i| needy[i] && mates[i].buffer.len() - mates[i].start < BUFFER_SIZE) {
                if cancel.load(Ordering::Relaxed) {
                    bail!("FASTQ stream cancelled")
                }
                if !paired.next().with_context(|| format!("chunk {chunk}"))? {
                    if paired.readers.iter().any(|r| r.get_ref().limit() != 0) {
                        bail!("truncated FASTQ chunk")
                    }
                    eof = true;
                    break;
                }
                for (mate, record) in mates.iter_mut().zip(&paired.records) {
                    mate.append(record, scratch)?;
                }
            }
            progressed = true;
        }
        for mate in &mut mates {
            if mate.done {
                continue;
            }
            if mate.writer.is_none() {
                match std::fs::OpenOptions::new()
                    .write(true)
                    .custom_flags(libc::O_NONBLOCK)
                    .open(&mate.path)
                {
                    Ok(file) => {
                        mate.writer = Some(file);
                        progressed = true;
                    }
                    Err(error) if error.raw_os_error() == Some(libc::ENXIO) => continue,
                    Err(error) => {
                        return Err(error).with_context(|| format!("open {}", mate.path.display()))
                    }
                }
            }
            if mate.start == mate.buffer.len() {
                if eof && mate.spill.is_none() {
                    mate.writer.take();
                    mate.done = true;
                    progressed = true;
                }
                continue;
            }
            match mate
                .writer
                .as_mut()
                .expect("writer opened")
                .write(&mate.buffer[mate.start..])
            {
                Ok(0) => bail!("zero-length write to FASTQ pipe"),
                Ok(n) => {
                    mate.start += n;
                    mate.written = mate
                        .written
                        .checked_add(n as u64)
                        .context("FASTQ byte count overflow")?;
                    progressed = true;
                }
                Err(error)
                    if matches!(
                        error.kind(),
                        io::ErrorKind::WouldBlock | io::ErrorKind::Interrupted
                    ) => {}
                Err(error) => return Err(error).context("write FASTQ pipe"),
            }
        }
        if mates.iter().all(|mate| mate.done) {
            return Ok([mates[0].written, mates[1].written]);
        }
        if !progressed {
            let mut fds: Vec<_> = mates
                .iter()
                .filter(|mate| mate.start < mate.buffer.len())
                .filter_map(|mate| mate.writer.as_ref())
                .map(|file| libc::pollfd {
                    fd: file.as_raw_fd(),
                    events: libc::POLLOUT,
                    revents: 0,
                })
                .collect();
            if fds.is_empty() {
                thread::sleep(Duration::from_millis(25));
            } else {
                let rc = unsafe { libc::poll(fds.as_mut_ptr(), fds.len() as libc::nfds_t, 25) };
                if rc < 0 {
                    let error = io::Error::last_os_error();
                    if error.kind() != io::ErrorKind::Interrupted {
                        return Err(error).context("poll FASTQ pipes");
                    }
                }
            }
        }
    }
}

fn read_fastq_for_interleaving<R: BufRead>(reader: &mut R, record: &mut Vec<u8>) -> Result<bool> {
    record.clear();
    let mut offsets = [0usize; 5];
    for i in 0..4 {
        let n = reader.read_until(b'\n', record)?;
        if n == 0 {
            if i == 0 {
                return Ok(false);
            }
            bail!("truncated FASTQ record")
        }
        offsets[i + 1] = record.len();
    }
    let line = |i: usize| {
        let value = &record[offsets[i]..offsets[i + 1]];
        let value = value.strip_suffix(b"\n").unwrap_or(value);
        value.strip_suffix(b"\r").unwrap_or(value)
    };
    if !line(0).starts_with(b"@") || !line(2).starts_with(b"+") || line(1).len() != line(3).len() {
        bail!("invalid FASTQ record")
    }
    if !record.ends_with(b"\n") {
        record.push(b'\n');
    }
    Ok(true)
}

fn fastq_id(record: &[u8]) -> &[u8] {
    let id = record
        .split(|b| b.is_ascii_whitespace())
        .next()
        .unwrap_or_default();
    id.strip_suffix(b"/1")
        .or_else(|| id.strip_suffix(b"/2"))
        .unwrap_or(id)
}

#[cfg(test)]
mod tests {
    use super::{fastq_id, read_fastq_for_interleaving, PairedRecords};
    use std::io::Cursor;

    #[test]
    fn interleaved_records_normalize_mate_ids_and_terminate_lines() {
        let mut first = Cursor::new(b"@read/1 extra\r\nACGT\r\n+\r\nIIII\r\n".as_slice());
        let mut second = Cursor::new(b"@read/2\nACGT\n+\nIIII".as_slice());
        let (mut a, mut b) = (Vec::new(), Vec::new());
        assert!(read_fastq_for_interleaving(&mut first, &mut a).unwrap());
        assert!(read_fastq_for_interleaving(&mut second, &mut b).unwrap());
        assert_eq!(fastq_id(&a), fastq_id(&b));
        assert!(b.ends_with(b"IIII\n"));
        assert!(!read_fastq_for_interleaving(&mut first, &mut a).unwrap());
        assert!(!read_fastq_for_interleaving(&mut second, &mut b).unwrap());
    }

    #[test]
    fn interleaved_records_reject_truncation_and_bad_quality() {
        for input in [
            b"@read\nACGT\n+\n".as_slice(),
            b"@read\nACGT\n+\nIII\n".as_slice(),
        ] {
            assert!(read_fastq_for_interleaving(&mut Cursor::new(input), &mut Vec::new()).is_err());
        }
    }

    #[test]
    fn paired_records_reject_internal_mismatches_and_missing_mates() {
        let first = b"@first/1\nAC\n+\nII\n@second/1\nGT\n+\nII\n";
        for second in [
            b"@first/2\nAC\n+\nII\n@wrong/2\nGT\n+\nII\n".as_slice(),
            b"@first/2\nAC\n+\nII\n".as_slice(),
        ] {
            let mut pairs =
                PairedRecords::new([Cursor::new(first.as_slice()), Cursor::new(second)]);
            assert!(pairs.next().unwrap());
            assert!(pairs.next().is_err());
        }
    }
}
