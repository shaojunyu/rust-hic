mod build;
mod cli;
mod manifest;
use build::{build_index, compact_existing};
use manifest::{fingerprint, Boundary, Checkpoint, Input, Manifest};
mod rapidgzip;
mod stream;
pub(crate) use cli::main;
pub(crate) use stream::{stream_paired, stream_split, FastqPair, PairedFastqReader};
mod sparse;
use anyhow::{bail, Context, Result};
use std::{
    env,
    fs::{self, File},
    io::{self, BufRead, BufReader, Read, Seek, SeekFrom},
    path::Path,
    time::Instant,
};

fn read_record<R: BufRead>(r: &mut R, offset: &mut u64, label: &str) -> Result<Option<Vec<u8>>> {
    let mut lines = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    for (i, line) in lines.iter_mut().enumerate() {
        let n = r.read_until(b'\n', line)?;
        *offset += n as u64;
        if n == 0 {
            if i == 0 {
                return Ok(None);
            }
            bail!("truncated {label} FASTQ record")
        }
    }
    if !lines[0].starts_with(b"@") || !lines[2].starts_with(b"+") {
        bail!("invalid {label} FASTQ at uncompressed offset {}", *offset)
    }
    fn trim(x: &[u8]) -> usize {
        x.strip_suffix(b"\n")
            .unwrap_or(x)
            .strip_suffix(b"\r")
            .unwrap_or(x.strip_suffix(b"\n").unwrap_or(x))
            .len()
    }
    if trim(&lines[1]) != trim(&lines[3]) {
        bail!("sequence/quality length mismatch in {label}")
    }
    let id = lines[0].split(|b| b.is_ascii_whitespace()).next().unwrap();
    let id = id
        .strip_suffix(b"/1")
        .or_else(|| id.strip_suffix(b"/2"))
        .unwrap_or(id);
    Ok(Some(id.to_vec()))
}
fn compact_index(path: &Path, offsets: &[u64]) -> Result<()> {
    use std::io::Write;
    let mut input = File::open(path)?;
    let mut head = [0u8; 35];
    input.read_exact(&mut head)?;
    if &head[..5] != b"GZIDX" || head[5] != 1 {
        bail!("unsupported gzip index format")
    }
    let count = u32::from_le_bytes(head[31..35].try_into()?) as usize;
    let window = u32::from_le_bytes(head[27..31].try_into()?) as usize;
    if count == 0 || count > 10_000_000 || window > 1_048_576 {
        bail!("invalid gzip index header")
    }
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        let mut row = [0u8; 18];
        input.read_exact(&mut row)?;
        rows.push(row);
    }
    let positions: Vec<u64> = rows
        .iter()
        .map(|r| u64::from_le_bytes(r[8..16].try_into().unwrap()))
        .collect();
    let selected = sparse::select_points(&positions, offsets)?;
    let selected_count = selected.iter().filter(|x| **x).count();
    head[31..35].copy_from_slice(&(selected_count as u32).to_le_bytes());
    let tmp = path.with_extension("gzidx.compacting");
    let mut output = File::create(&tmp)?;
    output.write_all(&head)?;
    for (row, keep) in rows.iter().zip(&selected) {
        if *keep {
            output.write_all(row)?;
        }
    }
    let mut data = vec![0u8; window];
    for (row, keep) in rows.iter().zip(&selected) {
        if row[17] != 0 {
            input.read_exact(&mut data)?;
            if *keep {
                output.write_all(&data)?;
            }
        }
    }
    drop(output);
    fs::rename(tmp, path)?;
    eprintln!(
        "compacted gzip index {}: {} -> {} checkpoints",
        path.display(),
        count,
        selected_count
    );
    Ok(())
}
fn checkpoints(path: &Path) -> Result<Vec<Checkpoint>> {
    let mut f = File::open(path)?;
    let mut header = [0u8; 35];
    f.read_exact(&mut header)?;
    if &header[..5] != b"GZIDX" || header[5] != 1 {
        bail!("unsupported gzip index format")
    }
    let count = u32::from_le_bytes(header[31..35].try_into()?) as usize;
    if count > 10_000_000 {
        bail!("invalid index checkpoint count")
    }
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        let mut row = [0u8; 18];
        f.read_exact(&mut row)?;
        out.push(Checkpoint {
            compressed_byte_offset: u64::from_le_bytes(row[0..8].try_into()?),
            uncompressed_offset: u64::from_le_bytes(row[8..16].try_into()?),
            compressed_bits_to_prime: row[16],
            has_window: row[17] != 0,
        });
    }
    if out
        .windows(2)
        .any(|w| w[0].uncompressed_offset > w[1].uncompressed_offset)
    {
        bail!("unsorted gzip checkpoints")
    }
    Ok(out)
}
fn nearest(points: &[Checkpoint], offset: u64) -> Option<Checkpoint> {
    points
        .partition_point(|p| p.uncompressed_offset <= offset)
        .checked_sub(1)
        .map(|i| points[i].clone())
}
fn validate_index_density(path: &Path, size: u64) -> Result<()> {
    let points = checkpoints(path)?;
    if points.first().map(|p| p.uncompressed_offset) != Some(0) {
        bail!("index has no initial checkpoint; rebuild the FASTQ index");
    }
    let mut previous = 0;
    for offset in points
        .iter()
        .map(|p| p.uncompressed_offset)
        .chain(std::iter::once(size))
    {
        let gap = offset
            .checked_sub(previous)
            .context("invalid index checkpoint order")?;
        if gap > sparse::MAX_CHECKPOINT_GAP {
            bail!("FASTQ index checkpoints are too sparse for rapidgzip (gap {} MiB); rebuild the index", gap / (1024 * 1024));
        }
        previous = offset;
    }
    Ok(())
}
fn open_indexed_reader(
    input: &Input,
    out: &Path,
    start: u64,
    threads: u32,
) -> Result<rapidgzip::RapidgzipReader> {
    let path = out.join(&input.index);
    validate_index_density(&path, input.uncompressed_size)?;
    let mut reader = rapidgzip::RapidgzipReader::open(Path::new(&input.path), threads)?;
    reader.import_index(&path)?;
    reader.seek(SeekFrom::Start(start))?;
    Ok(reader)
}
pub(crate) fn verify(out: &Path, _threads: u32) -> Result<()> {
    let manifest: Manifest = serde_json::from_reader(File::open(out.join("manifest.json"))?)?;
    if manifest.inputs.len() != 2 || manifest.boundaries.len() != manifest.chunks + 1 {
        bail!("invalid manifest")
    }
    for input in &manifest.inputs {
        validate_index_density(&out.join(&input.index), input.uncompressed_size)?;
        let f = fingerprint(
            Path::new(&input.path),
            input.index.clone(),
            input.uncompressed_size,
        )?;
        if (f.size, f.modified_ns, &f.edge_sha256)
            != (input.size, input.modified_ns, &input.edge_sha256)
        {
            bail!("input changed: {}", input.path)
        }
    }
    let mut readers = manifest
        .inputs
        .iter()
        .map(|input| {
            let mut reader = rapidgzip::RapidgzipReader::open(Path::new(&input.path), 1)?;
            reader.import_index(&out.join(&input.index))?;
            Ok(reader)
        })
        .collect::<Result<Vec<_>>>()?;
    for b in manifest.boundaries.iter().take(manifest.chunks) {
        for (i, offset) in [b.r1_uncompressed_offset, b.r2_uncompressed_offset]
            .into_iter()
            .enumerate()
        {
            readers[i].seek(SeekFrom::Start(offset))?;
            let mut buf = BufReader::new(&mut readers[i]);
            let mut pos = offset;
            let id = read_record(&mut buf, &mut pos, "indexed")?
                .context("unexpected EOF at boundary")?;
            if id != b.read_id.as_bytes() {
                bail!(
                    "boundary mismatch at pair {} for R{}",
                    b.pair
                        .map_or_else(|| "unknown".to_string(), |x| x.to_string()),
                    i + 1
                )
            }
        }
        eprintln!(
            "verified pair={} R1={} R2={} id={}",
            b.pair
                .map_or_else(|| "unknown".to_string(), |x| x.to_string()),
            b.r1_uncompressed_offset,
            b.r2_uncompressed_offset,
            b.read_id
        );
    }
    Ok(())
}
fn extract(out: &Path, chunk: usize, destinations: [&Path; 2], _threads: u32) -> Result<()> {
    let manifest: Manifest = serde_json::from_reader(File::open(out.join("manifest.json"))?)?;
    if manifest.version != 1
        || manifest.inputs.len() != 2
        || manifest.boundaries.len() != manifest.chunks + 1
        || chunk >= manifest.chunks
    {
        bail!("invalid manifest or chunk number")
    }
    for (i, dest) in destinations.into_iter().enumerate() {
        let input = &manifest.inputs[i];
        let fp = fingerprint(
            Path::new(&input.path),
            input.index.clone(),
            input.uncompressed_size,
        )?;
        if (fp.size, fp.modified_ns, &fp.edge_sha256)
            != (input.size, input.modified_ns, &input.edge_sha256)
        {
            bail!("input changed: {}", input.path)
        }
        let a = &manifest.boundaries[chunk];
        let b = &manifest.boundaries[chunk + 1];
        let (start, end) = if i == 0 {
            (a.r1_uncompressed_offset, b.r1_uncompressed_offset)
        } else {
            (a.r2_uncompressed_offset, b.r2_uncompressed_offset)
        };
        let count = end.checked_sub(start).context("invalid boundary order")?;
        let r = open_indexed_reader(input, out, start, 1)?;
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(dest)?;
        let copied = io::copy(&mut r.take(count), &mut file)?;
        if copied != count {
            bail!("truncated data for R{}", i + 1)
        }
        eprintln!(
            "chunk={} R{} start={} end={} bytes={}",
            chunk,
            i + 1,
            start,
            end,
            copied
        );
    }
    Ok(())
}
fn stream_chunk(out: &Path, chunk: usize, mate: usize) -> Result<()> {
    let m: Manifest = serde_json::from_reader(File::open(out.join("manifest.json"))?)?;
    if m.version != 1
        || m.inputs.len() != 2
        || m.boundaries.len() != m.chunks + 1
        || chunk >= m.chunks
        || !(1..=2).contains(&mate)
    {
        bail!("invalid manifest, chunk or mate")
    }
    let i = &m.inputs[mate - 1];
    let f = fingerprint(Path::new(&i.path), i.index.clone(), i.uncompressed_size)?;
    if (f.size, f.modified_ns, &f.edge_sha256) != (i.size, i.modified_ns, &i.edge_sha256) {
        bail!("input changed: {}", i.path)
    }
    let a = &m.boundaries[chunk];
    let b = &m.boundaries[chunk + 1];
    let (start, end) = if mate == 1 {
        (a.r1_uncompressed_offset, b.r1_uncompressed_offset)
    } else {
        (a.r2_uncompressed_offset, b.r2_uncompressed_offset)
    };
    let count = end.checked_sub(start).context("invalid boundary order")?;
    let reader = open_indexed_reader(i, out, start, 1)?;
    let copied = io::copy(&mut reader.take(count), &mut io::stdout().lock())?;
    if copied != count {
        bail!("short stream for R{mate}: expected {count}, got {copied}")
    }
    eprintln!("chunk={chunk} mate={mate} start={start} end={end} bytes={copied}");
    Ok(())
}
