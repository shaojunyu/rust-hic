use super::attr::{int_attr, str_attr as attr_str};
use crate::pairs::format::{read_header, Format};
use anyhow::{bail, Context, Result};
use hdf5::{
    types::{EnumMember, EnumType, FixedAscii, IntSize, TypeDescriptor},
    File, Group,
};
use rayon::prelude::*;
use std::{
    cmp::Reverse,
    collections::{BinaryHeap, HashMap},
    env,
    fs::{self, File as StdFile},
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    sync::mpsc::{sync_channel, SyncSender},
    thread,
    time::Instant,
};
#[derive(Clone)]
struct Chrom {
    name: String,
    len: u64,
    offset: u64,
}
const READ_BATCH: usize = 8 * 1024 * 1024;
fn produce_batches(path: &Path, tx: SyncSender<(Format, Vec<u8>)>) -> Result<()> {
    let mut input = BufReader::with_capacity(1024 * 1024, StdFile::open(path)?);
    let (_, format) = read_header(&mut input)?;
    let mut pending = Vec::with_capacity(READ_BATCH + 1024 * 1024);
    let mut buffer = vec![0u8; 1024 * 1024];
    loop {
        let n = input.read(&mut buffer)?;
        if n == 0 {
            if !pending.is_empty() {
                tx.send((format, pending))
                    .context("binning reader canceled")?;
            }
            break;
        }
        pending.extend_from_slice(&buffer[..n]);
        if pending.len() >= READ_BATCH {
            if let Some(last) = pending.iter().rposition(|b| *b == b'\n') {
                let rest = pending.split_off(last + 1);
                tx.send((format, pending))
                    .context("binning reader canceled")?;
                pending = rest;
            }
        }
    }
    Ok(())
}
fn parse_piece(
    data: &[u8],
    format: Format,
    lookup: &HashMap<String, usize>,
    chroms: &[Chrom],
    bin: u64,
) -> Result<Vec<u64>> {
    let text = std::str::from_utf8(data).context("pairs are not UTF-8")?;
    let mut rows = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            bail!("header inside body")
        }
        let f = format.fields(line)?;
        let c1 = f[1];
        let p1: u64 = f[2].parse()?;
        let c2 = f[3];
        let p2: u64 = f[4].parse()?;
        let a = &chroms[*lookup
            .get(c1)
            .with_context(|| format!("unknown chromosome {c1}"))?];
        let b = &chroms[*lookup
            .get(c2)
            .with_context(|| format!("unknown chromosome {c2}"))?];
        if p1 == 0 || p2 == 0 || p1 > a.len || p2 > b.len {
            bail!("coordinate outside chromosome")
        }
        let mut x = a.offset + (p1 - 1) / bin;
        let mut y = b.offset + (p2 - 1) / bin;
        if x > y {
            std::mem::swap(&mut x, &mut y)
        }
        rows.push((x << 32) | y);
    }
    Ok(rows)
}
fn parse_batch(
    data: &[u8],
    format: Format,
    workers: usize,
    lookup: &HashMap<String, usize>,
    chroms: &[Chrom],
    bin: u64,
) -> Result<Vec<Vec<u64>>> {
    let mut boundaries = vec![0usize];
    for i in 1..workers {
        let target = data.len() * i / workers;
        if let Some(delta) = data[target..].iter().position(|b| *b == b'\n') {
            let cut = target + delta + 1;
            if cut > *boundaries.last().unwrap() && cut < data.len() {
                boundaries.push(cut);
            }
        }
    }
    boundaries.push(data.len());
    let pieces: Vec<_> = boundaries.windows(2).map(|w| &data[w[0]..w[1]]).collect();
    pieces
        .into_par_iter()
        .map(|piece| parse_piece(piece, format, lookup, chroms, bin))
        .collect()
}
fn attr_u64(g: &Group, k: &str, v: u64) -> Result<()> {
    let v = i64::try_from(v).context("attribute exceeds i64")?;
    int_attr(g, k, v)
}
fn run_file(rows: &mut Vec<u64>, dir: &Path, n: usize) -> Result<PathBuf> {
    rows.par_sort_unstable();
    let path = dir.join(format!("run-{n:06}.bin"));
    let mut out = BufWriter::with_capacity(1024 * 1024, StdFile::create(&path)?);
    let mut i = 0;
    while i < rows.len() {
        let key = rows[i];
        let mut j = i + 1;
        while j < rows.len() && rows[j] == key {
            j += 1
        }
        out.write_all(&key.to_le_bytes())?;
        out.write_all(&((j - i) as u64).to_le_bytes())?;
        i = j
    }
    out.flush()?;
    rows.clear();
    Ok(path)
}
fn read_rec(r: &mut BufReader<StdFile>) -> Result<Option<(u64, u64)>> {
    let mut b = [0u8; 16];
    let mut n = 0;
    while n < 16 {
        let got = r.read(&mut b[n..])?;
        if got == 0 {
            if n == 0 {
                return Ok(None);
            }
            bail!("truncated run file")
        }
        n += got
    }
    Ok(Some((
        u64::from_le_bytes(b[..8].try_into()?),
        u64::from_le_bytes(b[8..].try_into()?),
    )))
}
struct PixelWriter {
    a: hdf5::Dataset,
    b: hdf5::Dataset,
    c: hdf5::Dataset,
    aa: Vec<i64>,
    bb: Vec<i64>,
    cc: Vec<i32>,
    written: usize,
}
impl PixelWriter {
    fn flush(&mut self) -> Result<()> {
        if self.aa.is_empty() {
            return Ok(());
        }
        let n = self.aa.len();
        let end = self.written + n;
        for d in [&self.a, &self.b, &self.c] {
            d.resize((end,))?
        }
        self.a.write_slice(&self.aa[..], self.written..end)?;
        self.b.write_slice(&self.bb[..], self.written..end)?;
        self.c.write_slice(&self.cc[..], self.written..end)?;
        self.written = end;
        self.aa.clear();
        self.bb.clear();
        self.cc.clear();
        Ok(())
    }
}
fn run() -> Result<()> {
    let args: Vec<_> = env::args().collect();
    if args.len() != 8 {
        bail!(
            "usage: bin_pairs CHROMSIZES BIN_SIZE ASSEMBLY INPUT.pairs OUTPUT.cool TMPDIR THREADS"
        )
    }
    let begin = Instant::now();
    let bin: u64 = args[2].parse()?;
    if bin == 0 || bin > i32::MAX as u64 {
        bail!("invalid bin size")
    };
    let threads: usize = args[7].parse()?;
    if threads == 0 || threads > 128 {
        bail!("invalid thread count")
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;
    let output = Path::new(&args[5]);
    if output.exists() {
        bail!("output exists")
    };
    let work = Path::new(&args[6]).join(format!("pairs-to-cool-{}", std::process::id()));
    fs::create_dir(&work).context("create scratch directory")?;
    let temp = output.with_extension(format!("cool.partial.{}", std::process::id()));
    if temp.exists() {
        bail!("temporary output already exists")
    }
    let mut chroms = Vec::new();
    let mut lookup = HashMap::new();
    let mut offset = 0u64;
    for l in BufReader::new(StdFile::open(&args[1])?).lines() {
        let l = l?;
        let f: Vec<_> = l.split_whitespace().collect();
        if f.is_empty() {
            continue;
        }
        if f.len() != 2 {
            bail!("invalid chromosome sizes")
        };
        let len: u64 = f[1].parse()?;
        if len == 0 || len > i32::MAX as u64 || lookup.contains_key(f[0]) {
            bail!("invalid chromosome")
        };
        lookup.insert(f[0].to_owned(), chroms.len());
        chroms.push(Chrom {
            name: f[0].into(),
            len,
            offset,
        });
        offset = offset
            .checked_add(len.div_ceil(bin))
            .context("bin overflow")?;
    }
    if chroms.is_empty() || offset > u32::MAX as u64 {
        bail!("chromosome table empty or exceeds u32 bin IDs")
    }
    let mut rows = Vec::with_capacity(3_000_000);
    let mut runs = Vec::new();
    let mut total = 0u64;
    let mut next_progress = 10_000_000u64;
    thread::scope(|scope| -> Result<()> {
        let (tx, rx) = sync_channel(2);
        let producer = scope.spawn(|| produce_batches(Path::new(&args[4]), tx));
        let processing = (|| -> Result<()> {
            while let Ok((format, batch)) = rx.recv() {
                for piece in parse_batch(&batch, format, threads, &lookup, &chroms, bin)? {
                    total = total
                        .checked_add(piece.len() as u64)
                        .context("pair count overflow")?;
                    rows.extend(piece);
                    if rows.len() >= 3_000_000 {
                        runs.push(run_file(&mut rows, &work, runs.len())?);
                    }
                }
                if total >= next_progress {
                    eprintln!(
                        "event=bin_progress pairs={total} runs={} elapsed_seconds={:.3}",
                        runs.len(),
                        begin.elapsed().as_secs_f64()
                    );
                    next_progress = (total / 10_000_000 + 1) * 10_000_000;
                }
            }
            Ok(())
        })();
        drop(rx);
        let reading = producer
            .join()
            .map_err(|_| anyhow::anyhow!("binning reader panicked"))?;
        processing?;
        reading?;
        Ok(())
    })?;
    if !rows.is_empty() {
        runs.push(run_file(&mut rows, &work, runs.len())?)
    }
    eprintln!(
        "event=bin_done pairs={total} runs={} elapsed_seconds={:.3}",
        runs.len(),
        begin.elapsed().as_secs_f64()
    );
    let file = File::create(&temp)?;
    let root = file.group("/")?;
    let cg = root.create_group("chroms")?;
    let bg = root.create_group("bins")?;
    let pg = root.create_group("pixels")?;
    let ig = root.create_group("indexes")?;
    let names: Vec<FixedAscii<256>> = chroms
        .iter()
        .map(|c| FixedAscii::<256>::from_ascii(&c.name))
        .collect::<std::result::Result<_, _>>()?;
    let lens: Vec<i32> = chroms.iter().map(|c| c.len as i32).collect();
    cg.new_dataset_builder().with_data(&names).create("name")?;
    cg.new_dataset_builder().with_data(&lens).create("length")?;
    let mut bc = Vec::with_capacity(offset as usize);
    let mut bs = Vec::with_capacity(offset as usize);
    let mut be = Vec::with_capacity(offset as usize);
    let mut co = Vec::with_capacity(chroms.len() + 1);
    for (i, c) in chroms.iter().enumerate() {
        co.push(c.offset as i64);
        let mut p = 0;
        while p < c.len {
            bc.push(i as i32);
            bs.push(p as i32);
            be.push((p + bin).min(c.len) as i32);
            p += bin
        }
    }
    co.push(offset as i64);
    let chrom_type = TypeDescriptor::Enum(EnumType {
        size: IntSize::U4,
        signed: true,
        members: chroms
            .iter()
            .enumerate()
            .map(|(i, c)| EnumMember {
                name: c.name.clone(),
                value: i as u64,
            })
            .collect(),
    });
    let chrom_ds = bg
        .new_dataset_builder()
        .empty_as(&chrom_type)
        .shape(bc.len())
        .create("chrom")?;
    let chrom_dtype = chrom_ds.dtype()?;
    // The in-memory values use the same native i32 representation as the enum base type.
    // HDF5 does not provide an integer-to-enum conversion, so write with the enum type.
    let status = hdf5::sync::sync(|| unsafe {
        hdf5_sys::h5d::H5Dwrite(
            chrom_ds.id(),
            chrom_dtype.id(),
            hdf5_sys::h5s::H5S_ALL,
            hdf5_sys::h5s::H5S_ALL,
            hdf5_sys::h5p::H5P_DEFAULT,
            bc.as_ptr().cast(),
        )
    });
    if status < 0 {
        bail!("failed to write chromosome enum values")
    }
    bg.new_dataset_builder().with_data(&bs).create("start")?;
    bg.new_dataset_builder().with_data(&be).create("end")?;
    ig.new_dataset_builder()
        .with_data(&co)
        .create("chrom_offset")?;
    let a = pg
        .new_dataset::<i64>()
        .shape((0..,))
        .chunk((65_536,))
        .create("bin1_id")?;
    let b = pg
        .new_dataset::<i64>()
        .shape((0..,))
        .chunk((65_536,))
        .create("bin2_id")?;
    let c = pg
        .new_dataset::<i32>()
        .shape((0..,))
        .chunk((65_536,))
        .create("count")?;
    let mut writer = PixelWriter {
        a,
        b,
        c,
        aa: Vec::with_capacity(65536),
        bb: Vec::with_capacity(65536),
        cc: Vec::with_capacity(65536),
        written: 0,
    };
    let mut offsets = vec![0i64; offset as usize + 1];
    let mut readers: Vec<_> = runs
        .iter()
        .map(|p| StdFile::open(p).map(|f| BufReader::with_capacity(64 * 1024, f)))
        .collect::<io::Result<_>>()?;
    let mut heap = BinaryHeap::new();
    for (i, r) in readers.iter_mut().enumerate() {
        if let Some((key, count)) = read_rec(r)? {
            heap.push(Reverse((key, i, count)))
        }
    }
    let (mut nnz, mut next_bin, mut sum) = (0u64, 0usize, 0u64);
    while let Some(Reverse((key, i, count))) = heap.pop() {
        let mut value = count;
        if let Some((k, n)) = read_rec(&mut readers[i])? {
            heap.push(Reverse((k, i, n)))
        }
        while heap.peek().map(|x| x.0 .0) == Some(key) {
            let Reverse((_, j, n)) = heap.pop().unwrap();
            value = value.checked_add(n).context("count overflow")?;
            if let Some((k, n)) = read_rec(&mut readers[j])? {
                heap.push(Reverse((k, j, n)))
            }
        }
        let x = (key >> 32) as usize;
        let y = (key & 0xffff_ffff) as i64;
        while next_bin <= x {
            offsets[next_bin] = nnz as i64;
            next_bin += 1
        }
        writer.aa.push(x as i64);
        writer.bb.push(y);
        writer
            .cc
            .push(i32::try_from(value).context("pixel count exceeds i32")?);
        sum += value;
        nnz += 1;
        if writer.aa.len() == 65536 {
            writer.flush()?
        }
    }
    while next_bin < offsets.len() {
        offsets[next_bin] = nnz as i64;
        next_bin += 1
    }
    writer.flush()?;
    if sum != total {
        bail!("count sum mismatch")
    };
    ig.new_dataset_builder()
        .with_data(&offsets)
        .create("bin1_offset")?;
    for (k, v) in [
        ("format", "HDF5::Cooler"),
        ("bin-type", "fixed"),
        ("storage-mode", "symmetric-upper"),
        ("genome-assembly", args[3].as_str()),
        ("generated-by", "rust-hic Rust"),
        ("metadata", "{}"),
    ] {
        attr_str(&root, k, v)?
    }
    for (k, v) in [
        ("format-version", 3),
        ("bin-size", bin),
        ("nchroms", chroms.len() as u64),
        ("nbins", offset),
        ("nnz", nnz),
        ("sum", sum),
    ] {
        attr_u64(&root, k, v)?
    }
    file.flush()?;
    drop(writer);
    drop(file);
    fs::hard_link(&temp, output).context("publish cooler without overwriting")?;
    fs::remove_file(&temp)?;
    for p in runs {
        fs::remove_file(p)?
    }
    fs::remove_dir(&work)?;
    eprintln!(
        "event=cool_done input_records={total} nnz={nnz} elapsed_seconds={:.3}",
        begin.elapsed().as_secs_f64()
    );
    Ok(())
}
pub(crate) fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e:#}");
        std::process::exit(1)
    }
}
