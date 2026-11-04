use super::format::{read_header, Format};
use super::Key;
use anyhow::{bail, Context, Result};
use std::{
    cmp::{Ordering, Reverse},
    collections::{BinaryHeap, HashSet},
    env,
    fs::{self, File},
    io::{self, BufRead, BufReader, BufWriter, Write},
    time::Instant,
};
#[derive(Eq, PartialEq)]
struct Item {
    key: Key,
    source: usize,
    line: String,
}
impl Ord for Item {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key
            .cmp(&other.key)
            .then(self.source.cmp(&other.source))
    }
}
impl PartialOrd for Item {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
fn next(r: &mut BufReader<File>, source: usize, format: Format) -> Result<Option<Item>> {
    let mut line = String::new();
    if r.read_line(&mut line)? == 0 {
        return Ok(None);
    }
    if line.starts_with('#') || line.trim().is_empty() {
        bail!("invalid body line in input {source}")
    }
    let parts = format
        .fields(&line)
        .with_context(|| format!("invalid row in input {source}"))?;
    let key = Key {
        c1: parts[1].into(),
        c2: parts[3].into(),
        p1: parts[2].parse()?,
        p2: parts[4].parse()?,
        kind: parts[7].into(),
    };
    Ok(Some(Item { key, source, line }))
}
pub(crate) fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e:#}");
        std::process::exit(1)
    }
}
fn run() -> Result<()> {
    let args: Vec<_> = env::args().skip(1).collect();
    if args.is_empty() {
        bail!("usage: merge_pair_runs SORTED.pairs [SORTED.pairs ...] > merged.pairs")
    };
    let start = Instant::now();
    let mut readers = Vec::new();
    let mut reference = None;
    let mut paths = HashSet::new();
    let mut heap = BinaryHeap::new();
    let mut formats = Vec::new();
    let mut out = BufWriter::with_capacity(1024 * 1024, io::stdout().lock());
    for (i, path) in args.iter().enumerate() {
        if !paths.insert(fs::canonicalize(path)?) {
            bail!("duplicate input {path}")
        };
        let mut r = BufReader::with_capacity(
            1024 * 1024,
            File::open(path).with_context(|| format!("open {path}"))?,
        );
        let (header, format) =
            read_header(&mut r).with_context(|| format!("invalid header in {path}"))?;
        if !header
            .iter()
            .any(|x| x.starts_with("#sorted: chr1-chr2-pos1-pos2"))
        {
            bail!("missing sorted header in {path}")
        };
        if let Some(ref h) = reference {
            if *h != header {
                bail!("header mismatch: {path}")
            }
        } else {
            for line in &header {
                out.write_all(line.as_bytes())?
            }
            reference = Some(header)
        };
        if let Some(item) = next(&mut r, i, format)? {
            heap.push(Reverse(item))
        }
        formats.push(format);
        readers.push(r)
    }
    let mut count = 0u64;
    while let Some(Reverse(item)) = heap.pop() {
        out.write_all(item.line.as_bytes())?;
        count += 1;
        if count.is_multiple_of(10_000_000) {
            eprintln!(
                "event=merge_progress records={count} elapsed_seconds={:.3}",
                start.elapsed().as_secs_f64()
            );
        }
        let i = item.source;
        if let Some(next) = next(&mut readers[i], i, formats[i])? {
            if next.key < item.key {
                bail!("input not sorted: {} after {} rows", args[i], count)
            }
            heap.push(Reverse(next))
        }
    }
    out.flush()?;
    eprintln!(
        "event=merge_done inputs={} records={} elapsed_seconds={:.3}",
        args.len(),
        count,
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
