use super::format::{read_header, Format};
use anyhow::{bail, Context, Result};
use std::{
    cmp::{Ordering, Reverse},
    collections::{BinaryHeap, HashSet},
    env,
    fs::{self, File},
    io::{self, BufRead, BufReader, BufWriter, Seek, SeekFrom, Write},
    ops::Range,
    time::Instant,
};
#[derive(Eq, PartialEq)]
struct Item {
    line: String,
    c1: Range<usize>,
    c2: Range<usize>,
    kind: Range<usize>,
    p1: u64,
    p2: u64,
    source: usize,
}
impl Ord for Item {
    fn cmp(&self, o: &Self) -> Ordering {
        self.line[self.c1.clone()]
            .cmp(&o.line[o.c1.clone()])
            .then_with(|| self.line[self.c2.clone()].cmp(&o.line[o.c2.clone()]))
            .then(self.p1.cmp(&o.p1))
            .then(self.p2.cmp(&o.p2))
            .then_with(|| self.line[self.kind.clone()].cmp(&o.line[o.kind.clone()]))
            .then(self.source.cmp(&o.source))
    }
}
impl PartialOrd for Item {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
fn next(
    r: &mut BufReader<File>,
    source: usize,
    mut line: String,
    pos: &mut u64,
    end: u64,
    format: Format,
) -> Result<Option<Item>> {
    if *pos >= end {
        return Ok(None);
    }
    line.clear();
    if r.read_line(&mut line)? == 0 {
        return Ok(None);
    }
    *pos += line.len() as u64;
    if *pos > end {
        bail!("range end is not a line boundary");
    }
    if line.starts_with('#') || line.trim().is_empty() {
        bail!("invalid body line in input {source}")
    }
    let fields = format.fields(&line)?;
    let base = line.as_ptr() as usize;
    let c1 = fields[1];
    let p1 = fields[2].parse()?;
    let c2 = fields[3];
    let p2 = fields[4].parse()?;
    let kind = fields[7];
    let range = |s: &str| {
        let a = s.as_ptr() as usize - base;
        a..a + s.len()
    };
    let (c1, c2, kind) = (range(c1), range(c2), range(kind));
    Ok(Some(Item {
        line,
        c1,
        c2,
        kind,
        p1,
        p2,
        source,
    }))
}
fn run() -> Result<()> {
    let args: Vec<_> = env::args().skip(1).collect();
    if args.is_empty() {
        bail!("usage: merge_pair_ranges SORTED.pairs ...")
    }
    let start = Instant::now();
    let mut readers = Vec::new();
    let mut reference = None;
    let mut paths = HashSet::new();
    let ranges: Vec<(u64, u64)> = env::var(crate::pairs::MERGE_RANGES_ENV)
        .context("missing RUST_HIC_MERGE_RANGES")?
        .split(';')
        .map(|v| -> Result<_> {
            let (a, b) = v.split_once(',').context("bad range")?;
            Ok((a.parse()?, b.parse()?))
        })
        .collect::<Result<_>>()?;
    if ranges.len() != args.len() {
        bail!("range count mismatch")
    }
    let mut formats = Vec::new();
    let mut positions = Vec::new();
    let mut heap = BinaryHeap::new();
    let mut out = BufWriter::with_capacity(1024 * 1024, io::stdout().lock());
    for (i, path) in args.iter().enumerate() {
        if !paths.insert(fs::canonicalize(path)?) {
            bail!("duplicate input {path}")
        }
        let mut r = BufReader::with_capacity(1024 * 1024, File::open(path)?);
        let (header, format) =
            read_header(&mut r).with_context(|| format!("invalid header in {path}"))?;
        if !header
            .iter()
            .any(|s| s.starts_with("#sorted: chr1-chr2-pos1-pos2"))
        {
            bail!("missing sorted header in {path}")
        }
        if let Some(ref h) = reference {
            if *h != header {
                bail!("header mismatch: {path}")
            }
        } else {
            for s in &header {
                out.write_all(s.as_bytes())?
            }
            reference = Some(header)
        }
        let (start, end) = ranges[i];
        if start > end || end > r.get_ref().metadata()?.len() {
            bail!("invalid range")
        }
        r.seek(SeekFrom::Start(start))?;
        let mut position = start;
        if let Some(item) = next(&mut r, i, String::new(), &mut position, end, format)? {
            heap.push(Reverse(item))
        }
        formats.push(format);
        positions.push(position);
        readers.push(r)
    }
    let mut count = 0u64;
    let mut spare: Vec<String> = (0..readers.len()).map(|_| String::new()).collect();
    while let Some(Reverse(item)) = heap.pop() {
        out.write_all(item.line.as_bytes())?;
        count += 1;
        let i = item.source;
        if let Some(n) = next(
            &mut readers[i],
            i,
            std::mem::take(&mut spare[i]),
            &mut positions[i],
            ranges[i].1,
            formats[i],
        )? {
            if n < item {
                bail!("input not sorted: {}", args[i])
            }
            spare[i] = item.line;
            heap.push(Reverse(n))
        }
        if count.is_multiple_of(10_000_000) {
            eprintln!(
                "event=merge_progress records={count} elapsed_seconds={:.3}",
                start.elapsed().as_secs_f64()
            )
        }
    }
    out.flush()?;
    eprintln!(
        "event=merge_done inputs={} records={count} elapsed_seconds={:.3}",
        args.len(),
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
pub(crate) fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e:#}");
        std::process::exit(1)
    }
}
