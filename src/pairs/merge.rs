use super::format::{read_header, Format};
use super::Key;
use anyhow::{bail, Context, Result};
use std::{
    collections::HashSet,
    env,
    fs::{self, File},
    io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
    process::Stdio,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc, Mutex,
    },
    thread,
    time::Instant,
};
fn parse(line: &str, format: Format) -> Result<Key> {
    let f = format.fields(line)?;
    Ok(Key {
        c1: f[1].into(),
        c2: f[3].into(),
        p1: f[2].parse()?,
        p2: f[4].parse()?,
        kind: f[7].into(),
    })
}
fn body_start(r: &mut BufReader<File>) -> Result<u64> {
    r.seek(SeekFrom::Start(0))?;
    loop {
        let pos = r.stream_position()?;
        let b = r.fill_buf()?;
        if b.first() != Some(&b'#') {
            return Ok(pos);
        }
        let mut line = String::new();
        r.read_line(&mut line)?;
    }
}
fn read_line(r: &mut BufReader<File>) -> Result<Option<String>> {
    let mut line = String::new();
    if r.read_line(&mut line)? == 0 {
        Ok(None)
    } else {
        Ok(Some(line))
    }
}
fn lower_bound(
    r: &mut BufReader<File>,
    begin: u64,
    end: u64,
    target: &Key,
    format: Format,
) -> Result<u64> {
    let (mut lo, mut hi) = (begin, end);
    while hi - lo > 4096 {
        let mid = lo + (hi - lo) / 2;
        r.seek(SeekFrom::Start(mid))?;
        read_line(r)?;
        let p = r.stream_position()?;
        if p >= hi {
            hi = mid;
            continue;
        }
        let line = read_line(r)?;
        if line
            .as_ref()
            .map(|s| parse(s, format))
            .transpose()?
            .as_ref()
            .is_none_or(|key| key >= target)
        {
            hi = p
        } else {
            lo = r.stream_position()?
        }
    }
    r.seek(SeekFrom::Start(lo))?;
    loop {
        let p = r.stream_position()?;
        if p >= end {
            return Ok(end);
        }
        let line = read_line(r)?.context("unexpected end finding boundary")?;
        if parse(&line, format)? >= *target {
            return Ok(p);
        }
    }
}
// Each worker checks its own input range. Check adjoining rows too, so an
// inversion exactly on a partition boundary is not missed.
fn validate_boundary(
    r: &mut BufReader<File>,
    begin: u64,
    end: u64,
    pos: u64,
    format: Format,
) -> Result<()> {
    if pos <= begin || pos >= end {
        return Ok(());
    }
    let mut width = (pos - begin).min(65536);
    let previous = loop {
        r.seek(SeekFrom::Start(pos - width))?;
        let mut bytes = vec![0u8; width as usize];
        r.read_exact(&mut bytes)?;
        if bytes.last() != Some(&b'\n') {
            bail!("partition is not at a line boundary")
        }
        let start = bytes[..bytes.len() - 1]
            .iter()
            .rposition(|b| *b == b'\n')
            .map(|i| i + 1);
        if let Some(i) = start {
            break String::from_utf8(bytes[i..].to_vec())?;
        }
        if width == pos - begin {
            break String::from_utf8(bytes)?;
        }
        width = (width.saturating_mul(2)).min(pos - begin);
    };
    r.seek(SeekFrom::Start(pos))?;
    let next = read_line(r)?.context("missing record at partition boundary")?;
    if parse(&previous, format)? > parse(&next, format)? {
        bail!("input is not sorted at a range boundary")
    }
    Ok(())
}
fn run() -> Result<()> {
    let args: Vec<_> = env::args_os().collect();
    if args.len() < 4 {
        bail!("usage: merge_pairs WORKERS TMPDIR SORTED.pairs [SORTED.pairs ...] > merged.pairs")
    }
    let workers: usize = args[1]
        .to_str()
        .context("non-UTF-8 worker count")?
        .parse()?;
    if !(1..=64).contains(&workers) {
        bail!("workers must be 1..64")
    };
    let tmp = Path::new(&args[2]);
    if !tmp.is_dir() {
        bail!("temporary directory does not exist")
    };
    let start = Instant::now();
    let paths: Vec<PathBuf> = args[3..]
        .iter()
        .map(fs::canonicalize)
        .collect::<io::Result<_>>()?;
    let mut unique = HashSet::new();
    if !paths.iter().all(|p| unique.insert(p)) {
        bail!("duplicate input")
    }
    let mut first = BufReader::with_capacity(1024 * 1024, File::open(&paths[0])?);
    let (reference_header, format) = read_header(&mut first)?;
    let begin = first.stream_position()?;
    let end = first.get_ref().metadata()?.len();
    let mut cuts = Vec::new();
    for i in 1..workers {
        first.seek(SeekFrom::Start(
            begin + (end - begin) * i as u64 / workers as u64,
        ))?;
        read_line(&mut first)?;
        if let Some(line) = read_line(&mut first)? {
            cuts.push(parse(&line, format)?);
        }
    }
    cuts.sort_unstable();
    cuts.dedup();
    let mut offsets = Vec::new();
    for path in &paths {
        let mut r = BufReader::with_capacity(1024 * 1024, File::open(path)?);
        let (header, input_format) =
            read_header(&mut r).with_context(|| format!("invalid header in {}", path.display()))?;
        if header != reference_header || input_format != format {
            bail!("header mismatch: {}", path.display())
        }
        let begin = r.stream_position()?;
        let end = r.get_ref().metadata()?.len();
        let mut o = vec![begin];
        for cut in &cuts {
            o.push(lower_bound(&mut r, begin, end, cut, format)?)
        }
        o.push(end);
        if o.windows(2).any(|w| w[0] > w[1]) {
            bail!("range boundaries are not ordered: {}", path.display())
        }
        for &boundary in &o[1..o.len() - 1] {
            validate_boundary(&mut r, begin, end, boundary, format)
                .with_context(|| format!("invalid input {}", path.display()))?;
        }
        offsets.push(o)
    }
    let segments = cuts.len() + 1;
    eprintln!(
        "event=range_index_done segments={segments} elapsed_seconds={:.3}",
        start.elapsed().as_secs_f64()
    );
    let dir = tempfile::Builder::new()
        .prefix("rust-hic-merge-")
        .tempdir_in(tmp)?;
    let next = AtomicUsize::new(0);
    let errors = Mutex::new(Vec::<String>::new());
    let n = workers.min(segments);
    let all = Arc::new(offsets);
    thread::scope(|scope| {
        for _ in 0..n {
            let errors = &errors;
            let next = &next;
            let offsets = Arc::clone(&all);
            let paths = &paths;
            let dir = dir.path();
            scope.spawn(move || loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= segments {
                    break;
                }
                let result = (|| -> Result<()> {
                    let output = File::create(dir.join(format!("{i:05}.pairs")))?;
                    let logpath = dir.join(format!("{i:05}.log"));
                    let log = File::create(&logpath)?;
                    let ranges = offsets
                        .iter()
                        .map(|o| format!("{},{}", o[i], o[i + 1]))
                        .collect::<Vec<_>>()
                        .join(";");
                    let s = crate::stage::stage_command(crate::stage::Stage::MergeRanges)?
                        .args(paths)
                        .env(crate::pairs::MERGE_RANGES_ENV, ranges)
                        .stdout(Stdio::from(output))
                        .stderr(Stdio::from(log))
                        .status()?;
                    if !s.success() {
                        let msg = fs::read_to_string(logpath).unwrap_or_default();
                        bail!("segment {i} exit {s}: {msg}")
                    }
                    Ok(())
                })();
                if let Err(e) = result {
                    errors.lock().unwrap().push(format!("{e:#}"))
                }
            });
        }
    });
    let errors = errors.into_inner().unwrap();
    if !errors.is_empty() {
        bail!("range merging failed: {}", errors.join("; "))
    }
    let bytes: u64 = (0..segments)
        .map(|i| fs::metadata(dir.path().join(format!("{i:05}.pairs"))).map(|m| m.len()))
        .collect::<io::Result<Vec<_>>>()?
        .iter()
        .sum();
    eprintln!(
        "event=range_merge_done segments={segments} bytes={bytes} elapsed_seconds={:.3}",
        start.elapsed().as_secs_f64()
    );
    let stdout = io::stdout();
    let mut out = BufWriter::with_capacity(1024 * 1024, stdout.lock());
    for i in 0..segments {
        let mut r = BufReader::with_capacity(
            1024 * 1024,
            File::open(dir.path().join(format!("{i:05}.pairs")))?,
        );
        if i > 0 {
            let pos = body_start(&mut r)?;
            r.seek(SeekFrom::Start(pos))?;
        }
        let mut limited = (&mut r).take(u64::MAX);
        io::copy(&mut limited, &mut out)?;
    }
    out.flush()?;
    eprintln!(
        "event=range_stream_done elapsed_seconds={:.3}",
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
