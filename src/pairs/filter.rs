use super::format::{read_header, Format};
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    env,
    fs::{self, File, OpenOptions},
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};
#[derive(Default, Deserialize, Serialize)]
struct Stats {
    input_records: u64,
    retained_records: u64,
    same_fragment: u64,
    dangling: u64,
    self_ligation: u64,
    unknown_mechanism: u64,
}

/// Combine independently filtered chunks into the workflow's usual totals.
pub(crate) fn combine_stats(inputs: &[PathBuf], output: &Path) -> Result<()> {
    let mut total = Stats::default();
    for input in inputs {
        let stats: Stats = serde_json::from_reader(File::open(input)?)
            .with_context(|| format!("read filter statistics {}", input.display()))?;
        macro_rules! add {
            ($($field:ident),+ $(,)?) => {
                $(total.$field = total.$field.checked_add(stats.$field)
                    .context("filter statistics overflow")?;)+
            };
        }
        add!(
            input_records,
            retained_records,
            same_fragment,
            dangling,
            self_ligation,
            unknown_mechanism
        );
    }
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)?;
    serde_json::to_writer_pretty(&mut file, &total)?;
    writeln!(file)?;
    Ok(())
}
fn load_fragments(path: &Path) -> Result<HashMap<String, Vec<u64>>> {
    let mut raw: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for (i, line) in BufReader::new(File::open(path)?).lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split_whitespace().collect();
        if f.len() != 3 {
            bail!("invalid BED line {}", i + 1)
        }
        let a: u64 = f[1].parse()?;
        let b: u64 = f[2].parse()?;
        if a > b {
            bail!("invalid fragment line {}", i + 1)
        }
        raw.entry(f[0].to_owned()).or_default().push((a, b));
    }
    if raw.is_empty() {
        bail!("empty fragment BED")
    }
    let mut result = HashMap::new();
    for (c, mut intervals) in raw {
        intervals.sort_unstable();
        let mut borders = Vec::with_capacity(intervals.len() + 1);
        borders.push(0);
        let mut prev = 0;
        for (start, end) in intervals {
            if start != prev {
                bail!("fragment BED for {c} is not contiguous from zero")
            };
            borders.push(end.checked_add(1).context("fragment overflow")?);
            prev = end
        }
        result.insert(c, borders);
    }
    Ok(result)
}

/// Validate the same BED format used by the filtering stage.
pub(crate) fn validate_fragments(path: &Path) -> Result<()> {
    load_fragments(path).map(|_| ())
}
fn find(frags: &HashMap<String, Vec<u64>>, c: &str, p: u64) -> Result<(usize, u64, u64)> {
    let borders = frags
        .get(c)
        .with_context(|| format!("chromosome {c} missing from fragment BED"))?;
    let idx = borders
        .partition_point(|x| *x <= p)
        .saturating_sub(1)
        .min(borders.len() - 2);
    Ok((idx, borders[idx], borders[idx + 1]))
}
fn run() -> Result<()> {
    let args: Vec<_> = env::args().collect();
    if args.len() != 3 {
        bail!("usage: filter_pairs FRAGMENTS.bed STATS.json < sorted.pairs > filtered.pairs")
    }
    let started = Instant::now();
    let frags = load_fragments(Path::new(&args[1]))?;
    eprintln!(
        "event=fragments_loaded chromosomes={} elapsed_seconds={:.3}",
        frags.len(),
        started.elapsed().as_secs_f64()
    );
    let stats_path = Path::new(&args[2]);
    if stats_path.exists() {
        bail!("stats output exists")
    };
    let partial = stats_path.with_extension(format!("partial.{}", std::process::id()));
    let stdin = io::stdin();
    let mut input = BufReader::with_capacity(1024 * 1024, stdin.lock());
    let mut out = BufWriter::with_capacity(1024 * 1024, io::stdout().lock());
    let mut stats = Stats::default();
    let (header, format) = read_header(&mut input)?;
    if format != Format::Standard {
        bail!(
            "restriction filtering requires 8-column pairs; input already has fragment annotations"
        )
    }
    for h in &header {
        if h.starts_with("#columns:") {
            out.write_all(Format::Fragments.columns_line().as_bytes())?;
        } else {
            out.write_all(h.as_bytes())?;
        }
    }
    let mut line = String::new();
    loop {
        line.clear();
        if input.read_line(&mut line)? == 0 {
            break;
        }
        let f: Vec<_> = line.trim_end_matches(['\r', '\n']).split('\t').collect();
        if f.len() != 8 {
            bail!("expected eight tab-separated columns")
        };
        let p1: u64 = f[2].parse()?;
        let p2: u64 = f[4].parse()?;
        if p1 == 0 || p2 == 0 {
            bail!("unmapped or zero position unsupported")
        };
        let a = find(&frags, f[1], p1)?;
        let b = find(&frags, f[3], p2)?;
        stats.input_records += 1;
        if f[1] == f[3] && a.0 == b.0 {
            stats.same_fragment += 1;
            match (f[5], f[6]) {
                ("+", "-") => stats.dangling += 1,
                ("-", "+") => stats.self_ligation += 1,
                _ => stats.unknown_mechanism += 1,
            }
        } else {
            stats.retained_records += 1;
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                line.trim_end_matches(['\r', '\n']),
                a.0,
                a.1,
                a.2,
                b.0,
                b.1,
                b.2
            )?
        }
        if stats.input_records % 10_000_000 == 0 {
            eprintln!(
                "event=filter_progress input_records={} same_fragment={} elapsed_seconds={:.3}",
                stats.input_records,
                stats.same_fragment,
                started.elapsed().as_secs_f64()
            )
        }
    }
    out.flush()?;
    let mut f = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&partial)?;
    serde_json::to_writer_pretty(&mut f, &stats)?;
    writeln!(f)?;
    f.sync_all()?;
    fs::hard_link(&partial, stats_path)?;
    fs::remove_file(&partial)?;
    eprintln!("event=filter_done input_records={} retained_records={} same_fragment={} elapsed_seconds={:.3}",stats.input_records,stats.retained_records,stats.same_fragment,started.elapsed().as_secs_f64());
    Ok(())
}
pub(crate) fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e:#}");
        std::process::exit(1)
    }
}
