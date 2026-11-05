use super::format::Format;
use anyhow::{bail, Context, Result};
use serde::Serialize;
use std::{
    collections::BTreeMap,
    env,
    fs::{self, OpenOptions},
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::Path,
    time::Instant,
};
#[derive(Default, Serialize)]
struct Stats {
    total_nodups: u64,
    cis: u64,
    trans: u64,
    cis_20kb_plus: u64,
    cis_shorter_than_20kb: u64,
    pair_types: BTreeMap<String, u64>,
    strand_orientations: BTreeMap<String, u64>,
    cis_distance_bins: Vec<DistanceBin>,
}

#[derive(Serialize)]
struct DistanceBin {
    start_bp: u64,
    end_bp: Option<u64>,
    count: u64,
}

const DISTANCE_STARTS: [u64; 20] = [
    0,
    1_000,
    2_000,
    5_000,
    10_000,
    20_000,
    50_000,
    100_000,
    200_000,
    500_000,
    1_000_000,
    2_000_000,
    5_000_000,
    10_000_000,
    20_000_000,
    50_000_000,
    100_000_000,
    200_000_000,
    500_000_000,
    1_000_000_000,
];
fn run() -> Result<()> {
    let args: Vec<_> = env::args().collect();
    if args.len() != 2 {
        bail!("usage: pairs_stats STATS.json < nodups.pairs > nodups.pairs")
    }
    let start = Instant::now();
    let dest = Path::new(&args[1]);
    if dest.exists() {
        bail!("stats output exists")
    };
    let temp = dest.with_extension(format!("partial.{}", std::process::id()));
    let stdin = io::stdin();
    let mut input = BufReader::with_capacity(1024 * 1024, stdin.lock());
    let mut out = BufWriter::with_capacity(1024 * 1024, io::stdout().lock());
    let mut stats = Stats::default();
    let mut orientations = [0u64; 5];
    let mut distance_counts = [0u64; DISTANCE_STARTS.len()];
    let mut line = String::new();
    let mut body = false;
    let mut format = None;
    loop {
        line.clear();
        if input.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with('#') {
            if body {
                bail!("header inside body")
            };
            if line.starts_with("#columns:") {
                if format.is_some() {
                    bail!("duplicate #columns header")
                }
                format = Some(Format::from_columns(&line)?);
            }
            out.write_all(line.as_bytes())?;
            continue;
        }
        body = true;
        let f = format.context("missing #columns header")?.fields(&line)?;
        let a: u64 = f[2].parse().context("pos1")?;
        let b: u64 = f[4].parse().context("pos2")?;
        if a == 0 || b == 0 || f[1] == "!" || f[3] == "!" {
            bail!("unmapped record unsupported")
        };
        stats.total_nodups += 1;
        if let Some(count) = stats.pair_types.get_mut(f[7]) {
            *count += 1;
        } else {
            stats.pair_types.insert(f[7].to_owned(), 1);
        }
        orientations[match (f[5], f[6]) {
            ("+", "+") => 0,
            ("+", "-") => 1,
            ("-", "+") => 2,
            ("-", "-") => 3,
            _ => 4,
        }] += 1;
        if f[1] == f[3] {
            stats.cis += 1;
            let distance = a.abs_diff(b);
            let bin = DISTANCE_STARTS.partition_point(|start| *start <= distance) - 1;
            distance_counts[bin] += 1;
            if distance >= 20_000 {
                stats.cis_20kb_plus += 1
            } else {
                stats.cis_shorter_than_20kb += 1
            }
        } else {
            stats.trans += 1
        }
        out.write_all(line.as_bytes())?;
        if stats.total_nodups % 10_000_000 == 0 {
            eprintln!(
                "event=stats_progress records={} elapsed_seconds={:.3}",
                stats.total_nodups,
                start.elapsed().as_secs_f64()
            )
        }
    }
    if format.is_none() {
        bail!("missing supported columns header")
    };
    for (name, count) in ["++", "+-", "-+", "--", "other"]
        .into_iter()
        .zip(orientations)
    {
        stats.strand_orientations.insert(name.into(), count);
    }
    stats.cis_distance_bins = DISTANCE_STARTS
        .iter()
        .enumerate()
        .map(|(i, &start_bp)| DistanceBin {
            start_bp,
            end_bp: DISTANCE_STARTS.get(i + 1).copied(),
            count: distance_counts[i],
        })
        .collect();
    out.flush()?;
    let mut f = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temp)?;
    serde_json::to_writer_pretty(&mut f, &stats)?;
    writeln!(f)?;
    f.sync_all()?;
    fs::hard_link(&temp, dest)?;
    fs::remove_file(&temp)?;
    eprintln!(
        "event=stats_done records={} elapsed_seconds={:.3}",
        stats.total_nodups,
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
