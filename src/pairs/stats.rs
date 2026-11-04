use super::format::Format;
use anyhow::{bail, Context, Result};
use serde::Serialize;
use std::{
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
}
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
        if f[1] == f[3] {
            stats.cis += 1;
            if a.abs_diff(b) >= 20_000 {
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
