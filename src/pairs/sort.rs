use super::format::Format;
use super::Key;
use anyhow::{bail, Context, Result};
use std::{
    collections::HashMap,
    env,
    fs::{self, File, OpenOptions},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    process::Stdio,
    thread,
    time::Instant,
};
struct Row {
    key: Key,
    line: String,
}
fn write_run(mut rows: Vec<Row>, path: PathBuf, header: String) -> Result<()> {
    rows.sort_by(|a, b| a.key.cmp(&b.key));
    let f = File::create(path)?;
    let mut out = BufWriter::with_capacity(1024 * 1024, f);
    out.write_all(header.as_bytes())?;
    for row in rows {
        out.write_all(row.line.as_bytes())?
    }
    out.flush()?;
    Ok(())
}
pub(crate) fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e:#}");
        std::process::exit(1)
    }
}
fn run() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 8 {
        bail!("usage: sort_pairs SIZES ASSEMBLY OUTPUT TMPDIR THREADS RUN_MIB INPUT.pairs")
    };
    let start = Instant::now();
    let output = Path::new(&args[3]);
    if output.exists() {
        bail!("output exists")
    };
    let threads: usize = args[5].parse()?;
    let limit: usize = args[6]
        .parse::<usize>()?
        .checked_mul(1024 * 1024)
        .context("run size overflow")?;
    if threads == 0 || limit < 1024 * 1024 {
        bail!("invalid thread count or run MiB")
    };
    let mut order = HashMap::new();
    let mut chroms = Vec::new();
    for line in BufReader::new(File::open(&args[1])?).lines() {
        let l = line?;
        let a: Vec<_> = l.split_whitespace().collect();
        if a.is_empty() {
            continue;
        }
        if a.len() != 2 {
            bail!("invalid chromsizes")
        };
        let n: u64 = a[1].parse()?;
        if n == 0 || order.insert(a[0].to_owned(), chroms.len()).is_some() {
            bail!("invalid chromsize")
        };
        chroms.push((a[0].to_owned(), n))
    }
    if chroms.is_empty() {
        bail!("empty chromsizes")
    };
    if args[2].contains(['\r', '\n']) {
        bail!("invalid assembly")
    };
    let mut header=format!("## pairs format v1.0.0\n#sorted: chr1-chr2-pos1-pos2\n#shape: upper triangle\n#genome_assembly: {}\n",args[2]);
    for (c, n) in &chroms {
        header += &format!("#chromsize: {c} {n}\n")
    }
    let scratch = Path::new(&args[4]).join(format!("rust-sort-{}", std::process::id()));
    fs::create_dir(&scratch).context("create exclusive scratch directory")?;
    let mut temp = output.to_path_buf();
    let name = output
        .file_name()
        .context("missing output filename")?
        .to_string_lossy();
    temp.set_file_name(format!(".{name}.{}.partial", std::process::id()));
    let result = (|| -> Result<()> {
        let mut jobs: Vec<thread::JoinHandle<Result<()>>> = Vec::new();
        let mut runs = Vec::new();
        let mut rows = Vec::new();
        let mut size = 0usize;
        let mut total = 0u64;
        let (mut body, mut columns) = (false, None);
        let mut reader = BufReader::with_capacity(1024 * 1024, File::open(&args[7])?);
        let mut line = String::new();
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            if line.starts_with('#') {
                if body {
                    bail!("header after body")
                };
                if let Some(meta) = line.strip_prefix("#chromsize:") {
                    let parts: Vec<_> = meta.split_whitespace().collect();
                    if parts.len() != 2 {
                        bail!("invalid #chromsize header");
                    }
                    let observed: u64 = parts[1].parse().context("invalid #chromsize length")?;
                    let idx = *order.get(parts[0]).with_context(|| {
                        format!("header chromosome absent from sizes: {}", parts[0])
                    })?;
                    if chroms[idx].1 != observed {
                        bail!("header length disagrees with reference for {}", parts[0]);
                    }
                }
                if line.starts_with("#columns:") {
                    if columns.is_some() {
                        bail!("duplicate #columns header")
                    }
                    let format = Format::from_columns(&line)?;
                    header.push_str(&format.columns_line());
                    columns = Some(format);
                }
                continue;
            }
            if line.trim().is_empty() {
                continue;
            }
            body = true;
            let format = columns.context("missing #columns header")?;
            let v = format.fields(&line)?;
            let p1: u64 = v[2].parse()?;
            let p2: u64 = v[4].parse()?;
            let a = *order
                .get(v[1])
                .with_context(|| format!("unknown chromosome {}", v[1]))?;
            let b = *order
                .get(v[3])
                .with_context(|| format!("unknown chromosome {}", v[3]))?;
            if p1 == 0
                || p2 == 0
                || p1 > chroms[a].1
                || p2 > chroms[b].1
                || !["+", "-"].contains(&v[5])
                || !["+", "-"].contains(&v[6])
            {
                bail!("invalid pair coordinates or strands")
            };
            let flip = (a, p1) > (b, p2);
            let (c1, x, c2, y, s1, s2) = if flip {
                (v[3], p2, v[1], p1, v[6], v[5])
            } else {
                (v[1], p1, v[3], p2, v[5], v[6])
            };
            let kind = v[7];
            if kind.len() != 2 || !kind.is_ascii() {
                bail!("pair_type must contain two ASCII characters")
            }
            let kind = if flip {
                kind.chars().rev().collect()
            } else {
                kind.to_owned()
            };
            let mut normalized = format!("{}\t{c1}\t{x}\t{c2}\t{y}\t{s1}\t{s2}\t{kind}", v[0]);
            if format == Format::Fragments {
                let (left, right) = if flip {
                    (&v[11..14], &v[8..11])
                } else {
                    (&v[8..11], &v[11..14])
                };
                for field in left.iter().chain(right) {
                    normalized.push('\t');
                    normalized.push_str(field);
                }
            }
            normalized.push('\n');
            size += normalized.len();
            rows.push(Row {
                key: Key {
                    c1: c1.into(),
                    c2: c2.into(),
                    p1: x,
                    p2: y,
                    kind,
                },
                line: normalized,
            });
            total += 1;
            if size >= limit {
                let path = scratch.join(format!("run-{:06}.pairs", runs.len()));
                runs.push(path.clone());
                eprintln!(
                    "event=run_started run={} records={total} elapsed_seconds={:.3}",
                    runs.len(),
                    start.elapsed().as_secs_f64()
                );
                let batch = std::mem::take(&mut rows);
                let h = header.clone();
                jobs.push(thread::spawn(move || write_run(batch, path, h)));
                size = 0;
                if jobs.len() >= threads {
                    jobs.remove(0)
                        .join()
                        .map_err(|_| anyhow::anyhow!("sort worker panicked"))??
                }
            }
        }
        if columns.is_none() {
            bail!("missing columns header")
        };
        if !rows.is_empty() || runs.is_empty() {
            let path = scratch.join(format!("run-{:06}.pairs", runs.len()));
            runs.push(path.clone());
            let h = header.clone();
            jobs.push(thread::spawn(move || write_run(rows, path, h)))
        }
        for job in jobs {
            job.join()
                .map_err(|_| anyhow::anyhow!("sort worker panicked"))??
        }
        eprintln!(
            "event=runs_done records={total} runs={} elapsed_seconds={:.3}",
            runs.len(),
            start.elapsed().as_secs_f64()
        );
        let file = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&temp)?;
        let status = crate::stage::stage_command(crate::stage::Stage::MergeRuns)?
            .args(&runs)
            .stdout(Stdio::from(file))
            .status()?;
        if !status.success() {
            bail!("run merge failed {status}")
        };
        fs::hard_link(&temp, output).context("publish output")?;
        eprintln!(
            "event=sort_done records={total} elapsed_seconds={:.3}",
            start.elapsed().as_secs_f64()
        );
        Ok(())
    })();
    if result.is_ok() {
        fs::remove_file(&temp)?;
        fs::remove_dir_all(&scratch).context("remove owned scratch directory")?;
    } else {
        eprintln!(
            "event=sort_failed scratch={} partial={}",
            scratch.display(),
            temp.display()
        );
    }
    result
}
