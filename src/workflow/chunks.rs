//! Prepare completed mapping chunks while other chunks are still aligning.
use super::{config::Config, log, RunClock};
use crate::{
    mapping::align,
    pairs::filter,
    stage::{stage_command, Stage},
};
use anyhow::{bail, Context, Result};
use crossbeam_channel::unbounded;
use std::{
    fs::{self, File, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
    process::Stdio,
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering},
        Mutex,
    },
    thread,
};

pub(super) fn sort_one(
    c: &Config,
    work: &Path,
    out: &Path,
    i: usize,
    input: &Path,
    output: &Path,
) -> Result<()> {
    let log_path = out.join(format!("sort-{i:04}.log"));
    let f = File::create(&log_path)?;
    let status = stage_command(Stage::SortPairs)?
        .arg(&c.chromsizes)
        .arg(&c.assembly)
        .arg(output)
        .arg(work)
        .arg(c.resources.sort_threads.to_string())
        .arg(c.resources.sort_run_mib.to_string())
        .arg(input)
        .stderr(Stdio::from(f))
        .status()?;
    if !status.success() {
        let contents = fs::read_to_string(&log_path).unwrap_or_default();
        let detail = contents
            .lines()
            .rev()
            .find(|line| line.starts_with("error:"))
            .unwrap_or("see sort log");
        bail!(
            "sort chunk {i}: {status}: {detail}; inspect {}",
            log_path.display()
        )
    }
    Ok(())
}

fn filter_one(
    c: &Config,
    out: &Path,
    i: usize,
    input: &Path,
    output: &Path,
    cpus: Option<&libc::cpu_set_t>,
) -> Result<()> {
    let fragments = c
        .fragments
        .as_ref()
        .context("missing restriction fragments")?;
    let log_path = out.join(format!("filter-{i:04}.log"));
    let stats_path = out.join(format!("filter-{i:04}-stats.json"));
    let mut command = stage_command(Stage::FilterPairs)?;
    super::process::limit_command(&mut command, cpus);
    let status = command
        .arg(fragments)
        .arg(stats_path)
        .stdin(Stdio::from(File::open(input)?))
        .stdout(Stdio::from(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(output)?,
        ))
        .stderr(Stdio::from(File::create(&log_path)?))
        .status()?;
    if !status.success() {
        let contents = fs::read_to_string(&log_path).unwrap_or_default();
        let detail = contents
            .lines()
            .rev()
            .find(|line| line.starts_with("error:"))
            .unwrap_or("see filter log");
        bail!(
            "filter chunk {i}: {status}: {detail}; inspect {}",
            log_path.display()
        )
    }
    Ok(())
}

fn chunk_event(
    out: &Path,
    lock: &Mutex<()>,
    name: &str,
    chunk: usize,
    start: RunClock,
) -> Result<()> {
    let _guard = lock
        .lock()
        .map_err(|_| anyhow::anyhow!("event log lock poisoned"))?;
    let event = serde_json::json!({"event":name,"chunk":chunk,"elapsed_seconds":start.elapsed().as_secs_f64()});
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(out.join("events.jsonl"))?;
    writeln!(file, "{event}")?;
    eprintln!("{event}");
    Ok(())
}

fn prepared_paths(c: &Config, work: &Path, chunks: usize) -> Vec<PathBuf> {
    (0..chunks)
        .map(|i| {
            let prefix = if c.fragments.is_some() {
                "filtered"
            } else {
                "sorted"
            };
            work.join(format!("{prefix}-{i:04}.pairs"))
        })
        .collect()
}

#[allow(clippy::too_many_arguments)] // Shared context for preparing one mapped chunk.
fn prepare_one(
    c: &Config,
    work: &Path,
    out: &Path,
    i: usize,
    input: &Path,
    output: &Path,
    started: RunClock,
    cpus: Option<&libc::cpu_set_t>,
    event_lock: &Mutex<()>,
) -> Result<()> {
    chunk_event(out, event_lock, "sort_chunk_start", i, started)?;
    let sorted = work.join(format!("sorted-{i:04}.pairs"));
    sort_one(c, work, out, i, input, &sorted)?;
    chunk_event(out, event_lock, "sort_chunk_done", i, started)?;
    if c.fragments.is_some() {
        chunk_event(out, event_lock, "filter_chunk_start", i, started)?;
        filter_one(c, out, i, &sorted, output, cpus)?;
        chunk_event(out, event_lock, "filter_chunk_done", i, started)?;
    }
    Ok(())
}

fn finish_preparation(c: &Config, out: &Path, chunks: usize, started: RunClock) -> Result<()> {
    if c.fragments.is_some() {
        let stats: Vec<_> = (0..chunks)
            .map(|i| out.join(format!("filter-{i:04}-stats.json")))
            .collect();
        filter::combine_stats(&stats, &out.join("filter-stats.json"))?;
        fs::write(out.join("filter.log"), format!("Restriction filtering completed independently for {chunks} chunks; see filter-NNNN.log.\n"))?;
    }
    log(out, "sort_done", started)
}

/// Sort and optionally filter already mapped chunks on the current node.
pub(super) fn prepare_mapped(
    c: &Config,
    inputs: &[PathBuf],
    work: &Path,
    out: &Path,
    started: RunClock,
    cpus: Option<&libc::cpu_set_t>,
) -> Result<Vec<PathBuf>> {
    let chunks = inputs.len();
    if chunks == 0 {
        bail!("no mapped chunks to prepare")
    }
    let paths = prepared_paths(c, work, chunks);
    let event_lock = Mutex::new(());
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let errors = Mutex::new(Vec::new());
    log(out, "sort_start", started)?;
    thread::scope(|scope| {
        for _ in 0..c.resources.sort_jobs.min(chunks) {
            let event_lock = &event_lock;
            let next = &next;
            let failed = &failed;
            let errors = &errors;
            let paths = &paths;
            scope.spawn(move || loop {
                if failed.load(Ordering::Acquire) {
                    break;
                }
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= chunks {
                    break;
                }
                if let Err(e) = prepare_one(
                    c, work, out, i, &inputs[i], &paths[i], started, cpus, event_lock,
                ) {
                    failed.store(true, Ordering::Release);
                    errors.lock().unwrap().push(format!("chunk {i}: {e:#}"));
                }
            });
        }
    });
    let errors = errors.into_inner().unwrap();
    if !errors.is_empty() {
        bail!("chunk preparation failed: {}", errors.join("; "))
    }
    finish_preparation(c, out, chunks, started)?;
    Ok(paths)
}

pub(super) fn map_and_prepare(
    c: &Config,
    fastq: &align::FastqConfig,
    index: &Path,
    work: &Path,
    out: &Path,
    started: RunClock,
    cpus: Option<&libc::cpu_set_t>,
) -> Result<Vec<PathBuf>> {
    let chunks = fastq.chunks.context("missing index chunk count")?;
    let paths = prepared_paths(c, work, chunks);
    log(out, "alignment_start", started)?;
    log(out, "sort_start", started)?;
    let event_lock = Mutex::new(());
    let failed = AtomicBool::new(false);
    let result = thread::scope(|scope| -> Result<()> {
        let (send, receive) = unbounded::<(usize, PathBuf)>();
        let mut workers = Vec::new();
        for _ in 0..c.resources.sort_jobs.min(chunks) {
            let receive = receive.clone();
            let failed = &failed;
            let event_lock = &event_lock;
            let paths = &paths;
            workers.push(scope.spawn(move || -> Result<()> {
                while let Ok((i, input)) = receive.recv() {
                    if failed.load(Ordering::Acquire) {
                        break;
                    }
                    let result = prepare_one(
                        c, work, out, i, &input, &paths[i], started, cpus, event_lock,
                    );
                    if let Err(e) = result {
                        failed.store(true, Ordering::Release);
                        return Err(e);
                    }
                }
                Ok(())
            }));
        }
        drop(receive);
        let mapping = align::map_with_completed(
            fastq,
            index,
            work,
            out,
            None,
            &failed,
            c.keep_bam,
            |i, path| {
                chunk_event(out, &event_lock, "alignment_chunk_done", i, started)?;
                if failed.load(Ordering::Acquire) {
                    bail!("chunk preparation failed")
                }
                send.send((i, path.to_path_buf()))
                    .context("chunk preparation workers stopped")?;
                Ok(())
            },
        );
        drop(send);
        if mapping.is_err() {
            failed.store(true, Ordering::Release);
        }
        let log_result = if mapping.is_ok() {
            let _guard = event_lock
                .lock()
                .map_err(|_| anyhow::anyhow!("event log lock poisoned"))?;
            log(out, "alignment_done", started)
        } else {
            Ok(())
        };
        let mut preparation_error = None;
        for worker in workers {
            match worker.join() {
                Ok(Ok(())) => {}
                Ok(Err(e)) if preparation_error.is_none() => preparation_error = Some(e),
                Err(_) if preparation_error.is_none() => {
                    preparation_error = Some(anyhow::anyhow!("chunk preparation worker panicked"))
                }
                _ => {}
            }
        }
        if let Some(e) = preparation_error {
            return Err(e);
        }
        mapping?;
        log_result?;
        Ok(())
    });
    result?;
    finish_preparation(c, out, chunks, started)?;
    Ok(paths)
}
