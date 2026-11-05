//! End-to-end workflow and its shared configuration.
mod chunks;
pub(crate) mod config;
pub(crate) mod digest;
pub(crate) mod pairix;
pub(crate) mod preflight;
pub(crate) mod process;
pub(crate) mod report;
pub(crate) mod slurm;
pub(crate) mod stage;

use crate::mapping::align;
use crate::stage::{stage_command, Stage};
use anyhow::{bail, Context, Result};
use config::Config;
use process::{affinity_mask, current_affinity, exec, pipe, set_affinity};
use std::{
    env,
    fs::{self, File, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc, Mutex,
    },
    thread,
    time::{Duration, Instant},
};

/// Monotonic within this process, with an offset for jobs started earlier.
#[derive(Clone, Copy)]
pub(super) struct RunClock {
    start: Instant,
    offset: Duration,
}
impl RunClock {
    pub(crate) fn new() -> Self {
        Self {
            start: Instant::now(),
            offset: Duration::ZERO,
        }
    }
    fn with_offset(offset: Duration) -> Self {
        Self {
            start: Instant::now(),
            offset,
        }
    }
    fn elapsed(self) -> Duration {
        self.offset.saturating_add(self.start.elapsed())
    }
}

pub(crate) fn log(dir: &Path, event: &str, start: RunClock) -> Result<()> {
    let record = serde_json::json!({"event":event,"elapsed_seconds":start.elapsed().as_secs_f64()});
    let mut f = OpenOptions::new()
        .create(true)
        .append(true)
        .open(dir.join("events.jsonl"))?;
    writeln!(f, "{}", record)?;
    eprintln!("{}", record);
    Ok(())
}
fn path_from(base: &Path, p: &Path) -> PathBuf {
    if p.is_absolute() {
        p.into()
    } else {
        base.join(p)
    }
}
pub(crate) fn run(path: &Path, check: bool, keep_pairs: bool, keep_bam: bool) -> Result<()> {
    let config_path =
        fs::canonicalize(path).with_context(|| format!("open configuration {}", path.display()))?;
    let base = config_path.parent().context("config parent")?;
    let mut c: Config = serde_yaml::from_reader(File::open(&config_path)?)
        .with_context(|| format!("parse configuration {}", config_path.display()))?;
    c.keep_pairs |= keep_pairs;
    c.keep_bam |= keep_bam;
    if check {
        let c = preflight::prepare(c, base, false)?;
        println!(
            "Preflight passed\n  Chromosome sizes: {}\n  Output: {}\n  Scratch: {}",
            c.chromsizes.display(),
            c.output_dir.display(),
            c.scratch_dir.display()
        );
        if let Some(fastq) = c.fastq {
            println!(
                "  FASTQ: {} and {}\n  Aligner: {}\n  Indexing threads per FASTQ: {}",
                fastq.read1.display(),
                fastq.read2.display(),
                fastq.aligner,
                fastq.index_threads_per_mate
            );
        } else {
            println!("  Pairs inputs: {}", c.input_pairs.len());
        }
        return Ok(());
    }
    run_config(c, base, true, false)
}
pub(crate) fn run_config(
    c: Config,
    base: &Path,
    retain_scratch: bool,
    pairs_only: bool,
) -> Result<()> {
    let started = RunClock::new();
    let mut c = preflight::prepare(c, base, pairs_only)?;
    let r = c.resources.clone();
    let postprocess_cpus = affinity_mask(r.postprocess_cpus.as_deref())?;
    fs::create_dir(&c.output_dir)
        .with_context(|| format!("create new output directory {}", c.output_dir.display()))?;
    let out = fs::canonicalize(&c.output_dir)?;
    let outcome = (|| -> Result<()> {
        serde_yaml::to_writer(File::create(out.join("config.resolved.yaml"))?, &c)?;
        log(&out, "pipeline_start", started)?;
        let work = tempfile::Builder::new()
            .prefix("rust-hic-")
            .tempdir_in(&c.scratch_dir)?
            .keep();
        fs::write(out.join("work-dir.txt"), format!("{}\n", work.display()))?;
        if let Some(digest) = &c.digest {
            log(&out, "digest_start", started)?;
            let bed = work.join("restriction-fragments.bed");
            digest::digest(digest, &c.chromsizes, &bed)?;
            c.fragments = Some(bed);
            log(&out, "digest_done", started)?;
        }
        let prefiltered = c.fastq.is_some() && c.fragments.is_some();
        let sorted = if let Some(fastq) = &mut c.fastq {
            let indexer = env::current_exe()?;
            log(&out, "fastq_index_start", started)?;
            let idx = align::index(
                fastq,
                &indexer,
                &work.join("fastq-index"),
                &out.join("fastq-index.log"),
                false,
                false,
            )?;
            log(&out, "fastq_index_done", started)?;
            let sorted = chunks::map_and_prepare(
                &c,
                c.fastq.as_ref().context("missing FASTQ configuration")?,
                &idx,
                &work,
                &out,
                started,
                postprocess_cpus.as_ref(),
            )?;
            if c.keep_bam {
                let count = c
                    .fastq
                    .as_ref()
                    .and_then(|f| f.chunks)
                    .context("missing chunk count")?;
                let parts: Vec<_> = (0..count)
                    .map(|i| work.join(format!("align-{i:04}/alignments.bam")))
                    .collect();
                crate::mapping::bam::merge(&parts, &out.join("alignments.bam"))?;
                log(&out, "bam_merge_done", started)?;
            }
            sorted
        } else if c.input_sorted {
            c.input_pairs.clone()
        } else {
            log(&out, "sort_start", started)?;
            let paths = Arc::new(
                (0..c.input_pairs.len())
                    .map(|i| work.join(format!("sorted-{i:04}.pairs")))
                    .collect::<Vec<_>>(),
            );
            let next = AtomicUsize::new(0);
            let errors = Mutex::new(Vec::new());
            thread::scope(|scope| {
                for _ in 0..r.sort_jobs.min(c.input_pairs.len()) {
                    let paths = Arc::clone(&paths);
                    let next = &next;
                    let errors = &errors;
                    let c = &c;
                    let out = &out;
                    let work = &work;
                    scope.spawn(move || loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= paths.len() {
                            break;
                        }
                        let result = (|| -> Result<()> {
                            chunks::sort_one(c, work, out, i, &c.input_pairs[i], &paths[i])?;
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
                bail!("{}", errors.join("; "))
            }
            log(&out, "sort_done", started)?;
            Arc::try_unwrap(paths).unwrap()
        };
        finish_stages(
            &c,
            &work,
            &out,
            &sorted,
            prefiltered,
            pairs_only,
            started,
            &postprocess_cpus,
        )?;
        if retain_scratch {
            fs::remove_dir_all(&work).with_context(|| {
                format!("remove completed scratch workspace {}", work.display())
            })?;
            fs::remove_file(out.join("work-dir.txt"))?;
        }
        log(&out, "pipeline_done", started)?;
        fs::write(
            out.join("status.json"),
            format!(
                "{{\"status\":\"complete\",\"elapsed_seconds\":{:.3}}}\n",
                started.elapsed().as_secs_f64()
            ),
        )?;
        Ok(())
    })();
    if let Err(error) = &outcome {
        let record = serde_json::json!({"status":"failed","elapsed_seconds":started.elapsed().as_secs_f64(),"error":format!("{error:#}")});
        let _ = fs::write(out.join("failure.json"), format!("{record}\n"));
        let _ = log(&out, "pipeline_failed", started);
    }
    outcome
}

#[allow(clippy::too_many_arguments)] // Shared workflow state for local and distributed execution.
pub(super) fn finish_stages(
    c: &Config,
    work: &Path,
    out: &Path,
    sorted: &[PathBuf],
    prefiltered: bool,
    pairs_only: bool,
    started: RunClock,
    postprocess_cpus: &Option<libc::cpu_set_t>,
) -> Result<()> {
    let r = c.resources.clone();
    log(out, "merge_filter_dedup_start", started)?;
    let mut merger = stage_command(Stage::MergePairs)?;
    merger
        .arg(r.merge_workers.to_string())
        .arg(c.merge_scratch_dir.as_ref().unwrap_or(&c.scratch_dir))
        .args(sorted);
    let mut cmds = vec![("merge".into(), merger)];
    if let Some(bed) = c.fragments.as_ref().filter(|_| !prefiltered) {
        let mut f = stage_command(Stage::FilterPairs)?;
        f.arg(bed).arg(out.join("filter-stats.json"));
        cmds.push(("filter".into(), f))
    }
    let mut d = stage_command(Stage::DedupPairs)?;
    d.arg(r.dedup_workers.to_string())
        .arg(c.dedup_max_mismatch.to_string());
    cmds.push(("dedup".into(), d));
    if !pairs_only {
        let mut stats = stage_command(Stage::PairsStats)?;
        stats.arg(out.join("contact-stats.json"));
        cmds.push(("contact".into(), stats));
    }
    let pairs = if pairs_only || c.keep_pairs {
        out.join("filtered-nodups.pairs")
    } else {
        work.join("filtered-nodups.pairs")
    };
    pipe(cmds, &pairs, out, postprocess_cpus.as_ref())?;
    log(out, "merge_filter_dedup_done", started)?;
    if !pairs_only {
        log(out, "cool_start", started)?;
        let base_cool = work.join("base.cool");
        let mut cmd = stage_command(Stage::BinPairs)?;
        cmd.arg(&c.chromsizes)
            .arg(c.resolutions[0].to_string())
            .arg(&c.assembly)
            .arg(&pairs)
            .arg(&base_cool)
            .arg(work)
            .arg(r.matrix_threads.to_string());
        exec(cmd, out, "pairs-to-cool", postprocess_cpus.as_ref())?;
        log(out, "cool_done", started)?;
        log(out, "mcool_start", started)?;
        let mut cmd = stage_command(Stage::BuildMcool)?;
        cmd.arg(&base_cool)
            .arg(&c.chromsizes)
            .arg(out.join("contacts.mcool"))
            .arg(&c.assembly)
            .arg(
                c.resolutions
                    .iter()
                    .map(u64::to_string)
                    .collect::<Vec<_>>()
                    .join(","),
            )
            .arg(r.matrix_threads.to_string())
            .env(
                crate::matrix::BALANCE_MAX_ITERS_ENV,
                c.balance_max_iterations.to_string(),
            );
        exec(cmd, out, "cool-to-mcool", postprocess_cpus.as_ref())?;
        log(out, "mcool_done", started)?;
        if c.compressed_pairs {
            log(out, "pairix_start", started)?;
            let original_affinity = if let Some(mask) = &postprocess_cpus {
                let original = current_affinity()?;
                set_affinity(mask)?;
                Some(original)
            } else {
                None
            };
            let result = pairix::pack(
                &pairs,
                &out.join("filtered-nodups.pairs.gz"),
                r.compression_threads,
                Some(&out.join("pairix.log")),
            );
            let restored = original_affinity.as_ref().map(set_affinity).transpose();
            result?;
            restored?;
            log(out, "pairix_done", started)?;
        }
        report::generate(out, c, started.elapsed().as_secs_f64())?;
        eprintln!("report: {}", out.join("report.html").display());
    }
    Ok(())
}
