//! Submit the index, alignment array, and finishing jobs.
use super::state::{hash_file, now_ms, record_failure, write_json, Jobs, Plan, SCHEMA};
use crate::{
    cli::run::SlurmOptions,
    mapping::align,
    workflow::{config::Config, preflight},
};
use anyhow::{bail, Context, Result};
use std::{
    env,
    ffi::OsStr,
    fs::{self, File, OpenOptions},
    io::Write,
    os::unix::fs::PermissionsExt,
    path::{Path, PathBuf},
    process::Command,
};

const WORK_PREFIX: &str = "rust-hic-slurm-";

fn shell_quote(value: &OsStr) -> Result<String> {
    let value = value.to_str().context("Slurm script paths must be UTF-8")?;
    Ok(format!("'{}'", value.replace('\'', "'\\''")))
}

fn log_pattern(path: &Path) -> Result<String> {
    let parent = path
        .parent()
        .context("Slurm log has no parent directory")?
        .to_str()
        .context("Slurm log paths must be UTF-8")?;
    if parent.contains(['\\', '\r', '\n']) {
        bail!("Slurm output paths cannot contain backslashes or line breaks")
    }
    let parent = parent.replace('%', "%%");
    let name = path
        .file_name()
        .and_then(OsStr::to_str)
        .context("invalid Slurm log name")?;
    Ok(format!("{parent}/{name}"))
}

fn script(work: &Path, phase: &str, array: bool) -> Result<PathBuf> {
    let path = work.join(format!("slurm-{phase}.sh"));
    let exe = shell_quote(work.join("rust-hic").as_os_str())?;
    let work = shell_quote(work.as_os_str())?;
    let chunk = if array {
        " \"${SLURM_ARRAY_TASK_ID:?missing array task id}\""
    } else {
        ""
    };
    let mut f = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)?;
    writeln!(
        f,
        "#!/bin/sh\nset -eu\nexec {exe} slurm-task {phase} {work}{chunk}"
    )?;
    f.flush()?;
    Ok(path)
}

fn validate_extra_args(args: &[String]) -> Result<()> {
    const MANAGED: &[&str] = &[
        "--array",
        "--dependency",
        "--kill-on-invalid-dep",
        "--output",
        "--error",
        "--open-mode",
        "--job-name",
        "--nodes",
        "--ntasks",
        "--ntasks-per-node",
        "--cpus-per-task",
        "--mem",
        "--mem-per-cpu",
        "--mem-per-gpu",
        "--time",
        "--account",
        "--partition",
        "--parsable",
        "--wrap",
        "--wait",
        "--chdir",
        "--export",
        "--clusters",
        "--cluster",
        "--hold",
        "--test-only",
        "--requeue",
        "--no-requeue",
        "--help",
        "--version",
    ];
    for arg in args {
        let Some((key, value)) = arg.split_once('=') else {
            if arg == "--exclusive" || arg == "--no-requeue" {
                continue;
            }
            bail!("--sbatch-arg must use --name=value (or --exclusive/--no-requeue): {arg}")
        };
        if !key.starts_with("--") || key.len() <= 2 || value.is_empty() || MANAGED.contains(&key) {
            bail!("--sbatch-arg conflicts with a managed Slurm option or is invalid: {arg}")
        }
    }
    Ok(())
}

fn job_cpus(c: &Config, args: &SlurmOptions) -> Result<(usize, usize, usize)> {
    let fastq = c
        .fastq
        .as_ref()
        .context("Slurm execution requires FASTQ input")?;
    let index = fastq
        .index_threads_per_mate
        .checked_mul(2)
        .and_then(|n| n.checked_add(1))
        .context("index CPU count overflow")?;
    let map = fastq
        .threads_per_job
        .checked_add(1)
        .context("map CPU count overflow")?;
    let sort = c
        .resources
        .sort_threads
        .checked_add(1)
        .and_then(|n| n.checked_mul(c.resources.sort_jobs))
        .context("sort CPU count overflow")?;
    let merge = c
        .resources
        .merge_workers
        .checked_add(c.resources.dedup_workers)
        .and_then(|n| n.checked_add(3))
        .context("merge CPU count overflow")?;
    let matrix = c
        .resources
        .matrix_threads
        .checked_add(1)
        .context("matrix CPU count overflow")?;
    let requested = (
        args.index_cpus.unwrap_or(index),
        args.map_cpus.unwrap_or(map),
        args.finish_cpus.unwrap_or(sort.max(merge).max(matrix)),
    );
    if requested.0 == 0 || requested.1 == 0 || requested.2 == 0 {
        bail!("Slurm CPU requests must be positive")
    }
    Ok(requested)
}

fn estimated_chunks(c: &align::FastqConfig) -> Result<usize> {
    if let Some(index) = &c.index_dir {
        return align::validate_reusable_index(c, index);
    }
    if let Some(chunks) = c.chunks {
        return Ok(chunks);
    }
    let bytes = fs::metadata(&c.read1)?
        .len()
        .max(fs::metadata(&c.read2)?.len());
    let chunks = usize::try_from(bytes.div_ceil(c.chunk_size_bytes))?;
    if !(1..=256).contains(&chunks) {
        bail!("estimated FASTQ chunks must be between 1 and 256")
    }
    Ok(chunks)
}

#[derive(Clone, Copy)]
enum JobStage {
    Index,
    Map,
    Finish,
}

impl JobStage {
    fn name(self, group_id: &str) -> String {
        let role = match self {
            Self::Index | Self::Finish => "main",
            Self::Map => "mapping",
        };
        format!("rust-hic-{group_id}-{role}")
    }

    fn description(self) -> &'static str {
        match self {
            Self::Index => "index",
            Self::Map => "alignment array",
            Self::Finish => "final processing",
        }
    }
}

#[allow(clippy::too_many_arguments)] // One scheduler request assembled from stage settings.
fn sbatch(
    args: &SlurmOptions,
    stage: JobStage,
    group_id: &str,
    log: &Path,
    cpus: usize,
    mem: &str,
    time: &str,
    dependency: Option<&str>,
    array: Option<&str>,
    script: &Path,
) -> Result<String> {
    let name = stage.name(group_id);
    let mut cmd = Command::new("sbatch");
    cmd.arg("--parsable")
        .arg(format!("--job-name={name}"))
        .arg("--nodes=1")
        .arg("--ntasks=1")
        .arg(format!("--cpus-per-task={cpus}"))
        .arg(format!("--mem={mem}"))
        .arg(format!("--time={time}"))
        .arg("--open-mode=append")
        .arg(format!("--output={}", log_pattern(log)?));
    if let Some(account) = &args.account {
        cmd.arg(format!("--account={account}"));
    }
    if let Some(partition) = &args.partition {
        cmd.arg(format!("--partition={partition}"));
    }
    if let Some(dependency) = dependency {
        cmd.arg(format!("--dependency={dependency}"));
        cmd.arg("--kill-on-invalid-dep=yes");
    }
    if let Some(array) = array {
        cmd.arg(format!("--array={array}"));
        // Hold the array until each pending task has its own chunk name.
        cmd.arg("--hold");
    }
    cmd.args(&args.sbatch_arg);
    // The final job writes several outputs and cannot safely resume in place.
    // Index and map tasks can rerun using their atomic completion markers.
    if matches!(stage, JobStage::Finish) {
        cmd.arg("--no-requeue");
    }
    cmd.arg(script);
    let output = cmd
        .output()
        .context("start sbatch; is Slurm available on PATH?")?;
    if !output.status.success() {
        bail!(
            "sbatch failed for {} ({name}): {}",
            stage.description(),
            String::from_utf8_lossy(&output.stderr).trim()
        )
    }
    let value = String::from_utf8(output.stdout).context("invalid sbatch output")?;
    let id = value.trim().split(';').next().unwrap_or_default();
    if id.is_empty() || !id.bytes().all(|b| b.is_ascii_digit()) {
        bail!("unexpected sbatch job ID: {:?}", value.trim())
    }
    Ok(id.to_owned())
}

fn name_array_tasks(array_id: &str, group_id: &str, chunks: usize) -> Result<()> {
    for chunk in 0..chunks {
        let output = Command::new("scontrol")
            .arg("update")
            .arg(format!("JobId={array_id}_{chunk}"))
            .arg(format!("JobName=rust-hic-{group_id}-mapping-{chunk}"))
            .output()
            .context("start scontrol to name Slurm array tasks")?;
        if !output.status.success() {
            bail!(
                "scontrol failed to name chunk {chunk}: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            )
        }
    }
    Ok(())
}

fn release_array(array_id: &str) -> Result<()> {
    let output = Command::new("scontrol")
        .arg("release")
        .arg(array_id)
        .output()
        .context("start scontrol to release the Slurm array")?;
    if !output.status.success() {
        bail!(
            "scontrol failed to release the alignment array: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        )
    }
    Ok(())
}

fn cancel_submitted(jobs: &[String]) -> Result<()> {
    let output = Command::new("scancel")
        .args(jobs)
        .output()
        .context("start scancel")?;
    if !output.status.success() {
        bail!("{}", String::from_utf8_lossy(&output.stderr).trim())
    }
    Ok(())
}

pub(crate) fn submit(
    config: &Path,
    args: SlurmOptions,
    check: bool,
    keep_pairs: bool,
    keep_bam: bool,
) -> Result<()> {
    validate_extra_args(&args.sbatch_arg)?;
    let config_path = fs::canonicalize(config)
        .with_context(|| format!("open configuration {}", config.display()))?;
    let base = config_path.parent().context("missing config parent")?;
    let mut c: Config = serde_yaml::from_reader(File::open(&config_path)?)
        .with_context(|| format!("parse configuration {}", config_path.display()))?;
    c.keep_pairs |= keep_pairs;
    c.keep_bam |= keep_bam;
    let mut c = preflight::prepare_for_submission(c, base)?;
    let fastq = c
        .fastq
        .as_mut()
        .context("Slurm execution requires FASTQ input")?;
    let chunks = estimated_chunks(fastq)?;
    if let Some(n) = args.max_parallel {
        if n == 0 {
            bail!("--max-parallel must be positive")
        }
        fastq.jobs = n;
    }
    fastq.chunks = Some(chunks);
    let parallel = fastq.jobs.min(chunks);
    let (index_cpus, map_cpus, finish_cpus) = job_cpus(&c, &args)?;
    log_pattern(&c.output_dir.join("slurm/map-%A_%a.log"))?;
    let version = Command::new("sbatch")
        .arg("--version")
        .output()
        .context("start sbatch; is Slurm available on PATH?")?;
    if !version.status.success() {
        bail!("sbatch --version failed")
    }
    let version = Command::new("scontrol")
        .arg("--version")
        .output()
        .context("start scontrol; is Slurm available on PATH?")?;
    if !version.status.success() {
        bail!("scontrol --version failed")
    }
    let time = args.time.as_deref().unwrap_or("04:00:00");
    let index_mem = args.index_mem.as_deref().unwrap_or("16G");
    let map_mem = args.map_mem.as_deref().unwrap_or("32G");
    let finish_mem = args.finish_mem.as_deref().unwrap_or("64G");
    if check {
        println!(
            "Slurm preflight passed\n  Output: {}\n  Scratch: {}\n  Chunks: {} (up to {} concurrent)\n  CPUs per job: index {}, map {}, finish {}\n  Memory per job: index {}, map {}, finish {}\n  Reused FASTQ indexes are verified on the compute node; ensure all paths are shared.",
            c.output_dir.display(),
            c.scratch_dir.display(),
            chunks,
            parallel,
            index_cpus,
            map_cpus,
            finish_cpus,
            index_mem,
            map_mem,
            finish_mem
        );
        return Ok(());
    }

    fs::create_dir(&c.output_dir)
        .with_context(|| format!("create output {}", c.output_dir.display()))?;
    let out = fs::canonicalize(&c.output_dir)?;
    let work = tempfile::Builder::new()
        .prefix(WORK_PREFIX)
        .rand_bytes(8)
        .tempdir_in(&c.scratch_dir)?
        .keep();
    let work = fs::canonicalize(&work)?;
    let run_id = work
        .file_name()
        .and_then(OsStr::to_str)
        .context("invalid work directory name")?
        .to_owned();
    let group_id = run_id
        .strip_prefix(WORK_PREFIX)
        .context("invalid work directory prefix")?
        .to_owned();
    if group_id.len() != 8 || !group_id.bytes().all(|b| b.is_ascii_alphanumeric()) {
        bail!("invalid random Slurm group ID")
    }
    fs::write(out.join("work-dir.txt"), format!("{}\n", work.display()))?;
    fs::create_dir(out.join("slurm"))?;
    fs::create_dir(work.join("chunks"))?;
    fs::create_dir(work.join("attempts"))?;
    let binary = work.join("rust-hic");
    fs::copy(env::current_exe()?, &binary)
        .context("copy executable to shared workflow directory")?;
    fs::set_permissions(&binary, fs::Permissions::from_mode(0o755))?;
    let resolved = out.join("config.resolved.yaml");
    serde_yaml::to_writer(File::create(&resolved)?, &c)?;
    let plan = Plan {
        schema: SCHEMA,
        run_id,
        output: out.clone(),
        work: work.clone(),
        config_sha256: hash_file(&resolved)?,
        executable_sha256: hash_file(&binary)?,
        submitted_ms: now_ms()?,
        chunks,
    };
    write_json(&work.join("plan.json"), &plan)?;
    let index_script = script(&work, "index", false)?;
    let map_script = script(&work, "map", true)?;
    let finish_script = script(&work, "finish", false)?;
    let mut jobs = Jobs {
        schema: SCHEMA,
        run_id: plan.run_id.clone(),
        group_id: Some(group_id.clone()),
        chunks,
        index: None,
        alignment: None,
        finish: None,
    };
    let jobs_file = out.join("slurm-jobs.json");
    write_json(&jobs_file, &jobs)?;
    let mut submitted = Vec::new();
    let result = (|| -> Result<()> {
        let id = sbatch(
            &args,
            JobStage::Index,
            &group_id,
            &out.join("slurm/index-%j.log"),
            index_cpus,
            index_mem,
            args.index_time.as_deref().unwrap_or(time),
            None,
            None,
            &index_script,
        )?;
        submitted.push(id.clone());
        jobs.index = Some(id.clone());
        write_json(&jobs_file, &jobs)?;
        let dep = format!("afterok:{id}");
        let array = format!("0-{}%{parallel}", chunks - 1);
        let id = sbatch(
            &args,
            JobStage::Map,
            &group_id,
            &out.join("slurm/map-%A_%a.log"),
            map_cpus,
            map_mem,
            args.map_time.as_deref().unwrap_or(time),
            Some(&dep),
            Some(&array),
            &map_script,
        )?;
        submitted.push(id.clone());
        jobs.alignment = Some(id.clone());
        write_json(&jobs_file, &jobs)?;
        name_array_tasks(&id, &group_id, chunks)?;
        // afterany lets the final job report missing chunks after an array failure.
        let dep = format!("afterany:{id}");
        let id = sbatch(
            &args,
            JobStage::Finish,
            &group_id,
            &out.join("slurm/finish-%j.log"),
            finish_cpus,
            finish_mem,
            args.finish_time.as_deref().unwrap_or(time),
            Some(&dep),
            None,
            &finish_script,
        )?;
        submitted.push(id.clone());
        jobs.finish = Some(id);
        write_json(&jobs_file, &jobs)?;
        release_array(
            jobs.alignment
                .as_deref()
                .context("missing alignment array ID")?,
        )?;
        Ok(())
    })();
    if let Err(e) = result {
        let e = if submitted.is_empty() {
            e
        } else if let Err(cancel_error) = cancel_submitted(&submitted) {
            e.context(format!(
                "unable to cancel submitted Slurm jobs {}; check them manually: {cancel_error:#}",
                submitted.join(",")
            ))
        } else {
            e
        };
        record_failure(&out, "submission", &e);
        return Err(e)
            .with_context(|| format!("Slurm submission failed; inspect {}", out.display()));
    }
    let status_path = shell_quote(out.as_os_str())?;
    println!("Submitted Slurm workflow\n  Group ID: {}\n  Output: {}\n  Index job: {}\n  Alignment array: {} ({} chunks, up to {} concurrent)\n  Finish job: {}\n  Status: rust-hic run status {}",
        group_id, out.display(), jobs.index.as_deref().unwrap_or("?"), jobs.alignment.as_deref().unwrap_or("?"),
        chunks, parallel, jobs.finish.as_deref().unwrap_or("?"), status_path);
    Ok(())
}
