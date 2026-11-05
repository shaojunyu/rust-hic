//! Read workflow progress and cancel its recorded Slurm jobs.
use super::state::{read_failure, read_json, Jobs, SCHEMA};
use anyhow::{bail, Context, Result};
use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

fn read_jobs(output_dir: &Path) -> Result<(PathBuf, Jobs)> {
    let out = fs::canonicalize(output_dir)
        .with_context(|| format!("open output {}", output_dir.display()))?;
    let jobs: Jobs = read_json(&out.join("slurm-jobs.json"))?;
    if jobs.schema != SCHEMA
        || !(1..=256).contains(&jobs.chunks)
        || jobs
            .group_id
            .as_ref()
            .is_some_and(|id| id.len() != 8 || !id.bytes().all(|b| b.is_ascii_alphanumeric()))
        || [&jobs.index, &jobs.alignment, &jobs.finish]
            .into_iter()
            .flatten()
            .any(|id| id.is_empty() || !id.bytes().all(|b| b.is_ascii_digit()))
    {
        bail!("invalid Slurm job record")
    }
    Ok((out, jobs))
}

fn job_ids(jobs: &Jobs) -> Vec<&str> {
    [&jobs.index, &jobs.alignment, &jobs.finish]
        .into_iter()
        .flatten()
        .map(String::as_str)
        .collect()
}

pub(crate) fn status(output_dir: &Path) -> Result<()> {
    let (out, jobs) = read_jobs(output_dir)?;
    println!("Output: {}", out.display());
    if let Some(group_id) = &jobs.group_id {
        println!("Group ID: {group_id}");
    }
    let complete = out.join("status.json").exists();
    let failure = if complete { None } else { read_failure(&out)? };
    let failed = failure.is_some();
    if complete {
        println!("Workflow: complete");
    } else if let Some(failure) = failure {
        println!(
            "Workflow: failed ({})",
            failure
                .get("phase")
                .and_then(|v| v.as_str())
                .unwrap_or("unknown stage")
        );
        if let Some(error) = failure.get("error").and_then(|v| v.as_str()) {
            println!("  {error}");
        }
    }
    if let Ok(value) = fs::read_to_string(out.join("work-dir.txt")) {
        let work = Path::new(value.trim());
        let complete = (0..jobs.chunks)
            .filter(|i| {
                work.join("chunks")
                    .join(format!("{i:04}/complete.json"))
                    .is_file()
            })
            .count();
        println!("Alignment chunks: {complete}/{}", jobs.chunks);
    }
    let ids = job_ids(&jobs);
    if !ids.is_empty() {
        println!("Slurm jobs: {}", ids.join(","));
        match Command::new("squeue")
            .arg("--noheader")
            .arg(format!("--jobs={}", ids.join(",")))
            .arg("--format=%i %j %T %M %R")
            .output()
        {
            Ok(output) if output.status.success() => {
                if !output.stdout.is_empty() {
                    if !complete && !failed {
                        println!("Workflow: active or waiting");
                    }
                    print!("{}", String::from_utf8_lossy(&output.stdout));
                } else if !complete && !failed {
                    println!("Workflow: not complete; no active jobs found. Inspect Slurm accounting and logs.");
                }
            }
            Ok(output) => {
                if !complete && !failed {
                    println!("Workflow: status unknown");
                }
                eprintln!(
                    "Unable to query Slurm: {}",
                    String::from_utf8_lossy(&output.stderr).trim()
                )
            }
            Err(e) => {
                if !complete && !failed {
                    println!("Workflow: status unknown");
                }
                eprintln!("Unable to query Slurm: {e}")
            }
        }
    } else if !complete && !failed {
        println!("Workflow: not complete; no job IDs recorded.");
    }
    Ok(())
}

pub(crate) fn cancel(output_dir: &Path) -> Result<()> {
    let (out, jobs) = read_jobs(output_dir)?;
    let ids = job_ids(&jobs);
    if ids.is_empty() {
        println!("No submitted jobs recorded in {}", out.display());
        return Ok(());
    }
    // Query first so a successfully completed index job does not cause scancel
    // to report an invalid ID while later jobs are still active.
    let output = Command::new("squeue")
        .arg("--noheader")
        .arg(format!("--jobs={}", ids.join(",")))
        .arg("--format=%A")
        .output()
        .context("query Slurm jobs before cancellation")?;
    if !output.status.success() {
        bail!(
            "squeue failed: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        )
    }
    let text = String::from_utf8(output.stdout).context("invalid squeue output")?;
    let mut active: Vec<_> = text
        .lines()
        .map(str::trim)
        .filter(|s| ids.contains(s))
        .collect();
    active.sort_unstable();
    active.dedup();
    if active.is_empty() {
        println!("No active Slurm jobs for {}", out.display());
        return Ok(());
    }
    let output = Command::new("scancel")
        .args(&active)
        .output()
        .context("start scancel")?;
    if !output.status.success() {
        bail!(
            "scancel failed: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        )
    }
    println!(
        "Cancelled Slurm jobs: {}\nOutput and temporary files remain in place.",
        active.join(",")
    );
    Ok(())
}
