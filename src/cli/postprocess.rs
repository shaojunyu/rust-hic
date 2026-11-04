//! Command-line entry points for cleaning pairs and building matrices.
use super::{CleanArgs, MatrixArgs};
use crate::workflow::config::{Config, Resources};
use anyhow::{bail, Context, Result};
use std::{
    env, fs,
    path::{Path, PathBuf},
};

pub fn clean(args: CleanArgs) -> Result<()> {
    run(
        args.inputs,
        args.chromsizes,
        args.output,
        args.assembly,
        args.threads,
        args.fragments,
        args.max_mismatch,
        None,
        true,
    )
}
pub fn matrix(args: MatrixArgs) -> Result<()> {
    run(
        vec![args.input],
        args.chromsizes,
        args.output,
        args.assembly,
        args.threads,
        None,
        3,
        args.resolutions,
        false,
    )
}
#[allow(clippy::too_many_arguments)]
fn run(
    inputs: Vec<PathBuf>,
    chromsizes: PathBuf,
    output: PathBuf,
    assembly: String,
    threads: usize,
    fragments: Option<PathBuf>,
    max_mismatch: i64,
    resolution_arg: Option<String>,
    pairs_only: bool,
) -> Result<()> {
    let command = if pairs_only { "clean" } else { "matrix" };
    if !(1..=64).contains(&threads) {
        bail!("--threads must be between 1 and 64")
    }
    if max_mismatch < 0 {
        bail!("--max-mismatch must be non-negative")
    }
    let resolutions = if let Some(v) = resolution_arg {
        v.split(',')
            .map(|s| {
                s.parse::<u64>()
                    .context("invalid --resolutions; use comma-separated integers")
            })
            .collect::<Result<Vec<_>>>()?
    } else {
        crate::workflow::config::resolutions()
    };
    crate::workflow::config::validate_resolutions(&resolutions)?;
    if assembly.is_empty() {
        bail!("--assembly cannot be empty")
    }
    if output.file_name().is_none() {
        bail!("--output must include a file name")
    }
    match fs::symlink_metadata(&output) {
        Ok(_) => bail!("output already exists: {}", output.display()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => return Err(e).context("check output path"),
    }
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    if !pairs_only {
        crate::workflow::preflight::validate_chromsizes(&chromsizes, Some(resolutions[0]))?;
        crate::workflow::preflight::check_pairs(&inputs[0])?;
    }
    let workspace = tempfile::Builder::new()
        .prefix(if pairs_only {
            ".rust-hic-clean-"
        } else {
            ".rust-hic-matrix-"
        })
        .tempdir_in(parent)
        .with_context(|| format!("create {command} workspace next to output"))?;
    let result = if pairs_only {
        let jobs = threads.min(inputs.len());
        let resources = Resources {
            merge_workers: threads,
            dedup_workers: threads,
            matrix_threads: threads,
            sort_jobs: jobs,
            sort_threads: (threads / jobs).max(1),
            ..Resources::default()
        };
        let c = Config {
            assembly,
            chromsizes,
            fragments,
            digest: None,
            input_pairs: inputs,
            fastq: None,
            input_sorted: false,
            output_dir: workspace.path().join("result"),
            scratch_dir: workspace.path().to_path_buf(),
            merge_scratch_dir: None,
            compressed_pairs: false,
            dedup_max_mismatch: max_mismatch,
            resources,
            resolutions,
            balance_max_iterations: crate::workflow::config::iterations(),
        };
        if c.fragments.is_none() {
            eprintln!("Restriction-fragment filtering skipped (no --fragments BED supplied)");
        }
        crate::workflow::run_config(c, &env::current_dir()?, false, true)
    } else {
        build_matrix(
            &inputs[0],
            &chromsizes,
            &assembly,
            &resolutions,
            threads,
            workspace.path(),
        )
    };
    if let Err(error) = result {
        let path = workspace.keep();
        return Err(error).with_context(|| {
            format!(
                "{command} failed; logs and intermediate files retained in {}",
                path.display()
            )
        });
    }
    let staged = workspace.path().join(if pairs_only {
        "result/filtered-nodups.pairs"
    } else {
        "result/contacts.mcool"
    });
    if let Err(error) = fs::hard_link(&staged, &output) {
        let path = workspace.keep();
        return Err(error).with_context(|| {
            format!(
                "save {command} output to {}; completed result retained in {}",
                output.display(),
                path.display()
            )
        });
    }
    eprintln!(
        "{} written to {}",
        if pairs_only { "Pairs" } else { "Matrix" },
        output.display()
    );
    Ok(())
}

fn build_matrix(
    input: &Path,
    chromsizes: &Path,
    assembly: &str,
    resolutions: &[u64],
    threads: usize,
    work: &Path,
) -> Result<()> {
    let input = fs::canonicalize(input).context("open cleaned pairs")?;
    let chromsizes = fs::canonicalize(chromsizes).context("open chromosome sizes")?;
    let logs = work.join("result");
    fs::create_dir(&logs)?;
    let started = std::time::Instant::now();
    crate::workflow::log(&logs, "cool_start", started)?;
    let base = work.join("base.cool");
    let mut first = crate::stage::stage_command(crate::stage::Stage::BinPairs)?;
    first
        .arg(&chromsizes)
        .arg(resolutions[0].to_string())
        .arg(assembly)
        .arg(&input)
        .arg(&base)
        .arg(work)
        .arg(threads.to_string());
    crate::workflow::process::exec(first, &logs, "pairs-to-cool", None)?;
    crate::workflow::log(&logs, "cool_done", started)?;
    crate::workflow::log(&logs, "mcool_start", started)?;
    let mut zoom = crate::stage::stage_command(crate::stage::Stage::BuildMcool)?;
    zoom.arg(&base)
        .arg(&chromsizes)
        .arg(logs.join("contacts.mcool"))
        .arg(assembly)
        .arg(
            resolutions
                .iter()
                .map(u64::to_string)
                .collect::<Vec<_>>()
                .join(","),
        )
        .arg(threads.to_string())
        .env(
            crate::matrix::BALANCE_MAX_ITERS_ENV,
            crate::workflow::config::iterations().to_string(),
        );
    crate::workflow::process::exec(zoom, &logs, "cool-to-mcool", None)?;
    crate::workflow::log(&logs, "mcool_done", started)?;
    Ok(())
}
