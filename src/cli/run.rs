//! Unified local and Slurm workflow entry point.
use crate::workflow;
use anyhow::{bail, Context, Result};
use clap::{Args, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Args)]
#[command(
    args_conflicts_with_subcommands = true,
    subcommand_negates_reqs = true,
    after_help = "Create a minimal example with: rust-hic example-config --minimal > config.yaml\nOmit --minimal for the complete example. Slurm mode requires shared FASTQ, reference and scratch paths."
)]
pub(crate) struct RunArgs {
    /// YAML workflow configuration file.
    #[arg(value_name = "CONFIG", required = true)]
    config: Option<PathBuf>,
    /// Run locally or submit a distributed Slurm workflow.
    #[arg(long, value_enum, default_value_t = RunMode::Local)]
    mode: RunMode,
    /// Check inputs and configuration without starting the pipeline.
    #[arg(long)]
    check: bool,
    /// Retain the final sorted, filtered and deduplicated uncompressed pairs file.
    #[arg(long)]
    keep_pairs: bool,
    /// Retain aligner-emitted records in one unsorted, unindexed BAM; requires FASTQ input.
    #[arg(long)]
    keep_bam: bool,
    #[command(flatten, next_help_heading = "Slurm options")]
    slurm: SlurmOptions,
    #[command(subcommand)]
    action: Option<RunAction>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum RunMode {
    Local,
    Slurm,
}

#[derive(Subcommand)]
enum RunAction {
    /// Show the progress and scheduler states of a submitted Slurm run.
    Status {
        /// Output directory of the submitted Slurm workflow.
        output_dir: PathBuf,
    },
    /// Cancel the recorded active jobs of a Slurm run.
    Cancel {
        /// Output directory of the submitted Slurm workflow.
        output_dir: PathBuf,
    },
}

#[derive(Args, Default, PartialEq, Eq)]
pub(crate) struct SlurmOptions {
    /// Slurm account, if the cluster requires one.
    #[arg(long)]
    pub(crate) account: Option<String>,
    /// Slurm partition, if the cluster requires one.
    #[arg(long)]
    pub(crate) partition: Option<String>,
    /// Default time limit for each stage (default: 04:00:00).
    #[arg(long)]
    pub(crate) time: Option<String>,
    /// Override the index job time limit.
    #[arg(long)]
    pub(crate) index_time: Option<String>,
    /// Override the per-chunk alignment time limit.
    #[arg(long)]
    pub(crate) map_time: Option<String>,
    /// Override the final sorting and matrix job time limit.
    #[arg(long)]
    pub(crate) finish_time: Option<String>,
    /// Memory per index job, in Slurm units (default: 16G).
    #[arg(long)]
    pub(crate) index_mem: Option<String>,
    /// Memory per alignment task; each loads its own reference (default: 32G).
    #[arg(long)]
    pub(crate) map_mem: Option<String>,
    /// Memory for sorting, merging, and matrix creation (default: 64G).
    #[arg(long)]
    pub(crate) finish_mem: Option<String>,
    /// Maximum concurrent alignment tasks (defaults to fastq.jobs).
    #[arg(long)]
    pub(crate) max_parallel: Option<usize>,
    /// Override CPUs requested for the index job.
    #[arg(long)]
    pub(crate) index_cpus: Option<usize>,
    /// Override CPUs requested per alignment task.
    #[arg(long)]
    pub(crate) map_cpus: Option<usize>,
    /// Override CPUs requested for the final job.
    #[arg(long)]
    pub(crate) finish_cpus: Option<usize>,
    /// Additional scheduler option in --name=value form; repeatable.
    #[arg(long, allow_hyphen_values = true)]
    pub(crate) sbatch_arg: Vec<String>,
}

pub(crate) fn execute(args: RunArgs) -> Result<()> {
    if let Some(action) = args.action {
        return match action {
            RunAction::Status { output_dir } => workflow::slurm::status(&output_dir),
            RunAction::Cancel { output_dir } => workflow::slurm::cancel(&output_dir),
        };
    }
    let config = args
        .config
        .context("a workflow configuration is required")?;
    match args.mode {
        RunMode::Local => {
            if args.slurm != SlurmOptions::default() {
                bail!("Slurm options require --mode slurm")
            }
            workflow::run(&config, args.check, args.keep_pairs, args.keep_bam)
        }
        RunMode::Slurm => workflow::slurm::submit(
            &config,
            args.slurm,
            args.check,
            args.keep_pairs,
            args.keep_bam,
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Command};
    use clap::Parser;

    #[test]
    fn run_modes_and_management_share_the_run_command() {
        let cli = Cli::try_parse_from(["rust-hic", "run", "config.yaml"]).unwrap();
        let Some(Command::Run(args)) = cli.command else {
            panic!("expected run")
        };
        assert_eq!(args.mode, RunMode::Local);
        assert!(args.action.is_none());
        let cli = Cli::try_parse_from([
            "rust-hic",
            "run",
            "config.yaml",
            "--mode",
            "slurm",
            "--check",
        ])
        .unwrap();
        let Some(Command::Run(args)) = cli.command else {
            panic!("expected run")
        };
        assert_eq!(args.mode, RunMode::Slurm);
        assert!(args.check);
        assert!(Cli::try_parse_from(["rust-hic", "run", "status", "result"]).is_ok());
        assert!(Cli::try_parse_from(["rust-hic", "run", "cancel", "result"]).is_ok());
        assert!(Cli::try_parse_from(["rust-hic", "run"]).is_err());
        assert!(
            Cli::try_parse_from(["rust-hic", "run", "config.yaml", "--mode", "other"]).is_err()
        );
        assert!(
            Cli::try_parse_from(["rust-hic", "run", "config.yaml", "status", "result"]).is_err()
        );
    }
}
