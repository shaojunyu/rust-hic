//! Public command-line arguments and entry points.
pub(crate) mod index;
pub(crate) mod map;
pub(crate) mod postprocess;

use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "rust-hic", version, about = "Hi-C processing pipeline")]
pub(crate) struct Cli {
    #[command(subcommand)]
    pub command: Option<Command>,
}

#[derive(Subcommand)]
pub(crate) enum Command {
    /// Print an example YAML configuration to standard output.
    ExampleConfig {
        /// Print only the required settings; omit to show all available settings.
        #[arg(long)]
        minimal: bool,
    },
    /// Generate an HTML and JSON report from a completed run.
    Report {
        /// Directory containing a completed run.
        result_dir: PathBuf,
        /// Existing output directory; defaults to the run directory.
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Re-read final pairs to calculate detailed contact distributions.
        #[arg(long)]
        recompute_contact_stats: bool,
    },
    /// Run the pipeline from a YAML configuration.
    #[command(
        after_help = "Create a minimal example with: rust-hic example-config --minimal > config.yaml\nOmit --minimal for the complete example."
    )]
    Run {
        /// YAML configuration file.
        config: PathBuf,
        /// Check inputs and configuration without running the pipeline.
        #[arg(long)]
        check: bool,
        /// Retain the final sorted, filtered and deduplicated pairs file.
        #[arg(long)]
        keep_pairs: bool,
        /// Retain all aligner-emitted records in one unsorted BAM file.
        #[arg(long)]
        keep_bam: bool,
    },
    /// Build an index for paired gzip FASTQ files.
    #[command(
        after_help = "Boundaries are aligned on read IDs. Chunk size accepts bytes, MB, GB, MiB and GiB. Relative paths use the current directory."
    )]
    Index(IndexArgs),
    /// Map indexed paired FASTQ reads to one unsorted pairs file.
    #[command(
        after_help = "Built-in bwa-mem3 uses a bwa-mem2 index and supports -k, -w, -T, -A, -B and -c through --aligner-arg. Failed runs retain temporary files and logs."
    )]
    Map(MapArgs),
    /// Sort, filter and deduplicate pairs.
    #[command(
        after_help = "Without --fragments, restriction filtering is skipped. --threads is per stage, not a global CPU limit."
    )]
    Clean(CleanArgs),
    /// Build a multi-resolution matrix from cleaned pairs.
    Matrix(MatrixArgs),
    /// Show version information.
    Version,
}

#[derive(Args)]
pub(crate) struct IndexArgs {
    /// Input read 1 gzip FASTQ.
    pub read1: PathBuf,
    /// Input read 2 gzip FASTQ.
    pub read2: PathBuf,
    /// Output directory; existing index files are never overwritten.
    #[arg(short, long)]
    pub output: PathBuf,
    /// Split into N chunks (1–256); conflicts with --chunk-size.
    #[arg(long, conflicts_with = "chunk_size")]
    pub chunks: Option<usize>,
    /// Target compressed chunk size, e.g. 2GiB (default 2GiB).
    #[arg(long)]
    pub chunk_size: Option<String>,
    /// Indexing threads per input file.
    #[arg(long, default_value_t = crate::mapping::align::DEFAULT_INDEX_THREADS_PER_MATE)]
    pub threads_per_mate: usize,
    /// Search window for matching chunk boundaries in MiB.
    #[arg(long, default_value_t = 128)]
    pub search_mib: usize,
}

#[derive(Args)]
pub(crate) struct MapArgs {
    /// Paired FASTQ index directory.
    pub index_dir: PathBuf,
    /// Reference FASTA.
    #[arg(short = 'r', long)]
    pub reference: PathBuf,
    /// Aligner index path.
    #[arg(short = 'x', long)]
    pub aligner_index: PathBuf,
    /// Output pairs file; all selected chunks are combined.
    #[arg(short, long)]
    pub output: PathBuf,
    /// Aligner: chromap, bwa-mem, bwa-mem2, or built-in bwa-mem3.
    #[arg(long, default_value = "bwa-mem3")]
    pub aligner: String,
    /// Path to an external aligner executable.
    #[arg(long)]
    pub aligner_executable: Option<PathBuf>,
    /// Additional shell-quoted aligner arguments; repeatable.
    #[arg(long, allow_hyphen_values = true)]
    pub aligner_arg: Vec<String>,
    /// Map only one chunk, numbered from zero.
    #[arg(long)]
    pub chunk: Option<usize>,
    /// Number of concurrent chunks.
    #[arg(long, default_value_t = 4)]
    pub jobs: usize,
    /// Aligner threads per chunk.
    #[arg(long, default_value_t = 8)]
    pub threads_per_job: usize,
    /// Minimum MAPQ for BWA unique classification.
    #[arg(long, default_value_t = 1)]
    pub min_mapq: u8,
    /// Maximum BWA molecule size in bases.
    #[arg(long, default_value_t = 750)]
    pub max_molecule_size: i64,
    /// Maximum BWA inter-alignment gap in bases.
    #[arg(long, default_value_t = 20)]
    pub max_inter_align_gap: i64,
}

#[derive(Args)]
pub(crate) struct CleanArgs {
    /// Input pairs files.
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,
    /// Chromosome names and lengths.
    #[arg(long)]
    pub chromsizes: PathBuf,
    /// Output sorted pairs file.
    #[arg(short, long)]
    pub output: PathBuf,
    /// Restriction fragments for filtering.
    #[arg(long)]
    pub fragments: Option<PathBuf>,
    /// Genome assembly name.
    #[arg(long, default_value = "unknown")]
    pub assembly: String,
    /// Threads per processing stage.
    #[arg(long, default_value_t = 4)]
    pub threads: usize,
    /// Maximum distance on each end in bp.
    #[arg(long, default_value_t = 3)]
    pub max_mismatch: i64,
}

#[derive(Args)]
pub(crate) struct MatrixArgs {
    /// Cleaned pairs file.
    pub input: PathBuf,
    /// Chromosome names and lengths.
    #[arg(long)]
    pub chromsizes: PathBuf,
    /// Output mcool file.
    #[arg(short, long)]
    pub output: PathBuf,
    /// Genome assembly name.
    #[arg(long, default_value = "unknown")]
    pub assembly: String,
    /// Comma-separated bin sizes, smallest first.
    #[arg(long)]
    pub resolutions: Option<String>,
    /// Matrix worker threads.
    #[arg(long, default_value_t = 4)]
    pub threads: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accepts_repeatable_hyphenated_aligner_args() {
        let cli = Cli::try_parse_from([
            "rust-hic",
            "map",
            "idx",
            "-r",
            "ref",
            "-x",
            "ref",
            "-o",
            "out",
            "--aligner-arg",
            "-k 19 -T 30",
            "--aligner-arg=-w 50",
        ])
        .unwrap();
        let Some(Command::Map(args)) = cli.command else {
            panic!("expected map")
        };
        assert_eq!(args.aligner_arg, ["-k 19 -T 30", "-w 50"]);
    }

    #[test]
    fn rejects_duplicate_and_conflicting_options() {
        assert!(Cli::try_parse_from([
            "rust-hic",
            "index",
            "r1",
            "r2",
            "-o",
            "out",
            "--chunks",
            "2",
            "--chunk-size",
            "1GiB"
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "rust-hic",
            "clean",
            "input",
            "-o",
            "out",
            "--chromsizes",
            "sizes",
            "--threads",
            "1",
            "--threads",
            "2"
        ])
        .is_err());
    }

    #[test]
    fn index_and_clean_defaults_are_stable() {
        let cli = Cli::try_parse_from(["rust-hic", "index", "r1", "r2", "-o", "out"]).unwrap();
        let Some(Command::Index(args)) = cli.command else {
            panic!("expected index")
        };
        assert_eq!((args.threads_per_mate, args.search_mib), (8, 128));
        let cli = Cli::try_parse_from([
            "rust-hic",
            "clean",
            "--chromsizes",
            "sizes",
            "-o",
            "out",
            "--",
            "-input",
        ])
        .unwrap();
        let Some(Command::Clean(args)) = cli.command else {
            panic!("expected clean")
        };
        assert_eq!((args.threads, args.max_mismatch), (4, 3));
        assert_eq!(args.inputs, [PathBuf::from("-input")]);
    }
}
