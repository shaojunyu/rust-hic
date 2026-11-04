//! Names and process invocation for internal executable stages.
use std::{io, os::unix::process::CommandExt, process::Command};

macro_rules! define_stages {
    ($( $variant:ident => $name:literal ),+ $(,)?) => {
        #[derive(Clone, Copy)]
        pub(crate) enum Stage { $( $variant ),+ }
        impl Stage {
            pub(crate) fn name(self) -> &'static str {
                match self { $( Self::$variant => $name ),+ }
            }
            pub(crate) fn from_name(name: &str) -> Option<Self> {
                match name { $( $name => Some(Self::$variant), )+ _ => None }
            }
        }
    };
}
define_stages! {
    SortPairs => "sort_pairs",
    MergeRuns => "merge_pair_runs",
    MergePairs => "merge_pairs",
    MergeRanges => "merge_pair_ranges",
    FilterPairs => "filter_pairs",
    DedupPairs => "dedup_pairs",
    PairsStats => "pairs_stats",
    BinPairs => "bin_pairs",
    BuildMcool => "build_mcool",
}

pub(crate) fn stage_command(stage: Stage) -> io::Result<Command> {
    let mut command = Command::new(std::env::current_exe()?);
    command.arg0(stage.name());
    Ok(command)
}
