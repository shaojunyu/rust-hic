//! Dispatch internal stages to their implementations.
use crate::{matrix, pairs, stage::Stage};

pub(crate) fn dispatch(invoked: &std::ffi::OsStr) -> bool {
    let name = std::path::Path::new(invoked)
        .file_name()
        .and_then(|s| s.to_str());
    let Some(stage) = name.and_then(Stage::from_name) else {
        return false;
    };
    match stage {
        Stage::SortPairs => pairs::sort::main(),
        Stage::MergeRuns => pairs::merge_runs::main(),
        Stage::MergePairs => pairs::merge::main(),
        Stage::MergeRanges => pairs::merge_ranges::main(),
        Stage::FilterPairs => pairs::filter::main(),
        Stage::DedupPairs => pairs::dedup::main(),
        Stage::PairsStats => pairs::stats::main(),
        Stage::BinPairs => matrix::bin_pairs::main(),
        Stage::BuildMcool => matrix::multiresolution::main(),
    }
    true
}
