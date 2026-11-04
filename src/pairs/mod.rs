pub(crate) mod dedup;
pub(crate) mod filter;
pub(crate) mod format;
pub(crate) mod merge;
pub(crate) mod merge_ranges;
pub(crate) mod merge_runs;
pub(crate) mod sort;
pub(crate) mod stats;

/// Internal byte ranges passed to parallel merge workers.
pub(crate) const MERGE_RANGES_ENV: &str = "RUST_HIC_MERGE_RANGES";

/// Shared ordering for sorted pairs and merge boundaries.
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub(super) struct Key {
    pub(super) c1: String,
    pub(super) c2: String,
    pub(super) p1: u64,
    pub(super) p2: u64,
    pub(super) kind: String,
}
