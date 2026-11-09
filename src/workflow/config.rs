//! Configuration shared by the YAML pipeline and command-line workflows.
use super::digest;
use crate::mapping::align;
use anyhow::{bail, Result};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Deserialize, Serialize, Clone)]
#[serde(deny_unknown_fields)]
pub(crate) struct Resources {
    #[serde(default = "six")]
    pub(crate) merge_workers: usize,
    #[serde(default = "six")]
    pub(crate) dedup_workers: usize,
    #[serde(default = "six")]
    pub(crate) matrix_threads: usize,
    #[serde(default = "four")]
    pub(crate) sort_jobs: usize,
    #[serde(default = "two")]
    pub(crate) sort_threads: usize,
    #[serde(default = "four")]
    pub(crate) compression_threads: usize,
    #[serde(default)]
    pub(crate) postprocess_cpus: Option<Vec<usize>>,
    #[serde(default = "run_mib")]
    pub(crate) sort_run_mib: usize,
}
fn six() -> usize {
    6
}
fn four() -> usize {
    4
}
fn two() -> usize {
    2
}
fn run_mib() -> usize {
    256
}
impl Default for Resources {
    fn default() -> Self {
        Self {
            merge_workers: six(),
            dedup_workers: six(),
            matrix_threads: six(),
            sort_jobs: four(),
            sort_threads: two(),
            compression_threads: four(),
            postprocess_cpus: None,
            sort_run_mib: run_mib(),
        }
    }
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Config {
    pub(crate) assembly: String,
    pub(crate) chromsizes: PathBuf,
    pub(crate) fragments: Option<PathBuf>,
    #[serde(default)]
    pub(crate) digest: Option<digest::DigestConfig>,
    #[serde(default)]
    pub(crate) input_pairs: Vec<PathBuf>,
    #[serde(default)]
    pub(crate) fastq: Option<align::FastqConfig>,
    #[serde(default)]
    pub(crate) input_sorted: bool,
    pub(crate) output_dir: PathBuf,
    pub(crate) scratch_dir: PathBuf,
    #[serde(default)]
    pub(crate) merge_scratch_dir: Option<PathBuf>,
    #[serde(default)]
    pub(crate) compressed_pairs: bool,
    #[serde(default)]
    pub(crate) keep_pairs: bool,
    #[serde(default)]
    pub(crate) keep_bam: bool,
    #[serde(default = "dedup_max_mismatch")]
    pub(crate) dedup_max_mismatch: i64,
    #[serde(default)]
    pub(crate) resources: Resources,
    #[serde(default = "resolutions")]
    pub(crate) resolutions: Vec<u64>,
    #[serde(default = "enabled")]
    pub(crate) sort_chromosomes: bool,
    #[serde(default = "iterations")]
    pub(crate) balance_max_iterations: usize,
}
pub(crate) fn dedup_max_mismatch() -> i64 {
    3
}
fn enabled() -> bool {
    true
}
pub(crate) fn resolutions() -> Vec<u64> {
    vec![
        5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000, 5000000,
    ]
}
pub(crate) fn iterations() -> usize {
    1000
}
pub(crate) fn validate_resolutions(resolutions: &[u64]) -> Result<()> {
    if resolutions.is_empty()
        || resolutions.contains(&0)
        || resolutions.iter().min() != resolutions.first()
        || resolutions.iter().any(|x| x % resolutions[0] != 0)
    {
        bail!("resolutions must start with the smallest, and be positive multiples of it")
    }
    let mut unique = resolutions.to_vec();
    unique.sort_unstable();
    unique.dedup();
    if unique.len() != resolutions.len() {
        bail!("duplicate resolution")
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn example_config_matches_defaults() {
        let example: Config = serde_yaml::from_str(include_str!("../../docs/example.config.yaml"))
            .expect("example config must parse");
        let minimal: Config =
            serde_yaml::from_str(include_str!("../../docs/example.minimal.config.yaml"))
                .expect("minimal example config must parse");
        assert_eq!(
            serde_yaml::to_value(example).unwrap(),
            serde_yaml::to_value(minimal).unwrap(),
            "example defaults diverged from the implementation"
        );
    }
}
