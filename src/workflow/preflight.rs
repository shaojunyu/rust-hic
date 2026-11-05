//! Resolve and validate workflow inputs before starting expensive work.
use super::{
    config::{validate_resolutions, Config},
    digest, path_from,
    process::affinity_mask,
};
use crate::{
    mapping::align,
    pairs::{
        filter::validate_fragments,
        format::{read_header, Format},
    },
};
use anyhow::{bail, Context, Result};
use std::{
    collections::{HashMap, HashSet},
    fs::{self, File},
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
};

fn input_file(base: &Path, path: &Path, label: &str) -> Result<PathBuf> {
    let path = path_from(base, path);
    let resolved =
        fs::canonicalize(&path).with_context(|| format!("open {label} {}", path.display()))?;
    if !resolved.is_file() {
        bail!("{label} must be a regular file: {}", resolved.display());
    }
    File::open(&resolved).with_context(|| format!("read {label} {}", resolved.display()))?;
    Ok(resolved)
}

fn writable_dir(base: &Path, path: &Path, label: &str) -> Result<PathBuf> {
    let path = path_from(base, path);
    let resolved =
        fs::canonicalize(&path).with_context(|| format!("open {label} {}", path.display()))?;
    if !resolved.is_dir() {
        bail!("{label} must be a directory: {}", resolved.display());
    }
    let mut probe = tempfile::tempfile_in(&resolved)
        .with_context(|| format!("{label} is not writable: {}", resolved.display()))?;
    probe.write_all(b"1").with_context(|| {
        format!(
            "cannot write to {label} {}; check permissions and free space",
            resolved.display()
        )
    })?;
    Ok(resolved)
}

/// Check the sizes file against the constraints used by sorting and matrix output.
pub(crate) fn validate_chromsizes(
    path: &Path,
    base_resolution: Option<u64>,
) -> Result<HashMap<String, u64>> {
    let file =
        File::open(path).with_context(|| format!("open chromosome sizes {}", path.display()))?;
    let mut sizes = HashMap::new();
    let mut bins = 0u64;
    for (i, line) in BufReader::new(file).lines().enumerate() {
        let line = line.with_context(|| format!("read {} line {}", path.display(), i + 1))?;
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.is_empty() && base_resolution.is_none() {
            continue;
        }
        if fields.len() != 2 {
            bail!(
                "{}:{}: expected chromosome name and length",
                path.display(),
                i + 1
            );
        }
        let length: u64 = fields[1].parse().with_context(|| {
            format!(
                "{}:{}: invalid chromosome length '{}'",
                path.display(),
                i + 1,
                fields[1]
            )
        })?;
        if length == 0 || sizes.insert(fields[0].to_owned(), length).is_some() {
            bail!(
                "{}:{}: chromosome length must be positive and names must be unique",
                path.display(),
                i + 1
            );
        }
        if let Some(resolution) = base_resolution {
            if length > i32::MAX as u64 || !fields[0].is_ascii() || fields[0].len() > 256 {
                bail!(
                    "{}:{}: chromosome name or length is not supported by matrix output",
                    path.display(),
                    i + 1
                );
            }
            bins = bins
                .checked_add(length.div_ceil(resolution))
                .context("chromosome bin count overflow")?;
        }
    }
    if sizes.is_empty() {
        bail!("chromosome sizes file is empty: {}", path.display());
    }
    if bins > u32::MAX as u64 {
        bail!("chromosome sizes and base resolution exceed the supported number of bins");
    }
    Ok(sizes)
}

/// Compare contigs shared with chromsizes, allowing intentionally different subsets.
fn check_contig(
    name: &str,
    length: u64,
    sizes: &HashMap<String, u64>,
    sizes_path: &Path,
    source: &Path,
    matched: &mut usize,
) -> Result<()> {
    if let Some(expected) = sizes.get(name) {
        *matched += 1;
        if *expected != length {
            bail!(
                "chromosome {name} has length {expected} in {} but {length} in {}",
                sizes_path.display(),
                source.display()
            )
        }
    }
    Ok(())
}

fn check_shared_contigs(matched: usize, sizes_path: &Path, source: &Path) -> Result<()> {
    if matched == 0 {
        bail!(
            "no chromosome names in {} match {}",
            source.display(),
            sizes_path.display()
        )
    }
    Ok(())
}

fn suffixed(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_os_string();
    name.push(suffix);
    name.into()
}

/// BWA's small text .ann file describes the contigs used by its index.
fn check_bwa_contigs(prefix: &Path, sizes_path: &Path, sizes: &HashMap<String, u64>) -> Result<()> {
    let path = suffixed(prefix, ".ann");
    let mut lines = BufReader::new(File::open(&path)?).lines();
    let header = lines.next().context("missing BWA index header")??;
    let count: usize = header
        .split_whitespace()
        .nth(1)
        .context("missing BWA index chromosome count")?
        .parse()
        .context("invalid BWA index chromosome count")?;
    let mut matched = 0;
    for _ in 0..count {
        let name_line = lines
            .next()
            .context("missing BWA index chromosome name")??;
        let name = name_line
            .split_whitespace()
            .nth(1)
            .context("invalid BWA index chromosome name")?;
        let length_line = lines
            .next()
            .context("missing BWA index chromosome length")??;
        let length = length_line
            .split_whitespace()
            .nth(1)
            .context("invalid BWA index chromosome length")?
            .parse()
            .context("invalid BWA index chromosome length")?;
        check_contig(name, length, sizes, sizes_path, &path, &mut matched)?;
    }
    check_shared_contigs(matched, sizes_path, &path)
}

/// Use a FASTA index if one already exists; do not scan or index the FASTA.
fn check_fasta_contigs(
    reference: &Path,
    sizes_path: &Path,
    sizes: &HashMap<String, u64>,
) -> Result<()> {
    let path = suffixed(reference, ".fai");
    let file = match File::open(&path) {
        Ok(file) => file,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(()),
        Err(e) => return Err(e).with_context(|| format!("open FASTA index {}", path.display())),
    };
    let mut matched = 0;
    for line in BufReader::new(file).lines() {
        let line = line?;
        let mut fields = line.split_whitespace();
        let name = fields
            .next()
            .context("missing FASTA index chromosome name")?;
        let length = fields
            .next()
            .context("missing FASTA index chromosome length")?
            .parse()
            .context("invalid FASTA index chromosome length")?;
        check_contig(name, length, sizes, sizes_path, &path, &mut matched)?;
    }
    check_shared_contigs(matched, sizes_path, &path)
}

pub(crate) fn check_pairs(path: &Path) -> Result<(Vec<String>, Format)> {
    let file = File::open(path).with_context(|| format!("open pairs input {}", path.display()))?;
    if !file.metadata()?.is_file() {
        bail!("pairs input must be a regular file: {}", path.display());
    }
    let mut reader = BufReader::new(file);
    if reader.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
        bail!("pairs input must be uncompressed: {}", path.display());
    }
    read_header(&mut reader).with_context(|| format!("read pairs header {}", path.display()))
}

pub(super) fn prepare(mut c: Config, base: &Path, pairs_only: bool) -> Result<Config> {
    if c.assembly.trim().is_empty() || c.assembly.contains(['\r', '\n']) {
        bail!("assembly must be a non-empty single line");
    }
    if c.input_pairs.is_empty() == c.fastq.is_none() {
        bail!("specify exactly one of input_pairs or fastq");
    }
    if c.fastq.is_some() && c.input_sorted {
        bail!("FASTQ input cannot already be sorted");
    }
    if c.keep_bam && c.fastq.is_none() {
        bail!("keep_bam requires FASTQ input; pairs files do not contain the original alignments");
    }
    if c.fragments.is_some() && c.digest.is_some() {
        bail!("specify fragments or digest, not both");
    }
    if c.dedup_max_mismatch < 0 {
        bail!("dedup_max_mismatch must be non-negative");
    }
    let r = &c.resources;
    for (name, value) in [
        ("merge_workers", r.merge_workers),
        ("dedup_workers", r.dedup_workers),
        ("matrix_threads", r.matrix_threads),
        ("compression_threads", r.compression_threads),
    ] {
        if !(1..=64).contains(&value) {
            bail!("resources.{name} must be between 1 and 64; got {value}");
        }
    }
    if r.sort_jobs == 0 || r.sort_threads == 0 {
        bail!("resources.sort_jobs and resources.sort_threads must be positive");
    }
    if r.sort_run_mib == 0 || r.sort_run_mib.checked_mul(1024 * 1024).is_none() {
        bail!("resources.sort_run_mib must be positive and fit in memory-size units");
    }
    if !(1..=10000).contains(&c.balance_max_iterations) {
        bail!("balance_max_iterations must be between 1 and 10000");
    }
    affinity_mask(r.postprocess_cpus.as_deref())?;
    validate_resolutions(&c.resolutions)?;
    c.chromsizes = input_file(base, &c.chromsizes, "chromosome sizes")?;
    let sizes = validate_chromsizes(&c.chromsizes, (!pairs_only).then_some(c.resolutions[0]))?;
    c.fragments = c
        .fragments
        .as_ref()
        .map(|p| input_file(base, p, "restriction fragment BED"))
        .transpose()?;
    if let Some(path) = &c.fragments {
        validate_fragments(path)
            .with_context(|| format!("invalid restriction fragment BED {}", path.display()))?;
    }
    c.input_pairs = c
        .input_pairs
        .iter()
        .map(|p| input_file(base, p, "pairs input"))
        .collect::<Result<_>>()?;
    let mut seen = HashSet::new();
    let mut formats = HashMap::new();
    let mut sorted_header = None;
    for path in &c.input_pairs {
        if !seen.insert(path) {
            bail!("duplicate pairs input: {}", path.display());
        }
        let (header, format) = check_pairs(path)?;
        if c.fragments.is_some() && format == Format::Fragments {
            bail!("restriction filtering requires 8-column pairs; {} already has fragment annotations", path.display());
        }
        formats.insert(format.len(), path);
        if c.input_sorted {
            if let Some(expected) = &sorted_header {
                if expected != &header {
                    bail!("sorted pairs headers do not match: {}", path.display());
                }
            } else {
                sorted_header = Some(header);
            }
        }
    }
    if formats.len() > 1 {
        bail!("pairs inputs use incompatible column layouts");
    }
    if let Some(fastq) = &mut c.fastq {
        align::resolve(fastq, base)?;
        if fastq.aligner != "chromap" {
            check_bwa_contigs(&fastq.aligner_index, &c.chromsizes, &sizes).with_context(|| {
                format!(
                    "check BWA index {}",
                    suffixed(&fastq.aligner_index, ".ann").display()
                )
            })?;
        }
        check_fasta_contigs(&fastq.reference, &c.chromsizes, &sizes).with_context(|| {
            format!(
                "check FASTA index {}",
                suffixed(&fastq.reference, ".fai").display()
            )
        })?;
        if let Some(index) = &fastq.index_dir {
            align::validate_reusable_index(fastq, index)?;
            crate::fastq_index::verify(index, 1)
                .with_context(|| format!("verify reusable FASTQ index {}", index.display()))?;
        }
    }
    if let Some(digest) = &mut c.digest {
        digest::resolve(digest, base)?;
        input_file(base, &digest.reference, "digest reference FASTA")?;
    }
    c.scratch_dir = writable_dir(base, &c.scratch_dir, "scratch directory")?;
    c.merge_scratch_dir = c
        .merge_scratch_dir
        .as_ref()
        .map(|p| writable_dir(base, p, "merge scratch directory"))
        .transpose()?;
    let output = path_from(base, &c.output_dir);
    match fs::symlink_metadata(&output) {
        Ok(_) => bail!("output already exists: {}", output.display()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => return Err(e).with_context(|| format!("check output {}", output.display())),
    }
    let name = output
        .file_name()
        .context("output_dir must include a directory name")?;
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    c.output_dir = writable_dir(base, parent, "output parent directory")?.join(name);
    Ok(c)
}
