//! Self-contained reports for completed workflow results.
mod html;

use super::config::Config;
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeMap,
    fs::{self, File},
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
    process::Stdio,
};

#[derive(Serialize)]
struct Report {
    schema_version: u32,
    report_generator_version: &'static str,
    assembly: String,
    input: Input,
    restriction: Restriction,
    dedup_max_mismatch: i64,
    counts: Counts,
    timings: Vec<Timing>,
    elapsed_seconds: f64,
    chunks: Vec<Chunk>,
    resolutions: Vec<Resolution>,
    outputs: Outputs,
    notes: Vec<String>,
}

#[derive(Serialize)]
struct Input {
    mode: &'static str,
    aligner: Option<String>,
    paths: Vec<PathBuf>,
    reference: Option<PathBuf>,
    chromsizes: PathBuf,
    fastq_index_reused: Option<bool>,
}

#[derive(Serialize)]
struct Restriction {
    mode: &'static str,
    source: Option<PathBuf>,
    stats: Option<FilterStats>,
}

#[derive(Clone, Deserialize, Serialize)]
struct FilterStats {
    input_records: u64,
    retained_records: u64,
    same_fragment: u64,
    dangling: u64,
    self_ligation: u64,
    unknown_mechanism: u64,
}

#[derive(Serialize)]
struct Counts {
    input_read_pairs: Option<u64>,
    pairs_before_restriction: u64,
    pairs_before_dedup: u64,
    duplicates: u64,
    final_contacts: u64,
    cis: u64,
    trans: u64,
    cis_20kb_plus: u64,
    cis_shorter_than_20kb: u64,
    pair_types: Option<BTreeMap<String, u64>>,
    strand_orientations: Option<BTreeMap<String, u64>>,
    cis_distance_bins: Option<Vec<DistanceBin>>,
}

#[derive(Deserialize, Serialize)]
struct DistanceBin {
    start_bp: u64,
    end_bp: Option<u64>,
    count: u64,
}

#[derive(Deserialize)]
struct ContactStats {
    total_nodups: u64,
    cis: u64,
    trans: u64,
    cis_20kb_plus: u64,
    cis_shorter_than_20kb: u64,
    #[serde(default)]
    pair_types: Option<BTreeMap<String, u64>>,
    #[serde(default)]
    strand_orientations: Option<BTreeMap<String, u64>>,
    #[serde(default)]
    cis_distance_bins: Option<Vec<DistanceBin>>,
}

#[derive(Deserialize, Serialize)]
struct Chunk {
    chunk: usize,
    input_pairs: u64,
    #[serde(default)]
    selected_pairs: Option<u64>,
    elapsed_seconds: f64,
}

#[derive(Serialize)]
struct Timing {
    stage: &'static str,
    seconds: f64,
}

#[derive(Serialize)]
struct Resolution {
    bp: u64,
    bins: usize,
    pixels: usize,
    usable_bins: usize,
    converged: Option<bool>,
    iterations: Option<f64>,
    variance: Option<f64>,
}

#[derive(Serialize)]
struct Outputs {
    result_directory: PathBuf,
    pairs: Option<&'static str>,
    bam: Option<&'static str>,
    mcool: &'static str,
    compressed_pairs: Option<&'static str>,
    pairix_index: Option<&'static str>,
}

fn read_json<T: for<'de> Deserialize<'de>>(path: &Path) -> Result<T> {
    serde_json::from_reader(File::open(path).with_context(|| format!("open {}", path.display()))?)
        .with_context(|| format!("parse {}", path.display()))
}

fn dedup_counts(path: &Path) -> Result<(u64, u64, u64)> {
    let reader =
        BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
    let mut result = None;
    for line in reader.lines() {
        let line = line?;
        if !line.starts_with("event=dedup_done ") {
            continue;
        }
        let fields: BTreeMap<_, _> = line
            .split_whitespace()
            .filter_map(|s| s.split_once('='))
            .collect();
        let get = |name: &str| -> Result<u64> {
            fields
                .get(name)
                .with_context(|| format!("missing {name} in {}", path.display()))?
                .parse()
                .with_context(|| format!("invalid {name} in {}", path.display()))
        };
        if result.is_some() {
            bail!("duplicate dedup summary in {}", path.display())
        }
        result = Some((
            get("input_records")?,
            get("duplicates")?,
            get("output_records")?,
        ));
    }
    result.with_context(|| format!("missing dedup summary in {}", path.display()))
}

fn chunk_stats(out: &Path) -> Result<Vec<Chunk>> {
    let mut chunks = Vec::new();
    for entry in fs::read_dir(out)? {
        let entry = entry?;
        let name = entry.file_name();
        let name = name.to_string_lossy();
        if !name.starts_with("align-")
            || !name.ends_with("-stats.json")
            || name.ends_with("-sam-stats.json")
        {
            continue;
        }
        let mut chunk: Chunk = read_json(&entry.path())?;
        if chunk.selected_pairs.is_none() {
            let sam = out.join(format!("align-{:04}-sam-stats.json", chunk.chunk));
            if sam.exists() {
                let stats: serde_json::Value = read_json(&sam)?;
                chunk.selected_pairs = stats.get("selected_pairs").and_then(|v| v.as_u64());
            }
        }
        chunks.push(chunk);
    }
    chunks.sort_by_key(|x| x.chunk);
    for pair in chunks.windows(2) {
        if pair[0].chunk == pair[1].chunk {
            bail!("duplicate alignment chunk {}", pair[0].chunk)
        }
    }
    Ok(chunks)
}

fn timings(out: &Path, fastq: bool) -> Result<Vec<Timing>> {
    let reader = BufReader::new(File::open(out.join("events.jsonl"))?);
    let mut events = BTreeMap::new();
    for line in reader.lines() {
        let record: serde_json::Value =
            serde_json::from_str(&line?).context("parse pipeline events")?;
        if record.get("chunk").is_some() {
            continue;
        }
        if let (Some(name), Some(seconds)) = (
            record.get("event").and_then(|v| v.as_str()),
            record.get("elapsed_seconds").and_then(|v| v.as_f64()),
        ) {
            events.insert(name.to_owned(), seconds);
        }
    }
    let mut result = Vec::new();
    let stages = [
        ("Restriction digest", "digest_start", "digest_done"),
        ("FASTQ index", "fastq_index_start", "fastq_index_done"),
        (
            if fastq {
                "Alignment (overlaps chunk preparation)"
            } else {
                "Sort pairs"
            },
            if fastq {
                "alignment_start"
            } else {
                "sort_start"
            },
            if fastq { "alignment_done" } else { "sort_done" },
        ),
        (
            "Alignment and chunk preparation",
            "alignment_start",
            "sort_done",
        ),
        (
            "Merge and deduplicate",
            "merge_filter_dedup_start",
            "merge_filter_dedup_done",
        ),
        ("Create base Cooler", "cool_start", "cool_done"),
        ("Create and balance mcool", "mcool_start", "mcool_done"),
        ("Compress and index pairs", "pairix_start", "pairix_done"),
    ];
    for (stage, start, end) in stages {
        if let (Some(a), Some(b)) = (events.get(start), events.get(end)) {
            if b < a {
                bail!("pipeline event {end} precedes {start}")
            }
            result.push(Timing {
                stage,
                seconds: b - a,
            });
        }
    }
    Ok(result)
}

fn finite_attr(dataset: &hdf5::Dataset, name: &str) -> Result<Option<f64>> {
    let value: f64 = dataset.attr(name)?.read_scalar()?;
    Ok(value.is_finite().then_some(value))
}

fn matrix_resolutions(path: &Path, resolutions: &[u64]) -> Result<Vec<Resolution>> {
    let file = hdf5::File::open(path).with_context(|| format!("open matrix {}", path.display()))?;
    let mut result = Vec::new();
    for &bp in resolutions {
        let group = file.group(&format!("resolutions/{bp}"))?;
        let weights = group.dataset("bins/weight")?;
        let values: Vec<f64> = weights.read_raw()?;
        let usable_bins = values
            .iter()
            .filter(|value| value.is_finite() && **value > 0.0)
            .count();
        let raw_converged: bool = weights.attr("converged")?.read_scalar()?;
        result.push(Resolution {
            bp,
            bins: values.len(),
            pixels: group.dataset("pixels/count")?.size(),
            usable_bins,
            converged: (usable_bins > 0).then_some(raw_converged),
            iterations: finite_attr(&weights, "iterations")?,
            variance: finite_attr(&weights, "var")?,
        });
    }
    Ok(result)
}

fn collect(out: &Path, c: &Config, elapsed_seconds: f64, contact_path: &Path) -> Result<Report> {
    let contact: ContactStats = read_json(contact_path)?;
    let (before_dedup, duplicates, after_dedup) = dedup_counts(&out.join("dedup.log"))?;
    if before_dedup.checked_sub(duplicates) != Some(after_dedup)
        || after_dedup != contact.total_nodups
        || contact.cis.checked_add(contact.trans) != Some(after_dedup)
        || contact
            .cis_20kb_plus
            .checked_add(contact.cis_shorter_than_20kb)
            != Some(contact.cis)
    {
        bail!("inconsistent deduplication or contact statistics")
    }
    for (name, counts) in [
        ("pair type", contact.pair_types.as_ref()),
        ("strand orientation", contact.strand_orientations.as_ref()),
    ] {
        if let Some(counts) = counts {
            let total = counts
                .values()
                .try_fold(0u64, |sum, value| sum.checked_add(*value));
            if total != Some(after_dedup) {
                bail!("inconsistent {name} statistics")
            }
        }
    }
    if let Some(bins) = &contact.cis_distance_bins {
        let total = bins
            .iter()
            .try_fold(0u64, |sum, bin| sum.checked_add(bin.count));
        if total != Some(contact.cis)
            || bins.first().map(|bin| bin.start_bp) != Some(0)
            || bins.last().and_then(|bin| bin.end_bp).is_some()
            || bins.windows(2).any(|pair| {
                pair[0].end_bp != Some(pair[1].start_bp) || pair[0].start_bp >= pair[1].start_bp
            })
        {
            bail!("inconsistent cis-distance statistics")
        }
    }
    let filter_path = out.join("filter-stats.json");
    let filter: Option<FilterStats> = if filter_path.exists() {
        Some(read_json(&filter_path)?)
    } else {
        None
    };
    if let Some(f) = &filter {
        if f.retained_records != before_dedup
            || f.retained_records.checked_add(f.same_fragment) != Some(f.input_records)
            || f.dangling
                .checked_add(f.self_ligation)
                .and_then(|n| n.checked_add(f.unknown_mechanism))
                != Some(f.same_fragment)
        {
            bail!("inconsistent restriction-filter statistics")
        }
    }
    let enabled = c.digest.is_some() || c.fragments.is_some();
    if enabled != filter.is_some() {
        bail!("restriction-filter statistics do not match the configuration")
    }
    let before_filter = filter.as_ref().map_or(before_dedup, |f| f.input_records);
    let chunks = chunk_stats(out)?;
    let read_pairs = if c.fastq.is_some() {
        if chunks.is_empty() || chunks.iter().enumerate().any(|(i, c)| c.chunk != i) {
            bail!("missing or non-contiguous alignment chunk statistics")
        }
        if let Some(n) = c.fastq.as_ref().and_then(|f| f.chunks) {
            if chunks.len() != n {
                bail!("alignment chunk statistics do not match configured chunk count")
            }
        }
        Some(
            chunks
                .iter()
                .try_fold(0u64, |a, x| a.checked_add(x.input_pairs))
                .context("input pair count overflow")?,
        )
    } else {
        None
    };
    if c.fastq.is_some() && chunks.iter().all(|c| c.selected_pairs.is_some()) {
        let selected = chunks
            .iter()
            .try_fold(0u64, |a, c| a.checked_add(c.selected_pairs.unwrap()))
            .context("selected pair count overflow")?;
        if selected != before_filter {
            bail!("alignment and downstream pair counts do not match")
        }
    }
    let fastq = c.fastq.as_ref();
    let input = Input {
        mode: if fastq.is_some() { "fastq" } else { "pairs" },
        aligner: fastq.map(|f| f.aligner.clone()),
        paths: fastq.map_or_else(
            || c.input_pairs.clone(),
            |f| vec![f.read1.clone(), f.read2.clone()],
        ),
        reference: fastq.map(|f| f.reference.clone()),
        chromsizes: c.chromsizes.clone(),
        fastq_index_reused: fastq.map(|f| f.index_dir.is_some()),
    };
    let restriction = Restriction {
        mode: if c.digest.is_some() {
            "arima_digest"
        } else if c.fragments.is_some() {
            "bed"
        } else {
            "disabled"
        },
        source: c
            .digest
            .as_ref()
            .map(|d| d.reference.clone())
            .or_else(|| c.fragments.clone()),
        stats: filter,
    };
    let resolutions = matrix_resolutions(&out.join("contacts.mcool"), &c.resolutions)?;
    let mut notes = vec!["Pairs entering downstream processing are not a standard mapping-rate measure; aligner-specific selection and filtering may already have occurred.".into(), "Timings for alignment and chunk preparation overlap and must not be added together.".into()];
    if resolutions.iter().any(|r| r.converged != Some(true)) {
        notes.push("At least one resolution did not converge or has no usable bins; inspect its balancing weights before using balanced values.".into());
    }
    Ok(Report {
        schema_version: 1,
        report_generator_version: env!("CARGO_PKG_VERSION"),
        assembly: c.assembly.clone(),
        input,
        restriction,
        dedup_max_mismatch: c.dedup_max_mismatch,
        counts: Counts {
            input_read_pairs: read_pairs,
            pairs_before_restriction: before_filter,
            pairs_before_dedup: before_dedup,
            duplicates,
            final_contacts: after_dedup,
            cis: contact.cis,
            trans: contact.trans,
            cis_20kb_plus: contact.cis_20kb_plus,
            cis_shorter_than_20kb: contact.cis_shorter_than_20kb,
            pair_types: contact.pair_types,
            strand_orientations: contact.strand_orientations,
            cis_distance_bins: contact.cis_distance_bins,
        },
        timings: timings(out, fastq.is_some())?,
        elapsed_seconds,
        chunks,
        resolutions,
        outputs: Outputs {
            result_directory: out.to_path_buf(),
            pairs: out
                .join("filtered-nodups.pairs")
                .is_file()
                .then_some("filtered-nodups.pairs"),
            bam: out
                .join("alignments.bam")
                .is_file()
                .then_some("alignments.bam"),
            mcool: "contacts.mcool",
            compressed_pairs: c.compressed_pairs.then_some("filtered-nodups.pairs.gz"),
            pairix_index: c.compressed_pairs.then_some("filtered-nodups.pairs.gz.px2"),
        },
        notes,
    })
}

fn write_report(
    out: &Path,
    destination: &Path,
    c: &Config,
    elapsed_seconds: f64,
    contact_path: &Path,
) -> Result<()> {
    if !destination.is_dir() {
        bail!(
            "report destination must be an existing directory: {}",
            destination.display()
        )
    }
    for name in ["report.html", "report.json"] {
        if destination.join(name).exists() {
            bail!(
                "report output already exists: {}",
                destination.join(name).display()
            )
        }
    }
    let report = collect(out, c, elapsed_seconds, contact_path)?;
    let json = serde_json::to_vec_pretty(&report)?;
    let page = html::render(&report);
    let mut json_file = tempfile::NamedTempFile::new_in(destination)?;
    json_file.write_all(&json)?;
    json_file.write_all(b"\n")?;
    json_file.flush()?;
    let mut html_file = tempfile::NamedTempFile::new_in(destination)?;
    html_file.write_all(page.as_bytes())?;
    html_file.flush()?;
    json_file.persist_noclobber(destination.join("report.json"))?;
    html_file.persist_noclobber(destination.join("report.html"))?;
    Ok(())
}

pub(super) fn generate(out: &Path, c: &Config, elapsed_seconds: f64) -> Result<()> {
    write_report(
        out,
        out,
        c,
        elapsed_seconds,
        &out.join("contact-stats.json"),
    )
    .context("generate workflow report")
}

pub(crate) fn regenerate(out: &Path, destination: Option<&Path>, recompute: bool) -> Result<()> {
    let out =
        fs::canonicalize(out).with_context(|| format!("open workflow result {}", out.display()))?;
    let status: serde_json::Value = read_json(&out.join("status.json"))?;
    if status.get("status").and_then(|v| v.as_str()) != Some("complete") {
        bail!("workflow result is not complete")
    }
    let elapsed = status
        .get("elapsed_seconds")
        .and_then(|v| v.as_f64())
        .context("missing elapsed time in status.json")?;
    let c: Config = serde_yaml::from_reader(File::open(out.join("config.resolved.yaml"))?)
        .context("parse resolved configuration")?;
    let destination = destination.unwrap_or(&out);
    let temporary = if recompute {
        if !out.join("filtered-nodups.pairs").is_file() {
            bail!("cannot recompute contact statistics: final pairs were not retained; rerun with --keep-pairs")
        }
        if !destination.is_dir() {
            bail!(
                "report destination must be an existing directory: {}",
                destination.display()
            )
        }
        for name in ["report.html", "report.json"] {
            if destination.join(name).exists() {
                bail!(
                    "report output already exists: {}",
                    destination.join(name).display()
                )
            }
        }
        eprintln!(
            "Recomputing contact statistics from final pairs; this reads the entire pairs file."
        );
        let temp = tempfile::tempdir_in(destination)?;
        let status = crate::stage::stage_command(crate::stage::Stage::PairsStats)?
            .arg(temp.path().join("contact-stats.json"))
            .stdin(File::open(out.join("filtered-nodups.pairs"))?)
            .stdout(Stdio::null())
            .status()
            .context("recompute contact statistics")?;
        if !status.success() {
            bail!("contact statistics recomputation failed: {status}")
        }
        Some(temp)
    } else {
        None
    };
    let stats_path = temporary
        .as_ref()
        .map(|dir| dir.path().join("contact-stats.json"))
        .unwrap_or_else(|| out.join("contact-stats.json"));
    write_report(&out, destination, &c, elapsed, &stats_path)?;
    println!(
        "Report written to {}",
        destination.join("report.html").display()
    );
    Ok(())
}
