//! Self-contained reports for completed workflow results.
mod html;

use super::{config::Config, layout::Layout};
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
    pairs: Option<String>,
    bam: Option<String>,
    mcool: String,
    compressed_pairs: Option<String>,
    pairix_index: Option<String>,
    report_html: String,
    report_json: String,
}

fn relative(root: &Path, path: &Path) -> String {
    path.strip_prefix(root)
        .unwrap_or(path)
        .display()
        .to_string()
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

fn chunk_stats(layout: &Layout) -> Result<Vec<Chunk>> {
    let mut chunks = Vec::new();
    for entry in fs::read_dir(layout.logs())? {
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
            let sam = layout
                .logs()
                .join(format!("align-{:04}-sam-stats.json", chunk.chunk));
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

fn timings(layout: &Layout, fastq: bool) -> Result<Vec<Timing>> {
    let distributed = layout.slurm_jobs().is_file();
    let reader = BufReader::new(File::open(layout.events())?);
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
                if distributed {
                    "Alignment (Slurm array)"
                } else {
                    "Alignment (overlaps chunk preparation)"
                }
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
            if fastq && distributed {
                "Chunk preparation (Slurm finish job)"
            } else {
                ""
            },
            "sort_start",
            "sort_done",
        ),
        (
            if distributed {
                "Alignment through chunk preparation (includes waits)"
            } else {
                "Alignment and chunk preparation"
            },
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
        if stage.is_empty() {
            continue;
        }
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

fn collect(
    layout: &Layout,
    c: &Config,
    elapsed_seconds: f64,
    contact_path: &Path,
) -> Result<Report> {
    let contact: ContactStats = read_json(contact_path)?;
    let (before_dedup, duplicates, after_dedup) = dedup_counts(&layout.log("dedup"))?;
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
    let filter_path = layout.filter_stats();
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
    let chunks = chunk_stats(layout)?;
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
    let resolutions = matrix_resolutions(&layout.mcool(), &c.resolutions)?;
    let mut notes = vec!["Pairs entering downstream processing are not a standard mapping-rate measure; aligner-specific selection and filtering may already have occurred.".into()];
    if layout.slurm_jobs().is_file() {
        notes.push("The Slurm alignment interval spans the first task start through the last task finish and can include scheduler waiting between tasks. The interval through chunk preparation also includes waiting for the finish job and any restriction digest. Total elapsed time includes queue waits.".into());
    } else {
        notes.push(
            "Timings for alignment and chunk preparation overlap and must not be added together."
                .into(),
        );
    }
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
        timings: timings(layout, fastq.is_some())?,
        elapsed_seconds,
        chunks,
        resolutions,
        outputs: Outputs {
            result_directory: layout.root().to_path_buf(),
            pairs: layout
                .pairs()
                .is_file()
                .then(|| relative(layout.root(), &layout.pairs())),
            bam: layout
                .bam()
                .is_file()
                .then(|| relative(layout.root(), &layout.bam())),
            mcool: relative(layout.root(), &layout.mcool()),
            compressed_pairs: c
                .compressed_pairs
                .then(|| relative(layout.root(), &layout.compressed_pairs())),
            pairix_index: c
                .compressed_pairs
                .then(|| relative(layout.root(), &layout.pairix())),
            report_html: relative(layout.root(), &layout.report_html()),
            report_json: relative(layout.root(), &layout.report_json()),
        },
        notes,
    })
}

#[allow(clippy::too_many_arguments)]
fn write_report(
    layout: &Layout,
    destination: &Path,
    html_path: &Path,
    json_path: &Path,
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
    for path in [html_path, json_path] {
        if path.exists() {
            bail!("report output already exists: {}", path.display())
        }
    }
    let report = collect(layout, c, elapsed_seconds, contact_path)?;
    let json = serde_json::to_vec_pretty(&report)?;
    let page = html::render(&report);
    let mut json_file = tempfile::NamedTempFile::new_in(destination)?;
    json_file.write_all(&json)?;
    json_file.write_all(b"\n")?;
    json_file.flush()?;
    let mut html_file = tempfile::NamedTempFile::new_in(destination)?;
    html_file.write_all(page.as_bytes())?;
    html_file.flush()?;
    json_file.persist_noclobber(json_path)?;
    html_file.persist_noclobber(html_path)?;
    Ok(())
}

pub(super) fn generate(layout: &Layout, c: &Config, elapsed_seconds: f64) -> Result<()> {
    write_report(
        layout,
        layout.root(),
        &layout.report_html(),
        &layout.report_json(),
        c,
        elapsed_seconds,
        &layout.contact_stats(),
    )
    .context("generate workflow report")
}

pub(crate) fn regenerate(out: &Path, destination: Option<&Path>, recompute: bool) -> Result<()> {
    let out =
        fs::canonicalize(out).with_context(|| format!("open workflow result {}", out.display()))?;
    // Control files always live at the result root, including for older runs.
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
    let layout = Layout::for_existing(&out, &c)?;

    let destination = destination.unwrap_or(&out);
    if !destination.is_dir() {
        bail!(
            "report destination must be an existing directory: {}",
            destination.display()
        )
    }
    let (html_path, json_path) = if destination == out.as_path() {
        (layout.report_html(), layout.report_json())
    } else {
        (
            destination.join("report.html"),
            destination.join("report.json"),
        )
    };

    let temporary = if recompute {
        if !layout.pairs().is_file() {
            bail!("cannot recompute contact statistics: final pairs were not retained; rerun with --keep-pairs")
        }
        for path in [&html_path, &json_path] {
            if path.exists() {
                bail!("report output already exists: {}", path.display())
            }
        }
        eprintln!(
            "Recomputing contact statistics from final pairs; this reads the entire pairs file."
        );
        let temp = tempfile::tempdir_in(destination)?;
        let status = crate::stage::stage_command(crate::stage::Stage::PairsStats)?
            .arg(temp.path().join("contact-stats.json"))
            .stdin(File::open(layout.pairs())?)
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
        .unwrap_or_else(|| layout.contact_stats());
    write_report(
        &layout,
        destination,
        &html_path,
        &json_path,
        &c,
        elapsed,
        &stats_path,
    )?;
    println!("Report written to {}", html_path.display());
    Ok(())
}
