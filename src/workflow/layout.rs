//! Resolve and expose every path written to a workflow result directory.
//!
//! A result directory is always grouped:
//!
//! ```text
//! <root>/matrix/<prefix>.mcool
//! <root>/pairs/<prefix>.filtered-nodups.pairs[.gz[.px2]]
//! <root>/alignments/<prefix>.bam
//! <root>/qc/<prefix>.report.html|json
//! <root>/qc/<prefix>.contact-stats.json
//! <root>/qc/<prefix>.filter-stats.json
//! <root>/logs/events.jsonl and stage logs
//! <root>/status.json, failure.json, config.resolved.yaml, slurm-jobs.json,
//!        work-dir.txt, MANIFEST.tsv
//! ```
//!
//! Older results used a flat directory. `for_existing` recognises those so that
//! `run report` keeps working, but new runs never write the flat layout.

use super::config::Config;
use anyhow::{bail, Context, Result};
use std::{
    fs,
    path::{Path, PathBuf},
};

pub(crate) const MATRIX_DIR: &str = "matrix";
pub(crate) const PAIRS_DIR: &str = "pairs";
pub(crate) const ALIGNMENTS_DIR: &str = "alignments";
pub(crate) const QC_DIR: &str = "qc";
pub(crate) const LOGS_DIR: &str = "logs";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Grouped,
    Flat,
}

#[derive(Clone, Debug)]
pub(crate) struct Layout {
    root: PathBuf,
    prefix: String,
    kind: Kind,
}

impl Layout {
    /// Grouped layout for a new run. `prefix` must already be resolved.
    pub(crate) fn new(root: &Path, prefix: &str) -> Result<Self> {
        validate_prefix(prefix)?;
        Ok(Self {
            root: root.to_path_buf(),
            prefix: prefix.to_owned(),
            kind: Kind::Grouped,
        })
    }

    /// Flat layout, used only to read results produced before grouping existed.
    fn flat(root: &Path) -> Self {
        Self {
            root: root.to_path_buf(),
            prefix: String::new(),
            kind: Kind::Flat,
        }
    }

    pub(crate) fn from_config(root: &Path, c: &Config) -> Result<Self> {
        let prefix = match &c.output.prefix {
            Some(prefix) => prefix.clone(),
            None => derive_prefix(c)?,
        };
        Self::new(root, &prefix)
    }

    /// Pick the grouped layout if this run wrote one, otherwise the legacy flat one.
    pub(crate) fn for_existing(root: &Path, c: &Config) -> Result<Self> {
        // Check legacy names first so old runs do not need a derivable prefix.
        if root.join("contacts.mcool").is_file()
            || root.join("report.json").is_file()
            || root.join("contact-stats.json").is_file()
        {
            return Ok(Self::flat(root));
        }
        Self::from_config(root, c)
    }

    /// Fill `output.prefix`, deriving it from the input filenames when omitted.
    pub(crate) fn resolve_config(c: &mut Config) -> Result<()> {
        if let Some(prefix) = &c.output.prefix {
            return validate_prefix(prefix);
        }
        let prefix = derive_prefix(c)?;
        validate_prefix(&prefix)?;
        c.output.prefix = Some(prefix);
        Ok(())
    }

    pub(crate) fn root(&self) -> &Path {
        &self.root
    }

    pub(crate) fn prefix(&self) -> &str {
        &self.prefix
    }

    pub(crate) fn is_grouped(&self) -> bool {
        self.kind == Kind::Grouped
    }

    fn in_dir(&self, dir: &str, name: String) -> PathBuf {
        match self.kind {
            Kind::Grouped => self.root.join(dir).join(name),
            Kind::Flat => self.root.join(name),
        }
    }

    pub(crate) fn logs(&self) -> PathBuf {
        match self.kind {
            Kind::Grouped => self.root.join(LOGS_DIR),
            Kind::Flat => self.root.clone(),
        }
    }

    pub(crate) fn slurm_logs(&self) -> PathBuf {
        self.logs().join("slurm")
    }

    pub(crate) fn log(&self, name: &str) -> PathBuf {
        self.logs().join(format!("{name}.log"))
    }

    pub(crate) fn events(&self) -> PathBuf {
        self.logs().join("events.jsonl")
    }

    pub(crate) fn mcool(&self) -> PathBuf {
        match self.kind {
            Kind::Grouped => self.in_dir(MATRIX_DIR, format!("{}.mcool", self.prefix)),
            Kind::Flat => self.root.join("contacts.mcool"),
        }
    }

    pub(crate) fn pairs(&self) -> PathBuf {
        match self.kind {
            Kind::Grouped => {
                self.in_dir(PAIRS_DIR, format!("{}.filtered-nodups.pairs", self.prefix))
            }
            Kind::Flat => self.root.join("filtered-nodups.pairs"),
        }
    }

    pub(crate) fn compressed_pairs(&self) -> PathBuf {
        match self.kind {
            Kind::Grouped => self.in_dir(
                PAIRS_DIR,
                format!("{}.filtered-nodups.pairs.gz", self.prefix),
            ),
            Kind::Flat => self.root.join("filtered-nodups.pairs.gz"),
        }
    }

    pub(crate) fn pairix(&self) -> PathBuf {
        match self.kind {
            Kind::Grouped => self.in_dir(
                PAIRS_DIR,
                format!("{}.filtered-nodups.pairs.gz.px2", self.prefix),
            ),
            Kind::Flat => self.root.join("filtered-nodups.pairs.gz.px2"),
        }
    }

    pub(crate) fn bam(&self) -> PathBuf {
        match self.kind {
            Kind::Grouped => self.in_dir(ALIGNMENTS_DIR, format!("{}.bam", self.prefix)),
            Kind::Flat => self.root.join("alignments.bam"),
        }
    }

    pub(crate) fn report_html(&self) -> PathBuf {
        match self.kind {
            Kind::Grouped => self.in_dir(QC_DIR, format!("{}.report.html", self.prefix)),
            Kind::Flat => self.root.join("report.html"),
        }
    }

    pub(crate) fn report_json(&self) -> PathBuf {
        match self.kind {
            Kind::Grouped => self.in_dir(QC_DIR, format!("{}.report.json", self.prefix)),
            Kind::Flat => self.root.join("report.json"),
        }
    }

    pub(crate) fn contact_stats(&self) -> PathBuf {
        match self.kind {
            Kind::Grouped => self.in_dir(QC_DIR, format!("{}.contact-stats.json", self.prefix)),
            Kind::Flat => self.root.join("contact-stats.json"),
        }
    }

    pub(crate) fn filter_stats(&self) -> PathBuf {
        match self.kind {
            Kind::Grouped => self.in_dir(QC_DIR, format!("{}.filter-stats.json", self.prefix)),
            Kind::Flat => self.root.join("filter-stats.json"),
        }
    }

    pub(crate) fn config_resolved(&self) -> PathBuf {
        self.root.join("config.resolved.yaml")
    }

    pub(crate) fn status(&self) -> PathBuf {
        self.root.join("status.json")
    }

    pub(crate) fn failure(&self) -> PathBuf {
        self.root.join("failure.json")
    }

    pub(crate) fn work_dir_file(&self) -> PathBuf {
        self.root.join("work-dir.txt")
    }

    pub(crate) fn slurm_jobs(&self) -> PathBuf {
        self.root.join("slurm-jobs.json")
    }

    pub(crate) fn manifest(&self) -> PathBuf {
        self.root.join("MANIFEST.tsv")
    }

    /// Create the grouped subdirectories. Safe to call once a new root exists.
    pub(crate) fn create_dirs(&self, keep_bam: bool) -> Result<()> {
        if !self.is_grouped() {
            return Ok(());
        }
        for dir in [MATRIX_DIR, PAIRS_DIR, QC_DIR, LOGS_DIR] {
            let path = self.root.join(dir);
            fs::create_dir_all(&path)
                .with_context(|| format!("create output directory {}", path.display()))?;
        }
        if keep_bam {
            let path = self.root.join(ALIGNMENTS_DIR);
            fs::create_dir_all(&path)
                .with_context(|| format!("create output directory {}", path.display()))?;
        }
        Ok(())
    }

    pub(crate) fn create_slurm_dirs(&self) -> Result<()> {
        if !self.is_grouped() {
            return Ok(());
        }
        let path = self.slurm_logs();
        fs::create_dir_all(&path)
            .with_context(|| format!("create Slurm log directory {}", path.display()))?;
        Ok(())
    }

    /// Write a size/path listing for every file under the result root.
    pub(crate) fn write_manifest(&self) -> Result<()> {
        if !self.is_grouped() {
            return Ok(());
        }
        let mut entries = Vec::new();
        collect_files(&self.root, &self.root, &mut entries)?;
        entries.sort_by(|a, b| a.1.cmp(&b.1));
        let path = self.manifest();
        let mut file = fs::File::create(&path)
            .with_context(|| format!("create manifest {}", path.display()))?;
        use std::io::Write;
        writeln!(file, "size_bytes\tpath")?;
        for (size, name) in entries {
            writeln!(file, "{size}\t{name}")?;
        }
        Ok(())
    }
}

fn collect_files(dir: &Path, root: &Path, out: &mut Vec<(u64, String)>) -> Result<()> {
    for entry in fs::read_dir(dir).with_context(|| format!("read {}", dir.display()))? {
        let entry = entry?;
        let path = entry.path();
        let metadata = fs::symlink_metadata(&path)?;
        if metadata.is_dir() {
            collect_files(&path, root, out)?;
        } else if metadata.is_file() {
            let relative = path
                .strip_prefix(root)
                .with_context(|| format!("manifest path {}", path.display()))?;
            if relative == Path::new("MANIFEST.tsv") {
                continue;
            }
            let name = relative
                .to_str()
                .with_context(|| format!("manifest path is not UTF-8: {}", path.display()))?;
            out.push((metadata.len(), name.to_owned()));
        }
    }
    Ok(())
}

/// Validate an explicit or derived sample prefix used as a filename component.
pub(crate) fn validate_prefix(prefix: &str) -> Result<()> {
    if prefix.is_empty() {
        bail!("output.prefix must not be empty");
    }
    if prefix.len() > 128 {
        bail!("output.prefix must be at most 128 characters");
    }
    if prefix.starts_with('.') || prefix.contains("..") || prefix.ends_with('.') {
        bail!("output.prefix must not start or end with '.' or contain '..'");
    }
    if !prefix
        .bytes()
        .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'.' | b'_' | b'-'))
    {
        bail!("output.prefix may only contain ASCII letters, digits, '.', '_' and '-'");
    }
    Ok(())
}

fn derive_prefix(c: &Config) -> Result<String> {
    if let Some(fastq) = &c.fastq {
        return derive_from_fastq(&fastq.read1);
    }
    let Some(first) = c.input_pairs.first() else {
        bail!("cannot derive output prefix: set output.prefix or provide an input");
    };
    if c.input_pairs.len() > 1 {
        eprintln!(
            "warning: deriving output prefix from {}; set output.prefix to override",
            first.display()
        );
    }
    derive_from_pairs(first)
}

fn file_name(path: &Path) -> Result<&str> {
    path.file_name()
        .and_then(|name| name.to_str())
        .with_context(|| format!("input filename is not valid UTF-8: {}", path.display()))
}

fn strip_suffix_ci<'a>(value: &'a str, suffix: &str) -> Option<&'a str> {
    let start = value.len().checked_sub(suffix.len())?;
    let tail = value.get(start..)?;
    tail.eq_ignore_ascii_case(suffix).then(|| &value[..start])
}

fn strip_fastq_suffix(name: &str) -> &str {
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"] {
        if let Some(stem) = strip_suffix_ci(name, suffix) {
            if !stem.is_empty() {
                return stem;
            }
        }
    }
    name
}

fn strip_mate_suffix(mut name: &str) -> &str {
    // Strip repeatedly so Illumina's _R1_001 form loses both tokens.
    loop {
        let mut changed = false;
        for suffix in [
            "_001", "_002", "_r1", "_r2", "-r1", "-r2", ".r1", ".r2", "_1", "_2", "-1", "-2", ".1",
            ".2",
        ] {
            if let Some(stem) = strip_suffix_ci(name, suffix) {
                if !stem.is_empty() {
                    name = stem;
                    changed = true;
                    break;
                }
            }
        }
        if !changed {
            return name;
        }
    }
}

fn derive_from_fastq(path: &Path) -> Result<String> {
    let name = file_name(path)?;
    let stem = strip_mate_suffix(strip_fastq_suffix(name))
        .trim_end_matches(['.', '_', '-'])
        .to_owned();
    if stem.is_empty() {
        bail!(
            "cannot derive output prefix from {}; set output.prefix",
            path.display()
        );
    }
    Ok(stem)
}

fn derive_from_pairs(path: &Path) -> Result<String> {
    let name = file_name(path)?;
    let name = strip_suffix_ci(name, ".gz").unwrap_or(name);
    let stem = strip_suffix_ci(name, ".pairs").unwrap_or(name);
    let stem = stem.trim_end_matches(['.', '_', '-']).to_owned();
    if stem.is_empty() {
        bail!(
            "cannot derive output prefix from {}; set output.prefix",
            path.display()
        );
    }
    Ok(stem)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config(yaml: &str) -> Config {
        serde_yaml::from_str(yaml).expect("test configuration must parse")
    }

    fn fastq(read1: &str) -> Config {
        let yaml = format!(
            "assembly: hg38\nchromsizes: s\noutput_dir: out\nscratch_dir: .\nfastq:\n  read1: {read1}\n  read2: R2.fastq.gz\n"
        );
        config(&yaml)
    }

    fn pairs(input: &str) -> Config {
        let yaml = format!(
            "assembly: hg38\nchromsizes: s\noutput_dir: out\nscratch_dir: .\ninput_pairs: [{input}]\n"
        );
        config(&yaml)
    }

    #[test]
    fn derives_prefix_from_common_fastq_names() {
        for (name, expected) in [
            ("GM12878_100M_R1.fastq.gz", "GM12878_100M"),
            ("GM12878_full_R1.fastq.gz", "GM12878_full"),
            ("sample_R1_001.fq.gz", "sample"),
            ("sample.1.fastq.gz", "sample"),
            ("sample-1.fastq", "sample"),
            ("R1.fastq.gz", "R1"),
        ] {
            let mut c = fastq(name);
            Layout::resolve_config(&mut c).unwrap();
            assert_eq!(c.output.prefix.as_deref(), Some(expected), "{name}");
        }
    }

    #[test]
    fn derives_prefix_from_pairs_names() {
        let mut c = pairs("sample1.pairs");
        Layout::resolve_config(&mut c).unwrap();
        assert_eq!(c.output.prefix.as_deref(), Some("sample1"));
    }

    #[test]
    fn explicit_prefix_wins() {
        let mut c = fastq("GM12878_100M_R1.fastq.gz");
        c.output.prefix = Some("tumor".into());
        Layout::resolve_config(&mut c).unwrap();
        assert_eq!(c.output.prefix.as_deref(), Some("tumor"));
    }

    #[test]
    fn rejects_unsafe_prefixes() {
        for prefix in [
            "", "a/b", "a\\b", "..", ".hidden", "a..b", "a b", "a:b", "a.", "ü",
        ] {
            assert!(
                validate_prefix(prefix).is_err(),
                "{prefix:?} must be rejected"
            );
        }
    }

    #[test]
    fn grouped_paths_use_prefix_and_directories() {
        let layout = Layout::new(Path::new("/r"), "sampleA").unwrap();
        assert_eq!(layout.mcool(), Path::new("/r/matrix/sampleA.mcool"));
        assert_eq!(
            layout.pairs(),
            Path::new("/r/pairs/sampleA.filtered-nodups.pairs")
        );
        assert_eq!(
            layout.compressed_pairs(),
            Path::new("/r/pairs/sampleA.filtered-nodups.pairs.gz")
        );
        assert_eq!(
            layout.pairix(),
            Path::new("/r/pairs/sampleA.filtered-nodups.pairs.gz.px2")
        );
        assert_eq!(layout.bam(), Path::new("/r/alignments/sampleA.bam"));
        assert_eq!(layout.report_html(), Path::new("/r/qc/sampleA.report.html"));
        assert_eq!(layout.report_json(), Path::new("/r/qc/sampleA.report.json"));
        assert_eq!(
            layout.contact_stats(),
            Path::new("/r/qc/sampleA.contact-stats.json")
        );
        assert_eq!(
            layout.filter_stats(),
            Path::new("/r/qc/sampleA.filter-stats.json")
        );
        assert_eq!(layout.events(), Path::new("/r/logs/events.jsonl"));
        assert_eq!(layout.log("dedup"), Path::new("/r/logs/dedup.log"));
        assert_eq!(layout.status(), Path::new("/r/status.json"));
        assert_eq!(layout.manifest(), Path::new("/r/MANIFEST.tsv"));
    }

    #[test]
    fn flat_layout_keeps_legacy_names() {
        let layout = Layout::flat(Path::new("/r"));
        assert_eq!(layout.mcool(), Path::new("/r/contacts.mcool"));
        assert_eq!(layout.pairs(), Path::new("/r/filtered-nodups.pairs"));
        assert_eq!(layout.bam(), Path::new("/r/alignments.bam"));
        assert_eq!(layout.report_json(), Path::new("/r/report.json"));
        assert_eq!(layout.events(), Path::new("/r/events.jsonl"));
        assert_eq!(layout.log("dedup"), Path::new("/r/dedup.log"));
    }

    #[test]
    fn existing_flat_result_is_detected() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("contacts.mcool"), b"").unwrap();
        let c = pairs("sample1.pairs");
        let layout = Layout::for_existing(dir.path(), &c).unwrap();
        assert_eq!(layout.mcool(), dir.path().join("contacts.mcool"));
        assert_eq!(layout.report_json(), dir.path().join("report.json"));
    }

    #[test]
    fn existing_grouped_result_is_detected() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::create_dir(dir.path().join("matrix")).unwrap();
        std::fs::write(dir.path().join("matrix/sample1.mcool"), b"").unwrap();
        let c = pairs("sample1.pairs");
        let layout = Layout::for_existing(dir.path(), &c).unwrap();
        assert_eq!(layout.mcool(), dir.path().join("matrix/sample1.mcool"));
        assert_eq!(
            layout.report_json(),
            dir.path().join("qc/sample1.report.json")
        );
    }
}
