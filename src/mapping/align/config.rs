//! FASTQ mapping configuration and validation.
use super::*;

pub(crate) const DEFAULT_INDEX_THREADS_PER_MATE: usize = 8;

#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct FastqConfig {
    pub read1: PathBuf,
    pub read2: PathBuf,
    #[serde(default)]
    pub reference: PathBuf,
    #[serde(default)]
    pub aligner_index: PathBuf,
    #[serde(default)]
    pub aligner_executable: PathBuf,
    #[serde(default = "default_aligner")]
    pub aligner: String,
    #[serde(default)]
    pub aligner_args: Vec<String>,
    #[serde(default)]
    pub chunks: Option<usize>,
    #[serde(default = "default_chunk_size")]
    pub chunk_size_bytes: u64,
    #[serde(default = "four")]
    pub jobs: usize,
    #[serde(default = "eight")]
    pub threads_per_job: usize,
    #[serde(default = "default_index_threads_per_mate")]
    pub index_threads_per_mate: usize,
    #[serde(default = "radius")]
    pub index_search_mib: usize,
    #[serde(default)]
    pub expected_pairs: Option<u64>,
    #[serde(default)]
    pub sam: SamConfig,
    /// Reuse a complete index after checking its inputs and boundaries.
    #[serde(default)]
    pub index_dir: Option<PathBuf>,
}
#[derive(Clone, Deserialize, Serialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct SamConfig {
    pub min_mapq: u8,
    pub max_molecule_size: i64,
    pub max_inter_align_gap: i64,
}
impl Default for SamConfig {
    fn default() -> Self {
        Self {
            min_mapq: 1,
            max_molecule_size: 750,
            max_inter_align_gap: 20,
        }
    }
}
fn default_aligner() -> String {
    "bwa-mem3".into()
}
fn default_chunk_size() -> u64 {
    2 * 1024 * 1024 * 1024
}
fn default_index_threads_per_mate() -> usize {
    DEFAULT_INDEX_THREADS_PER_MATE
}
fn four() -> usize {
    4
}
fn eight() -> usize {
    8
}
fn radius() -> usize {
    128
}
impl FastqConfig {
    pub fn for_index(read1: PathBuf, read2: PathBuf) -> Self {
        Self {
            read1,
            read2,
            reference: PathBuf::new(),
            aligner_index: PathBuf::new(),
            aligner_executable: PathBuf::new(),
            aligner: default_aligner(),
            aligner_args: Vec::new(),
            chunks: None,
            chunk_size_bytes: default_chunk_size(),
            jobs: four(),
            threads_per_job: eight(),
            index_threads_per_mate: default_index_threads_per_mate(),
            index_search_mib: radius(),
            expected_pairs: None,
            sam: SamConfig::default(),
            index_dir: None,
        }
    }
}

pub fn resolve_index(c: &mut FastqConfig, base: &Path) -> Result<()> {
    let fix = |p: &PathBuf, label: &str| {
        let path = if p.is_absolute() {
            p.clone()
        } else {
            base.join(p)
        };
        fs::canonicalize(&path).with_context(|| format!("open {label} {}", path.display()))
    };
    c.read1 = fix(&c.read1, "R1 FASTQ")?;
    c.read2 = fix(&c.read2, "R2 FASTQ")?;
    c.index_dir = c
        .index_dir
        .as_ref()
        .map(|path| fix(path, "FASTQ index directory"))
        .transpose()?;
    if !c.read1.is_file() || !c.read2.is_file() {
        bail!(
            "FASTQ mates must be regular files: {} and {}",
            c.read1.display(),
            c.read2.display()
        )
    }
    use std::os::unix::fs::MetadataExt;
    let r1_meta = fs::metadata(&c.read1)?;
    let r2_meta = fs::metadata(&c.read2)?;
    if (r1_meta.dev(), r1_meta.ino()) == (r2_meta.dev(), r2_meta.ino()) {
        bail!(
            "R1 and R2 refer to the same FASTQ file: {} and {}",
            c.read1.display(),
            c.read2.display()
        )
    }
    for (label, path) in [("R1", &c.read1), ("R2", &c.read2)] {
        let mut file =
            File::open(path).with_context(|| format!("open {label} FASTQ {}", path.display()))?;
        let mut magic = [0u8; 3];
        use std::io::Read;
        file.read_exact(&mut magic)
            .with_context(|| format!("read gzip header from {label} FASTQ {}", path.display()))?;
        if magic != [0x1f, 0x8b, 8] {
            bail!("{label} FASTQ is not gzip-compressed: {}", path.display())
        }
    }
    if c.chunks.is_some_and(|n| !(1..=256).contains(&n)) {
        bail!("FASTQ chunk count must be between 1 and 256")
    }
    if c.chunk_size_bytes == 0 {
        bail!("FASTQ chunk size must be positive")
    }
    if c.jobs == 0 || c.threads_per_job == 0 {
        bail!("FASTQ jobs and threads per job must be positive")
    }
    if c.index_threads_per_mate == 0 || c.index_search_mib == 0 {
        bail!("FASTQ indexing threads per mate and search window must be positive")
    }
    if c.sam.min_mapq == 0 {
        bail!("SAM min_mapq must be positive")
    }
    if !(0..=1_000_000_000).contains(&c.sam.max_molecule_size)
        || !(0..=1_000_000_000).contains(&c.sam.max_inter_align_gap)
    {
        bail!("SAM maximum molecule size and inter-alignment gap must be between 0 and 1000000000")
    }
    if u32::try_from(c.index_threads_per_mate).is_err()
        || u64::try_from(c.index_search_mib)
            .ok()
            .and_then(|n| n.checked_mul(1024 * 1024))
            .is_none()
    {
        bail!("FASTQ indexing thread count or search window is too large")
    }
    if c.index_dir.is_none() && c.chunks.is_none() {
        let largest = fs::metadata(&c.read1)?
            .len()
            .max(fs::metadata(&c.read2)?.len());
        let chunks = largest.div_ceil(c.chunk_size_bytes);
        if !(1..=256).contains(&chunks) {
            bail!("estimated FASTQ chunk count {chunks} is outside 1..=256; increase chunk size or set chunks")
        }
    }
    Ok(())
}

pub fn resolve(c: &mut FastqConfig, base: &Path) -> Result<()> {
    resolve_index(c, base)?;
    if !["chromap", "bwa-mem", "bwa-mem2", "bwa-mem3"].contains(&c.aligner.as_str()) {
        bail!(
            "unsupported aligner '{}'; choose chromap, bwa-mem, bwa-mem2 or bwa-mem3",
            c.aligner
        )
    }
    if c.reference.as_os_str().is_empty() || c.aligner_index.as_os_str().is_empty() {
        bail!("FASTQ alignment requires reference and aligner_index")
    }
    let fix = |p: &PathBuf| {
        fs::canonicalize(if p.is_absolute() {
            p.clone()
        } else {
            base.join(p)
        })
    };
    c.reference = fix(&c.reference)
        .with_context(|| format!("open reference FASTA {}", c.reference.display()))?;
    if !c.reference.is_file() {
        bail!(
            "reference must be a regular file: {}",
            c.reference.display()
        )
    }
    File::open(&c.reference)
        .with_context(|| format!("read reference FASTA {}", c.reference.display()))?;
    c.aligner_index = if c.aligner == "chromap" {
        let index = fix(&c.aligner_index)
            .with_context(|| format!("open Chromap index {}", c.aligner_index.display()))?;
        if !index.is_file() {
            bail!("Chromap index must be a regular file: {}", index.display())
        }
        File::open(&index).with_context(|| format!("read Chromap index {}", index.display()))?;
        index
    } else {
        let p = if c.aligner_index.is_absolute() {
            c.aligner_index.clone()
        } else {
            base.join(&c.aligner_index)
        };
        let parent = fs::canonicalize(p.parent().context("index prefix has no parent")?)?;
        parent.join(p.file_name().context("index prefix has no filename")?)
    };
    if c.aligner == "bwa-mem3" {
        #[cfg(target_arch = "x86_64")]
        if !std::is_x86_feature_detected!("sse4.2") {
            bail!("built-in bwa-mem3 requires an x86-64 CPU with SSE4.2")
        }
        if !c.aligner_executable.as_os_str().is_empty() {
            bail!("--aligner-executable does not apply to the built-in bwa-mem3")
        }
        super::super::bwa_mem3::options(c)?;
    }
    for arg in &c.aligner_args {
        if c.aligner == "bwa-mem3" {
            continue;
        }
        let conflicts = if c.aligner == "chromap" {
            [
                "--preset",
                "--low-mem",
                "--index",
                "--reference",
                "--threads",
                "--output",
            ]
            .iter()
            .any(|flag| arg == flag || arg.starts_with(&format!("{flag}=")))
                || ["-x", "-r", "-t", "-1", "-2", "-o"]
                    .iter()
                    .any(|flag| arg.starts_with(flag))
        } else {
            ["-t", "-o", "-p", "-S", "-P"]
                .iter()
                .any(|flag| arg.starts_with(flag))
        };
        if arg == "--" || conflicts {
            bail!("aligner argument '{arg}' conflicts with options managed by rust-hic")
        }
    }
    if c.aligner != "chromap" {
        let suffixes: &[&str] = if c.aligner == "bwa-mem" {
            &[".amb", ".ann", ".bwt", ".pac", ".sa"]
        } else if c.aligner == "bwa-mem2" {
            &[".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"]
        } else {
            &[".amb", ".ann", ".bwt.2bit.64", ".pac"]
        };
        for suffix in suffixes {
            let mut name = c.aligner_index.as_os_str().to_os_string();
            name.push(suffix);
            let path = PathBuf::from(name);
            if !path.is_file() {
                bail!("missing {} index file: {}", c.aligner, path.display())
            }
            File::open(&path)
                .with_context(|| format!("read {} index file {}", c.aligner, path.display()))?;
        }
    }
    if c.aligner == "bwa-mem3" {
        return Ok(());
    }
    let executable = if c.aligner_executable.as_os_str().is_empty() {
        PathBuf::from(match c.aligner.as_str() {
            "bwa-mem" => "bwa",
            "bwa-mem2" => "bwa-mem2",
            _ => "chromap",
        })
    } else {
        c.aligner_executable.clone()
    };
    let local = if executable.is_absolute() {
        executable.clone()
    } else {
        base.join(&executable)
    };
    let executable = if local.is_file() {
        fs::canonicalize(local)?
    } else if executable.components().count() == 1 {
        std::env::split_paths(&std::env::var_os("PATH").unwrap_or_default())
            .map(|dir| dir.join(&executable))
            .find(|path| {
                fs::metadata(path)
                    .map(|m| m.is_file() && m.permissions().mode() & 0o111 != 0)
                    .unwrap_or(false)
            })
            .with_context(|| {
                format!(
                    "aligner executable '{}' not found on PATH",
                    executable.display()
                )
            })?
            .canonicalize()?
    } else {
        bail!("aligner executable not found: {}", local.display())
    };
    if fs::metadata(&executable)?.permissions().mode() & 0o111 == 0 {
        bail!(
            "aligner executable is not executable: {}",
            executable.display()
        )
    }
    c.aligner_executable = executable;
    if c.aligner == "chromap" && c.sam != SamConfig::default() {
        bail!("SAM parsing options apply only to BWA aligners")
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_threads_default_and_override() {
        let default: FastqConfig = serde_yaml::from_str("read1: r1\nread2: r2\n").unwrap();
        let overridden: FastqConfig =
            serde_yaml::from_str("read1: r1\nread2: r2\nindex_threads_per_mate: 1\n").unwrap();
        assert_eq!(default.index_threads_per_mate, 8);
        assert_eq!(overridden.index_threads_per_mate, 1);
        assert_eq!(
            FastqConfig::for_index("r1".into(), "r2".into()).index_threads_per_mate,
            8
        );
    }
}
