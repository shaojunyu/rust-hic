//! Shared durable state for submitted Slurm jobs.
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File},
    io::{self, BufReader, Read, Write},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

pub(super) const SCHEMA: u32 = 1;

#[derive(Serialize, Deserialize)]
pub(super) struct Plan {
    pub(super) schema: u32,
    pub(super) run_id: String,
    pub(super) output: PathBuf,
    pub(super) work: PathBuf,
    pub(super) config_sha256: String,
    pub(super) executable_sha256: String,
    pub(super) submitted_ms: u64,
    pub(super) chunks: usize,
}

#[derive(Serialize, Deserialize)]
pub(super) struct Jobs {
    pub(super) schema: u32,
    pub(super) run_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(super) group_id: Option<String>,
    pub(super) chunks: usize,
    pub(super) index: Option<String>,
    pub(super) alignment: Option<String>,
    pub(super) finish: Option<String>,
}

#[derive(Serialize, Deserialize)]
pub(super) struct IndexDone {
    pub(super) schema: u32,
    pub(super) run_id: String,
    pub(super) chunks: usize,
    pub(super) start_ms: u64,
    pub(super) end_ms: u64,
    pub(super) manifest_sha256: String,
}

#[derive(Serialize, Deserialize)]
pub(super) struct ChunkDone {
    pub(super) schema: u32,
    pub(super) run_id: String,
    pub(super) chunk: usize,
    pub(super) start_ms: u64,
    pub(super) end_ms: u64,
    pub(super) manifest_sha256: String,
    pub(super) input_pairs: u64,
    pub(super) pairs_bytes: u64,
    pub(super) bam_bytes: Option<u64>,
}

pub(super) fn now_ms() -> Result<u64> {
    Ok(u64::try_from(
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_millis(),
    )?)
}

pub(super) fn hash_file(path: &Path) -> Result<String> {
    let mut f =
        BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
    let mut hasher = Sha256::new();
    let mut buf = [0; 128 * 1024];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

pub(super) fn write_json(path: &Path, value: &impl Serialize) -> Result<()> {
    let parent = path.parent().context("missing JSON parent directory")?;
    let mut tmp = tempfile::NamedTempFile::new_in(parent)?;
    serde_json::to_writer_pretty(&mut tmp, value)?;
    writeln!(tmp)?;
    tmp.flush()?;
    tmp.persist(path)
        .with_context(|| format!("publish {}", path.display()))?;
    Ok(())
}

pub(super) fn read_json<T: for<'de> Deserialize<'de>>(path: &Path) -> Result<T> {
    serde_json::from_reader(File::open(path).with_context(|| format!("open {}", path.display()))?)
        .with_context(|| format!("parse {}", path.display()))
}

/// Read a failure record even if a retry removes its path concurrently.
pub(super) fn read_failure(out: &Path) -> Result<Option<serde_json::Value>> {
    let path = out.join("failure.json");
    let file = match File::open(&path) {
        Ok(file) => file,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(None),
        Err(error) => return Err(error).with_context(|| format!("open {}", path.display())),
    };
    let record =
        serde_json::from_reader(file).with_context(|| format!("parse {}", path.display()))?;
    Ok(Some(record))
}

/// Clear only an earlier failure of the phase that is being retried.
pub(super) fn clear_failure(out: &Path, phase: &str) -> Result<()> {
    let path = out.join("failure.json");
    let record = match read_failure(out)? {
        Some(record) => record,
        None => return Ok(()),
    };
    if record.get("phase").and_then(|value| value.as_str()) == Some(phase) {
        match fs::remove_file(&path) {
            Ok(()) => {}
            Err(error) if error.kind() == io::ErrorKind::NotFound => {}
            Err(error) => return Err(error).with_context(|| format!("remove {}", path.display())),
        }
    }
    Ok(())
}

pub(super) fn record_failure(out: &Path, phase: &str, error: &anyhow::Error) {
    let path = out.join("failure.json");
    let Ok(mut file) = tempfile::NamedTempFile::new_in(out) else {
        return;
    };
    let value = serde_json::json!({"status":"failed", "phase":phase, "error":format!("{error:#}")});
    if writeln!(file, "{value}").is_ok() && file.flush().is_ok() {
        let _ = file.persist_noclobber(path);
    }
}
