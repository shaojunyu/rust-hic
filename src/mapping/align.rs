mod config;
mod external;
use anyhow::{bail, Context, Result};
pub(crate) use config::DEFAULT_INDEX_THREADS_PER_MATE;
pub use config::{resolve, resolve_index, FastqConfig};
use external::align_one;
use serde::{Deserialize, Serialize};
use std::{
    ffi::CString,
    fs::{self, File, OpenOptions},
    io::{self, BufRead, BufReader, BufWriter, Write},
    os::unix::{
        ffi::OsStrExt,
        fs::{OpenOptionsExt, PermissionsExt},
    },
    path::{Path, PathBuf},
    process::{Child, Command, Stdio},
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering},
        mpsc, Arc, Mutex,
    },
    thread,
    time::{Duration, Instant},
};
fn run_index_command(
    cmd: &mut Command,
    log: &Path,
    show: bool,
    append: bool,
) -> Result<std::process::ExitStatus> {
    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(append)
        .truncate(!append)
        .open(log)?;
    cmd.stdout(Stdio::null());
    if !show {
        return cmd
            .stderr(Stdio::from(file))
            .status()
            .context("start FASTQ index process");
    }
    let mut child = cmd
        .stderr(Stdio::piped())
        .spawn()
        .context("start FASTQ index process")?;
    let stderr = child.stderr.take().context("capture FASTQ index log")?;
    let mut reader = BufReader::new(stderr);
    let mut file = file;
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        file.write_all(&line)?;
        io::stderr().write_all(&line)?;
    }
    child.wait().context("wait for FASTQ index process")
}
/// Check reusable-index metadata before creating output or starting alignment.
pub(crate) fn validate_reusable_index(c: &FastqConfig, path: &Path) -> Result<usize> {
    #[derive(Deserialize)]
    struct IndexInput {
        path: PathBuf,
        index: PathBuf,
    }
    #[derive(Deserialize)]
    struct IndexManifest {
        version: u32,
        inputs: Vec<IndexInput>,
        pairs: Option<u64>,
        chunks: usize,
    }
    let manifest_path = path.join("manifest.json");
    let manifest: IndexManifest = serde_json::from_reader(
        File::open(&manifest_path)
            .with_context(|| format!("open FASTQ index manifest {}", manifest_path.display()))?,
    )
    .with_context(|| format!("read FASTQ index manifest {}", manifest_path.display()))?;
    if manifest.version != 1
        || manifest.inputs.len() != 2
        || !(1..=256).contains(&manifest.chunks)
        || c.chunks.is_some_and(|n| manifest.chunks != n)
    {
        bail!("FASTQ index version, mate count or chunk count does not match configuration")
    }
    for (input, expected) in manifest.inputs.iter().zip([&c.read1, &c.read2]) {
        let actual = fs::canonicalize(&input.path)
            .with_context(|| format!("open indexed FASTQ {}", input.path.display()))?;
        if &actual != expected {
            bail!("FASTQ index was built for different input mates or mate order")
        }
        let index_path = path.join(&input.index);
        let file = File::open(&index_path)
            .with_context(|| format!("open gzip index {}", index_path.display()))?;
        if !file.metadata()?.is_file() {
            bail!(
                "gzip index must be a regular file: {}",
                index_path.display()
            );
        }
    }
    if let (Some(actual), Some(expected)) = (manifest.pairs, c.expected_pairs) {
        if actual != expected {
            bail!("FASTQ index pair count {actual} does not match expected {expected}")
        }
    }
    Ok(manifest.chunks)
}

pub fn index(
    c: &mut FastqConfig,
    exe: &Path,
    new_index: &Path,
    log: &Path,
    show: bool,
    verify_reused: bool,
) -> Result<PathBuf> {
    let reused = c.index_dir.is_some();
    let path = if let Some(path) = &c.index_dir {
        validate_reusable_index(c, path)?;
        fs::write(log, format!("event=index_reuse path={}\n", path.display()))?;
        path.clone()
    } else {
        let path = new_index.to_path_buf();
        let mut cmd = Command::new(exe);
        cmd.arg("index").arg("--internal");
        cmd.arg("build")
            .arg(&c.read1)
            .arg(&c.read2)
            .arg(&path)
            .arg(
                c.chunks
                    .map(|n| n.to_string())
                    .unwrap_or_else(|| format!("size:{}", c.chunk_size_bytes)),
            )
            .arg(c.index_threads_per_mate.to_string())
            .arg(c.index_search_mib.to_string());
        let status = run_index_command(&mut cmd, log, show, show)?;
        if !status.success() {
            bail!("FASTQ index failed: {status}; inspect {}", log.display())
        }
        path
    };
    // A new index verifies itself during build. Workflow preflight verifies a
    // reusable index; standalone map requests verification here instead.
    if reused && verify_reused {
        let mut verify = Command::new(exe);
        verify
            .arg("index")
            .arg("--internal")
            .arg("verify")
            .arg(&path);
        let status = run_index_command(&mut verify, log, show, true)?;
        if !status.success() {
            bail!("FASTQ index verification failed: {status}")
        };
    }
    let manifest: serde_json::Value =
        serde_json::from_reader(File::open(path.join("manifest.json"))?)?;
    let chunks = manifest
        .get("chunks")
        .and_then(|n| n.as_u64())
        .context("missing index chunk count")?;
    let chunks = usize::try_from(chunks)?;
    if !(1..=256).contains(&chunks) {
        bail!("index chunk count must be between 1 and 256")
    }
    c.chunks = Some(chunks);
    if let Some(expected) = c.expected_pairs {
        let manifest: serde_json::Value =
            serde_json::from_reader(File::open(path.join("manifest.json"))?)?;
        if let Some(actual) = manifest.get("pairs").and_then(|v| v.as_u64()) {
            if actual != expected {
                bail!("FASTQ index pair count {actual} does not match expected {expected}")
            }
        }
    }
    Ok(path)
}
pub fn map(
    c: &FastqConfig,
    index: &Path,
    work: &Path,
    logs: &Path,
    chunk: Option<usize>,
) -> Result<Vec<PathBuf>> {
    let cancelled = AtomicBool::new(false);
    map_with_completed(
        c,
        index,
        work,
        logs,
        chunk,
        &cancelled,
        false,
        |_, _| Ok(()),
    )
}

/// Invoke `completed` as soon as each chunk has been completely written.
#[allow(clippy::too_many_arguments)] // Mapping context and completion callback.
pub(crate) fn map_with_completed<F>(
    c: &FastqConfig,
    index: &Path,
    work: &Path,
    logs: &Path,
    chunk: Option<usize>,
    cancelled: &AtomicBool,
    keep_bam: bool,
    completed: F,
) -> Result<Vec<PathBuf>>
where
    F: Fn(usize, &Path) -> Result<()> + Sync,
{
    let chunks = c
        .chunks
        .context("FASTQ index chunk count was not resolved")?;
    if let Some(i) = chunk {
        if i >= chunks {
            bail!(
                "chunk {i} is out of range; valid chunk numbers are 0..{}",
                chunks - 1
            )
        }
    }
    let selected_chunks = if chunk.is_some() { 1 } else { chunks };
    let mem3 = if c.aligner == "bwa-mem3" {
        eprintln!(
            "Loading built-in bwa-mem3 reference index: {}",
            c.aligner_index.display()
        );
        Some((
            Arc::new(bwa_mem3_rs::BwaIndex::load(&c.aligner_index)?),
            Arc::new(super::bwa_mem3::options(c)?),
        ))
    } else {
        None
    };
    let next = AtomicUsize::new(0);
    let output = Mutex::new(Vec::new());
    let errors = Mutex::new(Vec::new());
    thread::scope(|scope| {
        for _ in 0..c.jobs.min(selected_chunks) {
            let next = &next;
            let output = &output;
            let errors = &errors;
            let mem3 = &mem3;
            let completed = &completed;
            scope.spawn(move || loop {
                if cancelled.load(Ordering::Acquire) {
                    break;
                }
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= selected_chunks {
                    break;
                }
                let i = chunk.unwrap_or(i);
                match align_one(
                    c,
                    index,
                    work,
                    logs,
                    i,
                    mem3.as_ref().map(|(index, opts)| (index, opts)),
                    cancelled,
                    keep_bam,
                ) {
                    Ok((p, n)) => match completed(i, &p) {
                        Ok(()) => output.lock().unwrap().push((i, p, n)),
                        Err(e) => {
                            cancelled.store(true, Ordering::Release);
                            errors.lock().unwrap().push(format!("chunk={i}: {e:#}"));
                        }
                    },
                    Err(e) => {
                        cancelled.store(true, Ordering::Release);
                        errors.lock().unwrap().push(format!("chunk={i}: {e:#}"));
                    }
                }
            });
        }
    });
    let errors = errors.into_inner().unwrap();
    if !errors.is_empty() {
        bail!("alignment failed: {}", errors.join("; "))
    };
    if cancelled.load(Ordering::Acquire) {
        bail!("alignment cancelled")
    }
    let mut out = output.into_inner().unwrap();
    out.sort_by_key(|x| x.0);
    let total: u64 = out.iter().map(|x| x.2).sum();
    let manifest: serde_json::Value =
        serde_json::from_reader(File::open(index.join("manifest.json"))?)?;
    for expected in [
        manifest.get("pairs").and_then(|x| x.as_u64()),
        c.expected_pairs,
    ]
    .into_iter()
    .flatten()
    .filter(|_| chunk.is_none())
    {
        if total != expected {
            bail!("aligner read {total} pairs but expected {expected}")
        }
    }
    eprintln!(
        "event=alignment_done chunks={} input_pairs={total}",
        out.len()
    );
    Ok(out.into_iter().map(|(_, p, _)| p).collect())
}
