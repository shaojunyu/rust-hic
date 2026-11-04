//! CLI for building a paired FASTQ index.
use super::IndexArgs;
use crate::mapping::align;
use anyhow::{bail, Context, Result};
use std::{
    env,
    fs::{self, OpenOptions},
    io::{self, Write},
};

fn parse_chunk_size(value: &std::ffi::OsStr) -> Result<u64> {
    let value = value.to_str().context("invalid chunk size")?.trim();
    let split = value
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(value.len());
    let number: u64 = value[..split].parse().context("invalid chunk size")?;
    let multiplier = match value[split..].to_ascii_uppercase().as_str() {
        "" | "B" => 1,
        "K" | "KB" => 1_000,
        "M" | "MB" => 1_000_000,
        "G" | "GB" => 1_000_000_000,
        "KIB" => 1024,
        "MIB" => 1024 * 1024,
        "GIB" => 1024 * 1024 * 1024,
        _ => bail!("invalid chunk size unit; use B, KB, MB, GB, KiB, MiB, or GiB"),
    };
    let size = number
        .checked_mul(multiplier)
        .context("chunk size is too large")?;
    if size == 0 {
        bail!("chunk size must be positive")
    }
    Ok(size)
}
fn format_bytes(bytes: u64) -> String {
    for (unit, divisor) in [
        ("GiB", 1u64 << 30),
        ("MiB", 1u64 << 20),
        ("KiB", 1u64 << 10),
    ] {
        if bytes >= divisor {
            return format!("{:.2} {unit}", bytes as f64 / divisor as f64);
        }
    }
    format!("{bytes} bytes")
}
pub(crate) fn run(args: IndexArgs) -> Result<()> {
    let mut fastq = align::FastqConfig::for_index(args.read1, args.read2);
    let output = args.output.as_path();
    fastq.chunks = args.chunks;
    if let Some(size) = &args.chunk_size {
        fastq.chunk_size_bytes = parse_chunk_size(size.as_ref())?;
    }
    fastq.index_threads_per_mate = args.threads_per_mate;
    fastq.index_search_mib = args.search_mib;
    align::resolve_index(&mut fastq, &env::current_dir()?)?;
    let indexer = env::current_exe()?;
    match fs::create_dir(output) {
        Ok(()) => {}
        Err(error) if error.kind() == io::ErrorKind::AlreadyExists && output.is_dir() => {}
        Err(error) => return Err(error).context("create index directory"),
    }
    let output = fs::canonicalize(output)?;
    for name in [
        "manifest.json",
        "manifest.json.tmp",
        "r1.gzidx",
        "r2.gzidx",
        "r1.gzidx.pending",
        "r2.gzidx.pending",
    ] {
        let path = output.join(name);
        match fs::symlink_metadata(&path) {
            Ok(_) => bail!("index output already exists: {}", path.display()),
            Err(error) if error.kind() == io::ErrorKind::NotFound => {}
            Err(error) => return Err(error).with_context(|| format!("check {}", path.display())),
        }
    }
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output.join("index.log"))
        .context("create index.log (already exists or cannot be created)")?;
    let r1_bytes = fs::metadata(&fastq.read1)?.len();
    let r2_bytes = fs::metadata(&fastq.read2)?.len();
    let (strategy, expected_chunks) = match fastq.chunks {
        Some(n) => (format!("fixed chunk count: {n}"), n as u64),
        None => (
            format!(
                "target compressed chunk size: {} bytes ({})",
                fastq.chunk_size_bytes,
                format_bytes(fastq.chunk_size_bytes)
            ),
            r1_bytes.max(r2_bytes).div_ceil(fastq.chunk_size_bytes),
        ),
    };
    let summary = format!("Indexing paired FASTQ files\n  R1: {} ({r1_bytes} compressed bytes)\n  R2: {} ({r2_bytes} compressed bytes)\n  Splitting: {strategy}\n  Estimated chunks: {expected_chunks}\n  Indexing threads: {} per FASTQ\n  Output: {}\n", fastq.read1.display(), fastq.read2.display(), fastq.index_threads_per_mate, output.display());
    eprint!("{summary}");
    OpenOptions::new()
        .append(true)
        .open(output.join("index.log"))?
        .write_all(summary.as_bytes())?;
    align::index(
        &mut fastq,
        &indexer,
        &output,
        &output.join("index.log"),
        true,
        false,
    )?;
    let complete = format!(
        "Index complete: {} chunks\n",
        fastq.chunks.context("missing index chunk count")?
    );
    eprint!("{complete}");
    OpenOptions::new()
        .append(true)
        .open(output.join("index.log"))?
        .write_all(complete.as_bytes())?;
    println!("{}", output.display());
    Ok(())
}
