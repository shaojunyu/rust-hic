//! Command-line entry point for mapping an indexed paired FASTQ to pairs.
use super::MapArgs;
use crate::mapping::align::{self, FastqConfig};
use crate::pairs::format::read_header;
use anyhow::{bail, Context, Result};
use serde::Deserialize;
use std::{
    env,
    fs::{self, File, OpenOptions},
    io::{self, BufReader, BufWriter, Write},
    os::unix::fs::FileExt,
    path::{Path, PathBuf},
};

#[derive(Deserialize)]
struct ManifestInput {
    path: PathBuf,
}
#[derive(Deserialize)]
struct Manifest {
    inputs: Vec<ManifestInput>,
}

/// Concatenate unsorted chunk bodies in chunk order, with one common header.
fn combine_chunks(chunks: &[PathBuf], staged: &Path) -> Result<()> {
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(staged)?;
    let mut out = BufWriter::with_capacity(1024 * 1024, file);
    let mut expected_header = None;
    for path in chunks {
        let file = File::open(path).with_context(|| format!("open chunk {}", path.display()))?;
        let length = file.metadata()?.len();
        let mut last = [0];
        if length != 0 {
            file.read_exact_at(&mut last, length - 1)?;
        }
        let mut reader = BufReader::with_capacity(1024 * 1024, file);
        let (header, _) = read_header(&mut reader)
            .with_context(|| format!("invalid pairs header in {}", path.display()))?;
        let header: Vec<_> = header
            .into_iter()
            .filter(|line| !line.starts_with("#sorted:"))
            .collect();
        if let Some(expected) = &expected_header {
            if expected != &header {
                bail!("pairs header differs in {}", path.display())
            }
        } else {
            for line in &header {
                out.write_all(line.as_bytes())?;
            }
            expected_header = Some(header);
        }
        let copied = io::copy(&mut reader, &mut out)
            .with_context(|| format!("copy chunk {}", path.display()))?;
        if copied != 0 && last[0] != b'\n' {
            out.write_all(b"\n")?;
        }
    }
    if expected_header.is_none() {
        bail!("no mapping chunks to combine")
    }
    out.flush()?;
    Ok(())
}

pub fn run(args: MapArgs) -> Result<()> {
    let mut aligner_args = Vec::new();
    for value in &args.aligner_arg {
        let parsed = shlex::split(value).context("invalid quoting in --aligner-arg")?;
        if parsed.is_empty() {
            bail!("--aligner-arg must contain at least one argument")
        }
        aligner_args.extend(parsed);
    }
    let index = fs::canonicalize(&args.index_dir).context("open FASTQ index directory")?;
    let manifest: Manifest = serde_json::from_reader(File::open(index.join("manifest.json"))?)
        .context("read FASTQ index manifest")?;
    if manifest.inputs.len() != 2 {
        bail!("FASTQ index must contain two input files")
    }
    let mut config = FastqConfig::for_index(
        manifest.inputs[0].path.clone(),
        manifest.inputs[1].path.clone(),
    );
    config.aligner = args.aligner;
    config.aligner_args = aligner_args;
    config.index_dir = Some(index.clone());
    config.reference = args.reference;
    config.aligner_index = args.aligner_index;
    if let Some(value) = args.aligner_executable {
        config.aligner_executable = value;
    }
    config.jobs = args.jobs;
    let chunk = args.chunk;
    config.threads_per_job = args.threads_per_job;
    config.sam.min_mapq = args.min_mapq;
    config.sam.max_molecule_size = args.max_molecule_size;
    config.sam.max_inter_align_gap = args.max_inter_align_gap;
    align::resolve(&mut config, &env::current_dir()?)?;
    let output = args.output.as_path();
    if output.file_name().is_none() {
        bail!("--output must include a file name")
    }
    match fs::symlink_metadata(output) {
        Ok(_) => bail!("output already exists: {}", output.display()),
        Err(e) if e.kind() == io::ErrorKind::NotFound => {}
        Err(e) => return Err(e).context("check output path"),
    }
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let work = tempfile::Builder::new()
        .prefix(".rust-hic-map-")
        .tempdir_in(parent)
        .context("create temporary mapping directory next to output")?;
    let result = (|| -> Result<()> {
        let exe = env::current_exe()?;
        let idx = align::index(
            &mut config,
            &exe,
            &work.path().join("index"),
            &work.path().join("index.log"),
            false,
            true,
        )?;
        if let Some(i) = chunk {
            eprintln!("Mapping chunk {i} (chunk numbers start at 0)");
        } else {
            eprintln!(
                "Mapping {} chunks with {} concurrent jobs",
                config.chunks.unwrap_or(0),
                config.jobs
            );
        }
        let count = config.chunks.context("missing index chunk count")?;
        if let Some(i) = chunk {
            if i >= count {
                bail!(
                    "chunk {i} is out of range; valid chunk numbers are 0..{}",
                    count - 1
                )
            }
        }
        let chunks = align::map(&config, &idx, work.path(), work.path(), chunk)?;
        let staged = work.path().join("mapped.pairs");
        combine_chunks(&chunks, &staged)?;
        fs::hard_link(&staged, output)
            .with_context(|| format!("publish pairs output {}", output.display()))?;
        eprintln!("Pairs written to {}", output.display());
        Ok(())
    })();
    if let Err(error) = result {
        let path = work.keep();
        return Err(error).with_context(|| {
            format!(
                "mapping failed; logs and intermediate files retained in {}",
                path.display()
            )
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn combines_chunk_bodies_without_claiming_global_sort_order() {
        let dir = tempfile::tempdir().unwrap();
        let header = "## pairs format v1.0.0\n#sorted: chr1-chr2-pos1-pos2\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n";
        let first = dir.path().join("first.pairs");
        let second = dir.path().join("second.pairs");
        let output = dir.path().join("output.pairs");
        // A missing final newline must not join two records together.
        fs::write(&first, format!("{header}first\tchr1\t1\tchr1\t2\t+\t-\tUU")).unwrap();
        fs::write(
            &second,
            format!("{header}second\tchr1\t3\tchr1\t4\t+\t-\tUU\n"),
        )
        .unwrap();
        combine_chunks(&[first, second], &output).unwrap();
        assert_eq!(
            fs::read_to_string(output).unwrap(),
            "## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\nfirst\tchr1\t1\tchr1\t2\t+\t-\tUU\nsecond\tchr1\t3\tchr1\t4\t+\t-\tUU\n"
        );
    }

    #[test]
    fn rejects_incompatible_chunk_headers() {
        let dir = tempfile::tempdir().unwrap();
        let first = dir.path().join("first.pairs");
        let second = dir.path().join("second.pairs");
        let columns = "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n";
        fs::write(&first, format!("#shape: upper triangle\n{columns}")).unwrap();
        fs::write(&second, format!("#shape: whole matrix\n{columns}")).unwrap();
        let error = combine_chunks(&[first, second], &dir.path().join("output.pairs")).unwrap_err();
        assert!(error.to_string().contains("pairs header differs"));
    }
}
