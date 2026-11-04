//! Internal FASTQ index command dispatch.
use super::*;

fn probe(out: &Path, boundary: usize, mate: usize, count: u64) -> Result<()> {
    let m: Manifest = serde_json::from_reader(File::open(out.join("manifest.json"))?)?;
    if boundary >= m.chunks || !(1..=2).contains(&mate) {
        bail!("invalid boundary or mate")
    }
    let i = &m.inputs[mate - 1];
    let f = fingerprint(Path::new(&i.path), i.index.clone(), i.uncompressed_size)?;
    if (f.size, f.modified_ns, &f.edge_sha256) != (i.size, i.modified_ns, &i.edge_sha256) {
        bail!("input changed")
    }
    let b = &m.boundaries[boundary];
    let offset = if mate == 1 {
        b.r1_uncompressed_offset
    } else {
        b.r2_uncompressed_offset
    };
    let mut r = open_indexed_reader(i, out, offset, 1)?;
    let n = io::copy(&mut r.by_ref().take(count), &mut io::stdout().lock())?;
    if n != count {
        bail!("short read from gzip")
    };
    Ok(())
}
pub(crate) fn main() -> Result<()> {
    let args: Vec<String> = std::iter::once("rust-hic index".to_owned())
        .chain(env::args().skip(3))
        .collect();
    match args.get(1).map(String::as_str) {
        Some("stream") if args.len() == 5 => {
            stream_chunk(Path::new(&args[2]), args[3].parse()?, args[4].parse()?)
        }
        Some("probe") if args.len() == 6 => probe(
            Path::new(&args[2]),
            args[3].parse()?,
            args[4].parse()?,
            args[5].parse()?,
        ),
        Some("build") if (6..=8).contains(&args.len()) => build_index(
            Path::new(&args[2]),
            Path::new(&args[3]),
            Path::new(&args[4]),
            &args[5],
            args.get(6).map(|x| x.parse()).transpose()?.unwrap_or(1),
            args.get(7)
                .map(|x| x.parse::<u64>())
                .transpose()?
                .unwrap_or(128)
                .checked_mul(1024 * 1024)
                .context("window overflow")?,
        ),
        Some("extract") if (6..=7).contains(&args.len()) => extract(
            Path::new(&args[2]),
            args[3].parse()?,
            [Path::new(&args[4]), Path::new(&args[5])],
            args.get(6).map(|x| x.parse()).transpose()?.unwrap_or(2),
        ),
        Some("compact") if (3..=4).contains(&args.len()) => compact_existing(
            Path::new(&args[2]),
            args.get(3).map(|x| x.parse()).transpose()?.unwrap_or(2),
        ),
        Some("verify") if (3..=4).contains(&args.len()) => verify(
            Path::new(&args[2]),
            args.get(3).map(|x| x.parse()).transpose()?.unwrap_or(2),
        ),
        _ => {
            eprintln!("invalid internal FASTQ index command");
            std::process::exit(2)
        }
    }
}
