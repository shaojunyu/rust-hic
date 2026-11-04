//! Streaming in-silico Arima digestion compatible with runHiC's fragment BED.
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct DigestConfig {
    pub enzyme: String,
    pub reference: PathBuf,
}
pub fn resolve(c: &mut DigestConfig, base: &Path) -> Result<()> {
    if !c.enzyme.eq_ignore_ascii_case("arima") {
        bail!(
            "unsupported digest enzyme '{}'; currently supported: Arima",
            c.enzyme
        )
    }
    let path = if c.reference.is_absolute() {
        c.reference.clone()
    } else {
        base.join(&c.reference)
    };
    c.reference = fs::canonicalize(&path)
        .with_context(|| format!("open digest reference FASTA {}", path.display()))?;
    Ok(())
}
pub fn digest(c: &DigestConfig, sizes: &Path, output: &Path) -> Result<()> {
    let mut chroms: Vec<(String, u64, Vec<u64>, bool)> = Vec::new();
    let mut lookup = HashMap::new();
    for line in BufReader::new(File::open(sizes)?).lines() {
        let line = line?;
        let f: Vec<_> = line.split_whitespace().collect();
        if f.len() != 2 {
            bail!("invalid chromosome sizes row")
        }
        let len: u64 = f[1].parse()?;
        if len == 0 || lookup.insert(f[0].to_owned(), chroms.len()).is_some() {
            bail!("invalid or duplicate chromosome")
        }
        chroms.push((f[0].into(), len, Vec::new(), false));
    }
    if chroms.is_empty() {
        bail!("empty chromosome sizes")
    }
    let mut current = None;
    let mut len = 0u64;
    let mut window = [0u8; 5];
    let finish = |current: Option<usize>,
                  len: u64,
                  chroms: &mut Vec<(String, u64, Vec<u64>, bool)>|
     -> Result<()> {
        if let Some(i) = current {
            if len != chroms[i].1 {
                bail!(
                    "reference length for {} is {len}, expected {}",
                    chroms[i].0,
                    chroms[i].1
                )
            }
        }
        Ok(())
    };
    for line in BufReader::new(File::open(&c.reference)?).lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            finish(current, len, &mut chroms)?;
            let name = header
                .split_whitespace()
                .next()
                .context("empty FASTA name")?;
            current = lookup.get(name).copied();
            len = 0;
            window = [0; 5];
            if let Some(i) = current {
                if chroms[i].3 {
                    bail!("duplicate FASTA chromosome {name}")
                }
                chroms[i].3 = true;
            }
            continue;
        }
        if let Some(i) = current {
            for mut base in line.bytes() {
                if !base.is_ascii_alphabetic() {
                    bail!("invalid FASTA character on {}", chroms[i].0)
                }
                base.make_ascii_uppercase();
                window.rotate_left(1);
                window[4] = base;
                // Bio.Restriction excludes MboI sites with cuts at either end.
                // runHiC adds one to the returned 1-based cut coordinates.
                if len >= 4 && len - 4 > 0 && &window[0..4] == b"GATC" {
                    chroms[i].2.push(len - 2);
                }
                // HinfI GA^NTC has cuts inside the recognition sequence.
                if len >= 4
                    && window[0] == b'G'
                    && window[1] == b'A'
                    && window[3] == b'T'
                    && window[4] == b'C'
                {
                    chroms[i].2.push(len - 1);
                }
                len = len.checked_add(1).context("FASTA length overflow")?;
            }
        }
    }
    finish(current, len, &mut chroms)?;
    let mut out = BufWriter::with_capacity(
        1024 * 1024,
        File::options().write(true).create_new(true).open(output)?,
    );
    for (name, len, cuts, found) in &mut chroms {
        if !*found {
            bail!("chromosome {name} missing from FASTA")
        }
        cuts.sort_unstable();
        let mut prev = 0;
        for &cut in cuts.iter() {
            if cut > *len {
                bail!("cut outside chromosome {name}")
            }
            writeln!(out, "{name}\t{prev}\t{cut}")?;
            prev = cut;
        }
        writeln!(out, "{name}\t{prev}\t{len}")?;
    }
    out.flush()?;
    Ok(())
}
