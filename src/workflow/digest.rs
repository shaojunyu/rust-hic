//! Streaming in-silico restriction digestion compatible with runHiC's BED.
//!
//! `Arima` is treated as the MboI + HinfI cocktail. Other names and
//! `-`-separated cocktails use the generated `restriction` table.
//!
//! For a recognition site starting at 0-based position `s`, Biopython's
//! top-strand cut is `s + fst5 + 1` (1-based) and runHiC writes the fragment
//! boundary as `cut + 1`. Biopython's `_drop` filter also requires both the
//! Watson and Crick cuts to be strictly inside the chromosome, which is
//! reproduced with `ovhg`.

use super::restriction;
use anyhow::{anyhow, bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct DigestConfig {
    /// `Arima`, one enzyme name, or a `-`-separated enzyme cocktail.
    pub enzyme: String,
    pub reference: PathBuf,
}

#[derive(Clone, Copy, Debug)]
pub struct DigestStats {
    pub chromosomes: usize,
    pub cuts: u64,
    pub fragments: u64,
    pub elapsed_seconds: f64,
}

pub fn resolve(c: &mut DigestConfig, base: &Path) -> Result<()> {
    restriction::parse_spec(&c.enzyme)
        .map_err(|e| anyhow!("unsupported digest enzyme list '{}': {e}", c.enzyme))?;
    let path = if c.reference.is_absolute() {
        c.reference.clone()
    } else {
        base.join(&c.reference)
    };
    c.reference = fs::canonicalize(&path)
        .with_context(|| format!("open digest reference FASTA {}", path.display()))?;
    Ok(())
}

pub fn digest(c: &DigestConfig, sizes: &Path, output: &Path) -> Result<DigestStats> {
    let enzymes = restriction::parse_spec(&c.enzyme)
        .map_err(|e| anyhow!("unsupported digest enzyme list '{}': {e}", c.enzyme))?;
    let max_site = enzymes
        .iter()
        .map(|enzyme| enzyme.site.len())
        .max()
        .context("empty restriction enzyme list")?;
    if max_site == 0 {
        bail!("restriction enzymes must have non-empty recognition sites")
    }
    let started = Instant::now();

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
    let mut window = vec![0u8; max_site];
    let mut pos = 0u64;
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
            finish(current, pos, &mut chroms)?;
            let name = header
                .split_whitespace()
                .next()
                .context("empty FASTA name")?;
            current = lookup.get(name).copied();
            pos = 0;
            window.fill(0);
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
                window[max_site - 1] = base;
                for enzyme in &enzymes {
                    let site_len = enzyme.site.len();
                    if pos + 1 >= site_len as u64 {
                        let start = pos as i64 + 1 - site_len as i64;
                        let site_window = &window[max_site - site_len..];
                        if matches(enzyme.site.as_bytes(), site_window) {
                            let cut = start + i64::from(enzyme.fst5) + 1;
                            let crick = cut - i64::from(enzyme.ovhg);
                            let chrom_len = chroms[i].1 as i64;
                            if cut > 1 && cut <= chrom_len && crick > 1 && crick <= chrom_len {
                                let boundary = cut + 1;
                                if boundary <= chrom_len {
                                    chroms[i].2.push(boundary as u64);
                                }
                            }
                        }
                    }
                }
                pos = pos.checked_add(1).context("FASTA length overflow")?;
            }
        }
    }
    finish(current, pos, &mut chroms)?;

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
    let cuts: u64 = chroms.iter().map(|chrom| chrom.2.len() as u64).sum();
    let chromosomes = chroms.len();
    Ok(DigestStats {
        chromosomes,
        cuts,
        fragments: cuts + chromosomes as u64,
        elapsed_seconds: started.elapsed().as_secs_f64(),
    })
}

fn base_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        _ => 0,
    }
}

fn pattern_mask(code: u8) -> u8 {
    match code.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}

fn matches(pattern: &[u8], sequence: &[u8]) -> bool {
    pattern.len() == sequence.len()
        && pattern.iter().zip(sequence).all(|(&code, &base)| {
            let mask = pattern_mask(code);
            if mask == 0 {
                return false;
            }
            // Biopython turns an ambiguous pattern base into a wildcard that
            // matches any reference character, including N/R/Y. A specific
            // pattern base, however, does not match an ambiguous reference base.
            if mask == 15 {
                return base.is_ascii_alphabetic();
            }
            let base = base_mask(base);
            base != 0 && (mask & base) != 0
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SEQ1: &str = "AAAAGATCAAAGAGTCAAAAAGCTTAAAGAATTCAAAAGATCTAAAACATGAAAAACTGAAAAAGCTAAAAGGCCAAAATTAACAAAATCGAA";
    const SEQ2: &str = "TTTTGATCGGGGGAGTCAAAAAAGCTTTTTTGAATTCTTTTTAGATCTTTTTT";
    const SIZES: &str = "chr1\t93\nchr2\t53\n";

    fn fasta() -> String {
        format!(">chr1\n{SEQ1}\n>chr2\n{SEQ2}\n")
    }

    fn run(enzyme: &str, fasta: &str, sizes: &str) -> Vec<String> {
        let dir = tempfile::tempdir().unwrap();
        let reference = dir.path().join("reference.fa");
        let sizes_path = dir.path().join("chrom.sizes");
        let output = dir.path().join("fragments.bed");
        fs::write(&reference, fasta).unwrap();
        fs::write(&sizes_path, sizes).unwrap();
        let config = DigestConfig {
            enzyme: enzyme.to_owned(),
            reference,
        };
        digest(&config, &sizes_path, &output).unwrap();
        fs::read_to_string(output)
            .unwrap()
            .lines()
            .map(str::to_owned)
            .collect()
    }

    fn expected(text: &str) -> Vec<String> {
        text.lines().map(str::to_owned).collect()
    }

    #[test]
    fn arima_matches_runhic() {
        assert_eq!(
            run("Arima", &fasta(), SIZES),
            expected(
                "chr1\t0\t6\nchr1\t6\t14\nchr1\t14\t40\nchr1\t40\t93\n\
                 chr2\t0\t6\nchr2\t6\t15\nchr2\t15\t45\nchr2\t45\t53"
            )
        );
    }

    #[test]
    fn digest_stats_match_the_bed() {
        let dir = tempfile::tempdir().unwrap();
        let reference = dir.path().join("reference.fa");
        let sizes = dir.path().join("chrom.sizes");
        let output = dir.path().join("fragments.bed");
        fs::write(&reference, fasta()).unwrap();
        fs::write(&sizes, SIZES).unwrap();
        let config = DigestConfig {
            enzyme: "Arima".to_owned(),
            reference,
        };
        let stats = digest(&config, &sizes, &output).unwrap();
        assert_eq!(stats.chromosomes, 2);
        assert_eq!(stats.cuts, 6);
        assert_eq!(stats.fragments, 8);
        assert!(stats.elapsed_seconds >= 0.0);
    }

    #[test]
    fn common_enzymes_match_runhic() {
        for (enzyme, bed) in [
            (
                "MboI",
                "chr1\t0\t6\nchr1\t6\t40\nchr1\t40\t93\nchr2\t0\t6\nchr2\t6\t45\nchr2\t45\t53",
            ),
            (
                "HinfI",
                "chr1\t0\t14\nchr1\t14\t93\nchr2\t0\t15\nchr2\t15\t53",
            ),
            (
                "HindIII",
                "chr1\t0\t22\nchr1\t22\t93\nchr2\t0\t24\nchr2\t24\t53",
            ),
            (
                "EcoRI",
                "chr1\t0\t31\nchr1\t31\t93\nchr2\t0\t34\nchr2\t34\t53",
            ),
            (
                "BglII",
                "chr1\t0\t40\nchr1\t40\t93\nchr2\t0\t45\nchr2\t45\t53",
            ),
            ("NlaIII", "chr1\t0\t53\nchr1\t53\t93\nchr2\t0\t53"),
            (
                "AluI",
                "chr1\t0\t24\nchr1\t24\t67\nchr1\t67\t93\nchr2\t0\t26\nchr2\t26\t53",
            ),
            ("HaeIII", "chr1\t0\t75\nchr1\t75\t93\nchr2\t0\t53"),
            (
                "MseI",
                "chr1\t0\t26\nchr1\t26\t82\nchr1\t82\t93\nchr2\t0\t53",
            ),
            ("TaqI", "chr1\t0\t91\nchr1\t91\t93\nchr2\t0\t53"),
        ] {
            assert_eq!(run(enzyme, &fasta(), SIZES), expected(bed), "{enzyme}");
        }
    }

    #[test]
    fn cocktail_matches_arima() {
        assert_eq!(
            run("MboI-HinfI", &fasta(), SIZES),
            run("Arima", &fasta(), SIZES)
        );
    }

    #[test]
    fn non_palindromic_bsai_matches_runhic() {
        let seq = format!("{}GGTCTC{}", "A".repeat(20), "A".repeat(20));
        assert_eq!(
            run("BsaI", &format!(">x\n{seq}\n"), "x\t46\n"),
            expected("x\t0\t29\nx\t29\t46")
        );
    }

    #[test]
    fn sites_at_chromosome_edges_are_excluded() {
        let fasta = ">e1\nGATCAAAA\n>e2\nAAAAGATC\n>e3\nGATC\n";
        assert_eq!(
            run("MboI", fasta, "e1\t8\ne2\t8\ne3\t4\n"),
            expected("e1\t0\t8\ne2\t0\t8\ne3\t0\t4")
        );
    }

    #[test]
    fn lowercase_fasta_is_supported() {
        let lower = format!(
            ">chr1\n{}\n>chr2\n{}\n",
            SEQ1.to_ascii_lowercase(),
            SEQ2.to_ascii_lowercase()
        );
        assert_eq!(run("Arima", &lower, SIZES), run("Arima", &fasta(), SIZES));
    }

    #[test]
    fn unknown_enzyme_is_rejected() {
        assert!(restriction::parse_spec("NotAnEnzyme").is_err());
        assert!(restriction::parse_spec("").is_err());
    }

    #[test]
    fn enzyme_lookup_is_case_insensitive() {
        assert!(restriction::parse_spec("mboi").is_ok());
        assert!(restriction::parse_spec("MBOI").is_ok());
    }

    #[test]
    fn enzyme_names_are_case_insensitive_in_digestion() {
        let fasta = fasta();
        let baseline = run("MboI", &fasta, SIZES);
        for name in ["mboi", "MBOI", "Mboi"] {
            assert_eq!(run(name, &fasta, SIZES), baseline, "{name}");
        }
        assert_eq!(run("arima", &fasta, SIZES), run("Arima", &fasta, SIZES));
    }

    const SEQ_A: &str = "AAAAGATCAAACCGGAAACCCGGGAAAA";
    const SEQ_B: &str = "TTTTGATCGGGGCCGGTTTCCCGGGTTT";
    const ISO_SIZES: &str = "chrA\t28\nchrB\t28\n";

    fn iso_fasta() -> String {
        format!(">chrA\n{SEQ_A}\n>chrB\n{SEQ_B}\n")
    }

    #[test]
    fn isoschizomers_produce_identical_beds() {
        let fasta = iso_fasta();
        for group in [&["MboI", "DpnII", "Sau3AI"][..], &["HpaII", "MspI"][..]] {
            let first = run(group[0], &fasta, ISO_SIZES);
            for enzyme in &group[1..] {
                assert_eq!(run(enzyme, &fasta, ISO_SIZES), first, "{enzyme}");
            }
        }
    }

    #[test]
    fn neoschizomers_match_runhic_and_differ() {
        let fasta = iso_fasta();
        assert_eq!(
            run("SmaI", &fasta, ISO_SIZES),
            expected("chrA\t0\t23\nchrA\t23\t28\nchrB\t0\t24\nchrB\t24\t28")
        );
        assert_eq!(
            run("XmaI", &fasta, ISO_SIZES),
            expected("chrA\t0\t21\nchrA\t21\t28\nchrB\t0\t22\nchrB\t22\t28")
        );
        assert_ne!(
            run("SmaI", &fasta, ISO_SIZES),
            run("XmaI", &fasta, ISO_SIZES)
        );
    }

    #[test]
    fn ambiguous_reference_bases_match_biopython() {
        let fasta = ">n1\nAAAAGANTCAAAA\n>n2\nAAAAGANCAAAA\n";
        let sizes = "n1\t13\nn2\t12\n";
        assert_eq!(
            run("HinfI", fasta, sizes),
            expected("n1\t0\t7\nn1\t7\t13\nn2\t0\t12")
        );
        assert_eq!(run("MboI", fasta, sizes), expected("n1\t0\t13\nn2\t0\t12"));
    }

    #[test]
    #[ignore = "requires RUST_HIC_DIGEST_FASTA/SIZES/ENZYME/OUT"]
    fn write_bed_for_external_comparison() {
        let Ok(fasta) = std::env::var("RUST_HIC_DIGEST_FASTA") else {
            return;
        };
        let sizes = std::env::var("RUST_HIC_DIGEST_SIZES").unwrap();
        let enzyme = std::env::var("RUST_HIC_DIGEST_ENZYME").unwrap();
        let output = std::env::var("RUST_HIC_DIGEST_OUT").unwrap();
        let config = DigestConfig {
            enzyme,
            reference: fasta.into(),
        };
        digest(
            &config,
            std::path::Path::new(&sizes),
            std::path::Path::new(&output),
        )
        .unwrap();
    }
}
