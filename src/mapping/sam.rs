//! BWA SAM to selected Hi-C pairs (pairtools parse's 5unique policy).
use anyhow::{bail, Context, Result};
use std::io::{BufRead, Write};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Alignment {
    chrom: String,
    pos: i64,
    strand: char,
    mapped: bool,
    unique: bool,
    kind: char,
    dist5: i64,
    read_span: i64,
    query_start: i64,
}

/// Convert bwa-mem3's structured fields using the same 5unique semantics as SAM.
pub(crate) fn from_fields(
    fields: &bwa_mem3_rs::AlignedFields<'_>,
    index: &bwa_mem3_rs::BwaIndex,
    min_mapq: u8,
) -> Result<Alignment> {
    let chrom = if fields.flag & 4 == 0 {
        if fields.tid < 0 || fields.tid as usize >= index.n_contigs() {
            bail!("mapped BWA record lacks reference coordinates")
        }
        index.contig_name(fields.tid as usize)
    } else {
        "*"
    };
    from_fields_with_chrom(fields, chrom, min_mapq)
}

fn from_fields_with_chrom(
    fields: &bwa_mem3_rs::AlignedFields<'_>,
    chrom: &str,
    min_mapq: u8,
) -> Result<Alignment> {
    let flag = fields.flag;
    if flag & 1 == 0 {
        bail!("unpaired BWA alignment")
    }
    let mapped = flag & 4 == 0;
    let unique = fields.mapq >= min_mapq;
    let mut a = Alignment::empty();
    a.mapped = mapped;
    a.unique = unique;
    a.kind = if !mapped {
        'N'
    } else if !unique {
        'M'
    } else {
        'U'
    };
    if mapped {
        if fields.tid < 0 || fields.pos < 0 || chrom == "*" {
            bail!("mapped BWA record lacks reference coordinates")
        }
        let (mut matched, mut ref_span, mut read_span, mut leading, mut trailing, mut soft_start) =
            (0i64, 0i64, 0i64, 0i64, 0i64, 0i64);
        let mut before_alignment = true;
        for &c in fields.cigar {
            let n = i64::from(c >> 4);
            if n == 0 {
                bail!("invalid CIGAR")
            }
            match c & 15 {
                0 => {
                    matched += n;
                    ref_span += n;
                    read_span += n;
                    before_alignment = false;
                }
                1 => {
                    read_span += n;
                    before_alignment = false;
                }
                2 => {
                    ref_span += n;
                    before_alignment = false;
                }
                op @ (4 | 5) => {
                    if matched == 0 {
                        leading = n
                    } else {
                        trailing = n
                    }
                    if op == 4 && before_alignment {
                        soft_start += n
                    }
                }
                _ => bail!("unsupported BWA CIGAR operation"),
            }
        }
        if ref_span == 0 {
            bail!("invalid or empty mapped CIGAR")
        }
        let reverse = flag & 16 != 0;
        a.dist5 = if reverse { trailing } else { leading };
        a.read_span = read_span;
        a.query_start = soft_start;
        if unique {
            a.chrom = chrom.to_owned();
            let pos = i64::from(fields.pos) + 1;
            a.pos = if reverse { pos + ref_span - 1 } else { pos };
            a.strand = if reverse { '-' } else { '+' };
        }
    }
    Ok(a)
}
impl Alignment {
    fn empty() -> Self {
        Self {
            chrom: "!".into(),
            pos: 0,
            strand: '-',
            mapped: false,
            unique: false,
            kind: 'N',
            dist5: 0,
            read_span: 0,
            query_start: 0,
        }
    }
}
fn parse_alignment(f: &[&str], min_mapq: u8) -> Result<(bool, Alignment)> {
    if f.len() < 11 {
        bail!("SAM record has fewer than 11 fields")
    }
    let flag: u16 = f[1].parse().context("invalid SAM flag")?;
    if flag & 0x1 == 0 {
        bail!("unpaired SAM record")
    }
    let mapq: u8 = f[4].parse().context("invalid SAM MAPQ")?;
    let mapped = flag & 4 == 0;
    let unique = mapq >= min_mapq;
    let mut a = Alignment::empty();
    a.mapped = mapped;
    a.unique = unique;
    a.kind = if !mapped {
        'N'
    } else if !unique {
        'M'
    } else {
        'U'
    };
    if mapped {
        let pos: i64 = f[3].parse().context("invalid SAM position")?;
        if pos < 1 || f[2] == "*" || f[5] == "*" {
            bail!("mapped SAM record lacks coordinates or CIGAR")
        }
        let (mut n, mut matched, mut ref_span, mut read_span, mut leading, mut trailing) =
            (0i64, 0, 0, 0, 0, 0);
        let mut soft_start = 0;
        let mut before_alignment = true;
        for ch in f[5].chars() {
            if let Some(d) = ch.to_digit(10) {
                n = n
                    .checked_mul(10)
                    .and_then(|x| x.checked_add(d as i64))
                    .context("CIGAR overflow")?;
                continue;
            }
            if n == 0 {
                bail!("invalid CIGAR")
            }
            match ch {
                'M' => {
                    matched += n;
                    ref_span += n;
                    read_span += n;
                    before_alignment = false;
                }
                'I' => {
                    read_span += n;
                    before_alignment = false;
                }
                'D' => {
                    ref_span += n;
                    before_alignment = false;
                }
                'S' | 'H' => {
                    if matched == 0 {
                        leading = n
                    } else {
                        trailing = n
                    }
                    if ch == 'S' && before_alignment {
                        soft_start += n
                    }
                }
                _ => bail!("unsupported CIGAR operation {ch}; expected BWA M/I/D/S/H"),
            }
            n = 0;
        }
        if n != 0 || ref_span == 0 {
            bail!("invalid or empty mapped CIGAR")
        }
        let reverse = flag & 16 != 0;
        a.dist5 = if reverse { trailing } else { leading };
        a.read_span = read_span;
        a.query_start = soft_start;
        if unique {
            a.chrom = f[2].to_owned();
            a.pos = if reverse {
                pos.checked_add(ref_span - 1)
                    .context("coordinate overflow")?
            } else {
                pos
            };
            a.strand = if reverse { '-' } else { '+' };
        }
    }
    Ok((flag & 64 != 0, a))
}
// Order alignments by their 5' distance and represent internal uncovered
// segments as gaps. Unaligned sequence before the first alignment is clipping,
// not an additional alignment in the walk.
fn normalize(mut a: Vec<Alignment>, max_gap: i64) -> Vec<Alignment> {
    if a.is_empty() {
        return vec![Alignment::empty()];
    }
    a.sort_by_key(|x| x.query_start);
    a.sort_by_key(|x| x.dist5);
    if a.len() == 1 && !a[0].mapped {
        return a;
    }
    let mut alignments = a.into_iter();
    let first = alignments.next().expect("nonempty alignment list");
    let mut end = first.dist5 + first.read_span;
    let mut normalized = vec![first];
    for alignment in alignments {
        if alignment.dist5 - end > max_gap {
            let mut gap = Alignment::empty();
            gap.dist5 = end;
            gap.read_span = alignment.dist5 - end;
            normalized.push(gap);
        }
        end = end.max(alignment.dist5 + alignment.read_span);
        normalized.push(alignment);
    }
    normalized
}
fn classify(
    mut a: Vec<Alignment>,
    mut b: Vec<Alignment>,
    max_molecule: i64,
    max_gap: i64,
) -> Option<(Alignment, Alignment)> {
    if a.is_empty() || b.is_empty() {
        return None;
    }
    a = normalize(a, max_gap);
    b = normalize(b, max_gap);
    let mut rescued = false;
    if (a.len() == 2 && b.len() == 1) || (a.len() == 1 && b.len() == 2) {
        let (chim, linear) = if a.len() == 2 {
            (&a, &b[0])
        } else {
            (&b, &a[0])
        };
        if linear.mapped && linear.unique {
            let c5 = &chim[0];
            let c3 = &chim[1];
            // Rescue reports the first chimeric alignment. If it cannot be
            // used, fall back to 5unique selection instead of discarding a
            // later uniquely mapped alignment.
            let mut can = c5.mapped && c5.unique;
            if can && c3.mapped {
                can = c3.chrom == linear.chrom && c3.strand != linear.strand;
                let delta = if linear.strand == '+' {
                    c3.pos - linear.pos
                } else {
                    linear.pos - c3.pos
                };
                can &= delta > 0 && delta + c3.dist5 + linear.dist5 <= max_molecule;
            }
            if can {
                if a.len() == 2 {
                    b[0].kind = 'R';
                } else {
                    a[0].kind = 'R';
                }
                rescued = true;
            }
        }
    }
    let pick = |v: Vec<Alignment>| -> Alignment {
        if rescued || v.len() == 1 {
            v[0].clone()
        } else {
            v.iter()
                .find(|x| x.mapped && x.unique)
                .unwrap_or(&v[0])
                .clone()
        }
    };
    let x = pick(a);
    let y = pick(b);
    if !matches!((x.kind, y.kind), ('U', 'U') | ('U', 'R') | ('R', 'U')) {
        return None;
    }
    Some((x, y))
}
#[derive(Default, Debug, serde::Serialize)]
pub struct Counts {
    pub input_pairs: u64,
    pub selected_pairs: u64,
}
fn emit<W: Write>(
    name: &str,
    a: &mut Vec<Alignment>,
    b: &mut Vec<Alignment>,
    out: &mut W,
    counts: &mut Counts,
    max_molecule: i64,
    max_gap: i64,
) -> Result<()> {
    if name.is_empty() {
        return Ok(());
    }
    counts.input_pairs += 1;
    if let Some((x, y)) = classify(std::mem::take(a), std::mem::take(b), max_molecule, max_gap) {
        writeln!(
            out,
            "{name}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}",
            x.chrom, x.pos, y.chrom, y.pos, x.strand, y.strand, x.kind, y.kind
        )?;
        counts.selected_pairs += 1;
    }
    Ok(())
}

pub(crate) fn write_pair<W: Write>(
    name: &str,
    a: &mut Vec<Alignment>,
    b: &mut Vec<Alignment>,
    out: &mut W,
    counts: &mut Counts,
    max_molecule: i64,
    max_gap: i64,
) -> Result<()> {
    emit(name, a, b, out, counts, max_molecule, max_gap)
}
pub fn convert<R: BufRead, W: Write>(
    input: R,
    mut out: W,
    min_mapq: u8,
    max_molecule: i64,
    max_gap: i64,
) -> Result<Counts> {
    if min_mapq == 0 || max_molecule < 0 || max_gap < 0 {
        bail!("invalid parsing thresholds")
    }
    writeln!(out,"## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type")?;
    let (mut a, mut b) = (Vec::new(), Vec::new());
    let mut name = String::new();
    let mut counts = Counts::default();
    let mut body = false;
    for (i, line) in input.lines().enumerate() {
        let line = line.with_context(|| format!("read SAM line {}", i + 1))?;
        if line.starts_with('@') && !body {
            continue;
        }
        body = true;
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 11 || f[0].is_empty() || f[0] == "*" {
            bail!("invalid SAM record on line {}", i + 1)
        }
        if f[0] != name {
            emit(
                &name,
                &mut a,
                &mut b,
                &mut out,
                &mut counts,
                max_molecule,
                max_gap,
            )?;
            name.clear();
            name.push_str(f[0]);
        }
        let (first, aln) =
            parse_alignment(&f, min_mapq).with_context(|| format!("SAM line {}", i + 1))?;
        if first {
            a.push(aln)
        } else {
            b.push(aln)
        }
    }
    emit(
        &name,
        &mut a,
        &mut b,
        &mut out,
        &mut counts,
        max_molecule,
        max_gap,
    )?;
    out.flush()?;
    Ok(counts)
}

#[cfg(test)]
mod structured_tests {
    use super::*;
    use bwa_mem3_rs::AlignedFields;

    fn sam_record(name: &str, flag: u16, pos: i64, mapq: u8, cigar: &str) -> String {
        format!("{name}\t{flag}\tchr1\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t*\t*\n")
    }

    #[test]
    fn clipped_single_alignments_remain_usable_on_either_mate_and_strand() {
        let mut sam = String::new();
        for (name, first_flag, first_cigar, second_flag, second_cigar) in [
            ("leading", 65, "21S79M", 129, "100M"),
            ("reverse", 81, "79M21S", 129, "100M"),
            ("mate", 65, "100M", 129, "30S70M"),
        ] {
            sam.push_str(&sam_record(name, first_flag, 100, 60, first_cigar));
            sam.push_str(&sam_record(name, second_flag, 1000, 60, second_cigar));
        }
        let mut output = Vec::new();
        let counts = convert(sam.as_bytes(), &mut output, 1, 750, 20).unwrap();
        let pairs = String::from_utf8(output).unwrap();
        assert_eq!(counts.input_pairs, 3);
        assert_eq!(counts.selected_pairs, 3);
        assert!(pairs.contains("leading\tchr1\t100\tchr1\t1000\t+\t+\tUU\n"));
        assert!(pairs.contains("reverse\tchr1\t178\tchr1\t1000\t-\t+\tUU\n"));
        assert!(pairs.contains("mate\tchr1\t100\tchr1\t1000\t+\t+\tUU\n"));
    }

    #[test]
    fn every_internal_gap_is_detected_without_a_leading_gap() {
        let alignment = |start, length| Alignment {
            chrom: "chr1".into(),
            pos: 100,
            strand: '+',
            mapped: true,
            unique: true,
            kind: 'U',
            dist5: start,
            read_span: length,
            query_start: start,
        };
        let normalized = normalize(
            vec![alignment(90, 10), alignment(30, 10), alignment(60, 10)],
            20,
        );
        // Gaps equal to the threshold do not become separate alignments.
        assert_eq!(
            normalized
                .iter()
                .map(|a| (a.kind, a.dist5, a.read_span))
                .collect::<Vec<_>>(),
            [('U', 30, 10), ('U', 60, 10), ('U', 90, 10)]
        );
        let normalized = normalize(
            vec![alignment(0, 10), alignment(40, 10), alignment(80, 10)],
            20,
        );
        assert_eq!(
            normalized
                .iter()
                .map(|a| (a.kind, a.dist5, a.read_span))
                .collect::<Vec<_>>(),
            [
                ('U', 0, 10),
                ('N', 10, 30),
                ('U', 40, 10),
                ('N', 50, 30),
                ('U', 80, 10),
            ]
        );
    }

    #[test]
    fn unusable_first_chimeric_alignment_falls_back_to_unique_selection() {
        let mut sam = String::new();
        // The first alignment is non-unique, while the second is unique.
        // A rescue must not choose the non-unique first alignment.
        sam.push_str(&sam_record("fallback", 65, 100, 0, "40M60S"));
        sam.push_str(&sam_record("fallback", 2113, 300, 60, "40S60M"));
        sam.push_str(&sam_record("fallback", 129, 1000, 60, "100M"));
        // The same fallback applies when the chimeric read is mate 2.
        sam.push_str(&sam_record("mate-fallback", 65, 1000, 60, "100M"));
        sam.push_str(&sam_record("mate-fallback", 129, 100, 0, "40M60S"));
        sam.push_str(&sam_record("mate-fallback", 2177, 300, 60, "40S60M"));
        // A valid rescue must continue to produce U/R.
        sam.push_str(&sam_record("rescue", 65, 5000, 60, "40M60S"));
        sam.push_str(&sam_record("rescue", 2129, 1100, 60, "60M40S"));
        sam.push_str(&sam_record("rescue", 129, 1000, 60, "100M"));
        // A non-unique second alignment still falls back to 5unique.
        sam.push_str(&sam_record("second-nonunique", 65, 5000, 60, "40M60S"));
        sam.push_str(&sam_record("second-nonunique", 2113, 1100, 0, "40S60M"));
        sam.push_str(&sam_record("second-nonunique", 129, 1000, 60, "100M"));
        let mut output = Vec::new();
        let counts = convert(sam.as_bytes(), &mut output, 1, 750, 20).unwrap();
        let pairs = String::from_utf8(output).unwrap();
        assert_eq!(counts.input_pairs, 4);
        assert_eq!(counts.selected_pairs, 4);
        assert!(pairs.contains("fallback\tchr1\t300\tchr1\t1000\t+\t+\tUU\n"));
        assert!(pairs.contains("mate-fallback\tchr1\t1000\tchr1\t300\t+\t+\tUU\n"));
        assert!(pairs.contains("rescue\tchr1\t5000\tchr1\t1000\t+\t+\tUR\n"));
        assert!(pairs.contains("second-nonunique\tchr1\t5000\tchr1\t1000\t+\t+\tUU\n"));
    }

    #[test]
    fn structured_fields_match_sam_for_clips_indels_and_mapping_quality() {
        let cases: &[(u16, u8, &str, &[u32])] = &[
            (65, 60, "100M", &[(100 << 4)]),
            (65, 60, "21S79M", &[(21 << 4) | 4, (79 << 4)]),
            (81, 60, "79M21S", &[(79 << 4), (21 << 4) | 4]),
            (
                81,
                60,
                "5S45M2I3M4D45M5S",
                &[
                    (5 << 4) | 4,
                    (45 << 4),
                    (2 << 4) | 1,
                    (3 << 4),
                    (4 << 4) | 2,
                    (45 << 4),
                    (5 << 4) | 4,
                ],
            ),
            (65, 60, "10H90M", &[(10 << 4) | 5, (90 << 4)]),
            (81, 0, "100M", &[(100 << 4)]),
            (69, 0, "*", &[]),
        ];
        for &(flag, mapq, cigar, encoded) in cases {
            let fields = AlignedFields {
                tid: 0,
                pos: 99,
                next_tid: -1,
                next_pos: -1,
                tlen: 0,
                flag,
                bin: 0,
                mapq,
                query_start: 0,
                query_end: 100,
                cigar: encoded,
                nm: None,
                md: None,
                mc: None,
                mq: None,
                score: None,
                sub: None,
                rg: None,
                sa: None,
                pa: None,
                xa: None,
                hn: None,
                meth: None,
            };
            let flag_text = flag.to_string();
            let mapq_text = mapq.to_string();
            let sam = [
                "read", &flag_text, "chr1", "100", &mapq_text, cigar, "*", "0", "0", "*", "*",
            ];
            let (_, expected) = parse_alignment(&sam, 1).unwrap();
            assert_eq!(
                from_fields_with_chrom(&fields, "chr1", 1).unwrap(),
                expected,
                "{cigar} flag={flag}"
            );
        }
    }
}
