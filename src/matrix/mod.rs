mod attr;
pub(crate) mod bin_pairs;
pub(crate) mod multiresolution;

/// Internal setting passed to the matrix balancing stage.
pub(crate) const BALANCE_MAX_ITERS_ENV: &str = "RUST_HIC_BALANCE_MAX_ITERS";
const SORT_CHROMOSOMES_ENV: &str = "RUST_HIC_SORT_CHROMOSOMES";

pub(crate) fn chromosome_options(cmd: &mut std::process::Command, sort: bool) {
    cmd.env(SORT_CHROMOSOMES_ENV, if sort { "1" } else { "0" });
}

use anyhow::{bail, Context, Result};
use std::{
    cmp::Ordering,
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

struct Chromosome {
    pub name: String,
    pub len: u64,
}

fn canonical_rank(name: &str) -> Option<(u8, u32)> {
    let plain = name.strip_prefix("chr").unwrap_or(name);
    match plain {
        "X" => Some((1, 0)),
        "Y" => Some((2, 0)),
        "Z" => Some((3, 0)),
        "W" => Some((4, 0)),
        "M" | "MT" => Some((5, 0)),
        _ => {
            let n = plain.parse::<u32>().ok()?;
            (n > 0).then_some((0, n))
        }
    }
}

fn natural_cmp(a: &str, b: &str) -> Ordering {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i].is_ascii_digit() && b[j].is_ascii_digit() {
            let (start_i, start_j) = (i, j);
            while i < a.len() && a[i].is_ascii_digit() {
                i += 1;
            }
            while j < b.len() && b[j].is_ascii_digit() {
                j += 1;
            }
            let (mut x, mut y) = (start_i, start_j);
            while x < i && a[x] == b'0' {
                x += 1;
            }
            while y < j && b[y] == b'0' {
                y += 1;
            }
            let order = (i - x)
                .cmp(&(j - y))
                .then_with(|| a[x..i].cmp(&b[y..j]))
                .then_with(|| (i - start_i).cmp(&(j - start_j)));
            if order != Ordering::Equal {
                return order;
            }
        } else {
            let order = a[i].cmp(&b[j]);
            if order != Ordering::Equal {
                return order;
            }
            i += 1;
            j += 1;
        }
    }
    (a.len() - i).cmp(&(b.len() - j))
}

fn read_chromosomes(path: &Path, sort: bool) -> Result<Vec<Chromosome>> {
    let mut all = Vec::new();
    let mut seen = HashSet::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.is_empty() {
            continue;
        }
        if fields.len() != 2 {
            bail!("invalid chromosome sizes")
        }
        let len: u64 = fields[1].parse().context("invalid chromosome length")?;
        if len == 0 || len > i32::MAX as u64 || !seen.insert(fields[0].to_owned()) {
            bail!("invalid chromosome")
        }
        all.push(Chromosome {
            name: fields[0].to_owned(),
            len,
        });
    }
    if all.is_empty() {
        bail!("empty chromosome sizes")
    }
    if sort {
        all.sort_by(|a, b| {
            canonical_rank(&a.name)
                .unwrap_or((u8::MAX, 0))
                .cmp(&canonical_rank(&b.name).unwrap_or((u8::MAX, 0)))
                .then_with(|| natural_cmp(&a.name, &b.name))
        });
    }
    Ok(all)
}

#[cfg(test)]
mod tests {
    use super::{canonical_rank, natural_cmp};
    use std::cmp::Ordering;

    #[test]
    fn natural_order_compares_numbers_without_overflow() {
        assert_eq!(natural_cmp("chrUn_2", "chrUn_10"), Ordering::Less);
        assert_eq!(
            natural_cmp("scaffold_9", "scaffold_999999999999999999999"),
            Ordering::Less
        );
        assert_eq!(natural_cmp("chr_2", "chr_02"), Ordering::Less);
    }

    #[test]
    fn chromosome_rank_depends_only_on_name() {
        assert_eq!(canonical_rank("chr19"), Some((0, 19)));
        assert_eq!(canonical_rank("chr22"), Some((0, 22)));
        assert_eq!(canonical_rank("chr02"), Some((0, 2)));
        assert_eq!(canonical_rank("23"), Some((0, 23)));
        assert_eq!(canonical_rank("X"), Some((1, 0)));
        assert_eq!(canonical_rank("chrY"), Some((2, 0)));
        assert_eq!(canonical_rank("MT"), Some((5, 0)));
        assert_eq!(canonical_rank("chrM"), Some((5, 0)));
        assert_eq!(canonical_rank("chrUn_1"), None);
    }
}
