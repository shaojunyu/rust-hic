//! Supported pairs column layouts shared by every processing stage.
use anyhow::{bail, Context, Result};
use std::io::BufRead;

pub(crate) const BASE_COLUMNS: [&str; 8] = [
    "readID",
    "chrom1",
    "pos1",
    "chrom2",
    "pos2",
    "strand1",
    "strand2",
    "pair_type",
];
pub(crate) const FRAGMENT_COLUMNS: [&str; 6] = [
    "rfrag1",
    "rfrag_start1",
    "rfrag_end1",
    "rfrag2",
    "rfrag_start2",
    "rfrag_end2",
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Format {
    Standard,
    Fragments,
}
impl Format {
    pub fn from_columns(line: &str) -> Result<Self> {
        let columns: Vec<_> = line
            .strip_prefix("#columns:")
            .context("missing #columns header")?
            .split_whitespace()
            .collect();
        if columns == BASE_COLUMNS {
            Ok(Self::Standard)
        } else if columns.len() == 14
            && columns[..8] == BASE_COLUMNS
            && columns[8..] == FRAGMENT_COLUMNS
        {
            Ok(Self::Fragments)
        } else {
            bail!("unsupported pairs columns: expected standard 8 columns or 14 columns with restriction-fragment annotations")
        }
    }
    pub fn len(self) -> usize {
        match self {
            Self::Standard => 8,
            Self::Fragments => 14,
        }
    }
    pub fn columns_line(self) -> String {
        let mut s = format!("#columns: {}", BASE_COLUMNS.join(" "));
        if self == Self::Fragments {
            s.push(' ');
            s.push_str(&FRAGMENT_COLUMNS.join(" "));
        }
        s.push('\n');
        s
    }
    pub fn fields(self, line: &str) -> Result<Vec<&str>> {
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() != self.len() {
            bail!(
                "pairs row has {} fields but #columns declares {}",
                fields.len(),
                self.len()
            )
        }
        Ok(fields)
    }
}

/// Read a pairs header and leave the stream at the beginning of its body.
pub(crate) fn read_header<R: BufRead>(reader: &mut R) -> Result<(Vec<String>, Format)> {
    let mut header = Vec::new();
    let mut format = None;
    while reader.fill_buf()?.first() == Some(&b'#') {
        let mut line = String::new();
        reader.read_line(&mut line)?;
        if line.starts_with("#columns:") {
            if format.is_some() {
                bail!("duplicate #columns header")
            }
            format = Some(Format::from_columns(&line)?);
        }
        header.push(line);
    }
    Ok((header, format.context("missing #columns header")?))
}
