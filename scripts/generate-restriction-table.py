#!/usr/bin/env python3
"""Generate src/workflow/restriction.rs from Biopython's Bio.Restriction.

This is a development-time generator. rust-hic does not depend on Biopython at
build or run time; the generated table is checked in.

Only enzymes with integer top-strand cut and overhang offsets and a valid
IUPAC recognition site are included. Biopython's `search` returns one Watson
cut per site occurrence, including for non-palindromic enzymes.
"""

from pathlib import Path
import sys

from Bio import Restriction as R
import Bio

IUPAC = set("ACGTURYSWKMBDHVN")
OUT = Path(__file__).resolve().parent.parent / "src" / "workflow" / "restriction.rs"

entries = []
skipped = []
for enzyme in R.AllEnzymes:
    name = enzyme.__name__
    if not isinstance(enzyme.fst5, int) or not isinstance(enzyme.ovhg, int):
        skipped.append((name, "non-integer cut offset"))
        continue
    site = str(enzyme.site).upper()
    if not site or any(base not in IUPAC for base in site):
        skipped.append((name, "unsupported recognition site"))
        continue
    entries.append((name, site, enzyme.fst5, enzyme.ovhg))

entries.sort()

if not entries:
    sys.exit("no enzymes generated")

lines = [
    "//! Restriction enzyme definitions for in-silico digestion.",
    "//!",
    f"//! Generated from Biopython {Bio.__version__} `Bio.Restriction` with",
    "//! `scripts/generate-restriction-table.py`. Do not edit by hand.",
    "//!",
    "//! Only enzymes with integer top-strand cut and overhang offsets and a",
    "//! valid IUPAC recognition site are included. `fst5` is Biopython's 1-based",
    "//! top-strand cut offset; `ovhg` is its overhang offset. The scanner",
    "//! converts these to runHiC's fragment BED boundaries.",
    "",
    "#[derive(Clone, Copy, Debug)]",
    "pub(crate) struct Enzyme {",
    "    pub(crate) name: &'static str,",
    "    pub(crate) site: &'static str,",
    "    pub(crate) fst5: i32,",
    "    pub(crate) ovhg: i32,",
    "}",
    "",
    "pub(crate) const ENZYMES: &[Enzyme] = &[",
]
for name, site, fst5, ovhg in entries:
    lines.append(
        f'    Enzyme {{ name: "{name}", site: "{site}", fst5: {fst5}, ovhg: {ovhg} }},'
    )
lines.extend(
    [
        "];",
        "",
        "pub(crate) fn lookup(name: &str) -> Option<&'static Enzyme> {",
        "    ENZYMES",
        "        .iter()",
        "        .find(|enzyme| enzyme.name.eq_ignore_ascii_case(name))",
        "}",
        "",
        "/// Parse `Arima`, a single enzyme, or a `-`-separated cocktail.",
        "pub(crate) fn parse_spec(spec: &str) -> Result<Vec<&'static Enzyme>, String> {",
        "    let spec = spec.trim();",
        "    if spec.is_empty() {",
        '        return Err("digest enzyme must not be empty".into());',
        "    }",
        "    if spec.eq_ignore_ascii_case(\"arima\") {",
        '        return Ok(vec![lookup("MboI").unwrap(), lookup("HinfI").unwrap()]);',
        "    }",
        "    let mut enzymes = Vec::new();",
        "    for part in spec.split('-') {",
        "        let part = part.trim();",
        "        if part.is_empty() {",
        '            return Err(format!("invalid digest enzyme list: {spec}"));',
        "        }",
        "        let enzyme = lookup(part)",
        '            .ok_or_else(|| format!("unknown or unsupported restriction enzyme \'{part}\'"))?;',
        "        if !enzymes",
        "            .iter()",
        "            .any(|existing: &&Enzyme| existing.name.eq_ignore_ascii_case(enzyme.name))",
        "        {",
        "            enzymes.push(enzyme);",
        "        }",
        "    }",
        "    Ok(enzymes)",
        "}",
        "",
    ]
)

OUT.write_text("\n".join(lines))
print(f"wrote {len(entries)} enzymes to {OUT}; skipped {len(skipped)}")
for name, reason in skipped[:10]:
    print(f"  skipped {name}: {reason}")
