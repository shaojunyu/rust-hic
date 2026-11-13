"""Locate artifacts in a grouped rust-hic result directory.

Tests use globs instead of recomputing the output prefix, so they stay valid
for both explicit and auto-derived prefixes.
"""

from pathlib import Path


def _one(run_dir, pattern):
    matches = sorted(Path(run_dir).glob(pattern))
    if len(matches) != 1:
        raise AssertionError(f"expected one {pattern} in {run_dir}, found {matches}")
    return matches[0]


def mcool(run_dir):
    return _one(run_dir, "matrix/*.mcool")


def pairs(run_dir):
    return _one(run_dir, "pairs/*.filtered-nodups.pairs")


def compressed_pairs(run_dir):
    return _one(run_dir, "pairs/*.filtered-nodups.pairs.gz")


def pairix(run_dir):
    return _one(run_dir, "pairs/*.filtered-nodups.pairs.gz.px2")


def bam(run_dir):
    return _one(run_dir, "alignments/*.bam")


def qc(run_dir, name):
    return _one(run_dir, f"qc/*.{name}")


def logs(run_dir):
    return Path(run_dir) / "logs"


def event_log(run_dir):
    return logs(run_dir) / "events.jsonl"


def has_mcool(run_dir):
    return any(Path(run_dir).glob("matrix/*.mcool"))


def has_pairs(run_dir):
    return any(Path(run_dir).glob("pairs/*.filtered-nodups.pairs"))


def has_bam(run_dir):
    return any(Path(run_dir).glob("alignments/*.bam"))
