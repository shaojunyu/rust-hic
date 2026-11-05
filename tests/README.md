# Core tests

The suite exercises the public commands with a deterministic miniature dataset:

1. **FASTQ index:** verify the default and three-chunk indexes, checking that
   both mates cover every original read exactly once and chunk IDs match; also
   check inputs containing concatenated gzip members.
2. **Mapping:** build a BWA-MEM2 reference index, compare its combined pairs
   output to built-in bwa-mem3, and verify single-chunk mapping.
3. **Pairs:** sort, filter and deduplicate two pairs files, then verify the
   contact counts and default resolutions in the resulting `.mcool`.

Requirements: a built `target/release/rust-hic`, Python 3.9 or later,
`h5py`, and `bwa-mem2`. The BWA executable can be on `PATH`, in the local
`envs/runhic/bin` environment, or specified via `RUST_HIC_BWA_MEM2`.
The binary can be overridden using `RUST_HIC_BIN`.

```sh
python3 tests/test_core.py
```

`tests/test_chromap_pairing.py` checks combined multi-chunk output, valid input,
an internal mate-ID mismatch, and a missing mate record. It also checks consumers
that read one complete mate before the other. It requires a built binary and Chromap on
`PATH` or in `RUST_HIC_CHROMAP`:

```sh
python3 tests/test_chromap_pairing.py
```

`tests/test_dedup_pairtools.py` compares deduplicated rows against pairtools
1.1.3's default scipy backend, including a 10,000-row batch boundary. Set
`RUST_HIC_PAIRTOOLS` if `pairtools` is not on `PATH`.

`tests/test_stage_processes.py` checks nested sorting and merging processes
using multiple temporary runs and parallel merge ranges.

`tests/test_example_config.py` checks the embedded example and reports, and runs
it with pairs plus Arima digestion, without restriction filtering, and with
built-in BWA-MEM3 FASTQ input. It also runs the minimal example using defaults
and checks report regeneration and detailed statistics for older outputs. The
BWA-MEM3 case requires `bwa-mem2` to build
the test reference index and is skipped if it is unavailable.

`tests/test_workflow_overlap.py` checks that chunk sorting and restriction
filtering overlap mapping, match sequential processing and keep deduplication
global. It also checks that a failed chunk stops the workflow.

`tests/test_bam_retention.py` checks optional unsorted BAM retention for
BWA-MEM3, BWA-MEM2 and Chromap across concurrent chunks, compares pairs and
matrix outputs, and checks Chromap's emitted SAM records against the BAM. It
requires `bwa-mem2`, `chromap` and `samtools`.

`data/` contains two 25 kb chromosomes, 1,024 paired reads (including four
unmapped pairs), chromosome sizes, restriction fragments, and two pairs files
with known duplicates and one same-fragment contact. Regenerate them with:

```sh
python3 tests/data/generate.py
```

Tests create all indexes and outputs in a temporary directory and delete it
afterward. The reference comparison checks this dataset only; it does not
assert that the two aligners agree on all biological inputs.

`tests/test_preflight.py` exercises `run --check`, invalid input and output
paths, chromosome sizes, fragment BEDs, sorted pairs headers, gzip signatures,
reference index metadata, and reusable-index fingerprint failures. It requires
only Python 3 and a built executable; set `RUST_HIC_BIN` to select the binary.
