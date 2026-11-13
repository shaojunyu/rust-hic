# Usage

This guide describes the version in the same source checkout or release
archive. When using a downloaded binary, consult the documentation bundled
with it; commands and options can differ between versions.

`rust-hic run CONFIG.yaml` runs locally by default. Explicit `--mode local`
uses the same single-machine workflow; `--mode slurm` submits it to a cluster.
For Slurm configuration and resource requests, see [Slurm execution](slurm.md).

The workflow indexes paired gzip FASTQs, aligns chunks, sorts,
merges, filters and deduplicates contacts, and produces a multi-resolution
`.mcool` file. FASTQ alignment supports built-in BWA-MEM3 and external
Chromap, BWA-MEM or BWA-MEM2 executables.
In local mode, as each FASTQ chunk finishes mapping, `run` starts its sorting
and optional restriction filtering while other chunks continue mapping. It then merges
all chunks and deduplicates globally, so duplicates spanning chunks are found.
Mapping and chunk preparation can run at the same time; their configured
thread and job counts are not a shared CPU limit.
The pipeline is experimental and does not produce the complete runHiC statistics cache.

Download the Linux executable as described in the [README](../README.md). The
examples below assume it is named `rust-hic` in the current directory. See
[building from source](building.md) for development and release builds.

## Commands

| Command | Purpose |
| --- | --- |
| `example-config` | Print the complete YAML example; use `--minimal` for required settings only. |
| `run` | Run locally or submit a FASTQ workflow to Slurm. |
| `run status` / `run cancel` | Inspect or cancel a submitted Slurm workflow. |
| `run report` | Generate HTML and JSON reports for a completed run. |
| `index` | Build an index for paired gzip FASTQs. |
| `map` | Map indexed FASTQs to one unsorted pairs file. |
| `clean` | Sort, optionally filter, and deduplicate pairs. |
| `matrix` | Build a multi-resolution matrix from cleaned pairs. |
| `version` | Show the executable version. |

Run `./rust-hic COMMAND --help` for arguments and options. Pairs inputs to
`run`, `clean` and `matrix` must be uncompressed.

## Prepare a reference

The default built-in BWA-MEM3 aligner needs a reference FASTA, a two-column
chromosome sizes file and an existing BWA-MEM2-compatible reference index. It
can reuse an existing index; it does not build one itself. To prepare a new
index and sizes file, for example:

```sh
bwa-mem2 index hg38.fa
samtools faidx hg38.fa
cut -f 1,2 hg38.fa.fai > hg38.chrom.sizes
```

`bwa-mem2` and `samtools` are external programs needed for these preparation
commands. They are not needed to run the built-in aligner after the reference
files have been created. In the example configuration, set `reference` and
`aligner_index` to `hg38.fa`, and `chromsizes` to `hg38.chrom.sizes`. Use
matching reference, index and chromosome files for the chosen assembly.

## Configure and run

Create a minimal configuration and edit its input paths:

```sh
./rust-hic example-config --minimal > config.yaml
```

It contains only the settings required for FASTQ input. Other settings use their
defaults, including the built-in BWA-MEM3 aligner. Its reference index must be
created beforehand; a BWA-MEM2-compatible index can be used. Restriction
fragment filtering is disabled unless explicitly configured.

For an annotated example of the available YAML settings, use:

```sh
./rust-hic example-config > config.yaml
```

The files are also available as [minimal](example.minimal.config.yaml) and
[complete](example.config.yaml) examples. The complete example documents
alternatives for pairs input and optional fragment filtering. Relative paths are
resolved against the YAML file's directory. Run
`./rust-hic run config.yaml --check` before starting a local pipeline. This
checks settings, inputs, available aligner and index files, and whether the
output parent and scratch directories can create temporary files. It reads the
chromosome sizes and fragment BED, checks pairs headers and gzip signatures,
and verifies fingerprints and boundary IDs when reusing a FASTQ index. It does
not fully validate FASTQ or pairs records, verify all reference-index contents,
or guarantee sufficient disk space. The check does not create the result directory.

For BWA, preflight also compares shared chromosome names and lengths with the
index's `.ann` file. If the reference FASTA has a `.fai` index, it compares
that too for any aligner. Different chromosome subsets are allowed, but at
least one name must match. These checks do not compare the FASTA sequence with
the aligner index; the reference and its indexes must be kept consistent.

The main FASTQ options are:

| Setting | Default | Purpose |
| --- | ---: | --- |
| `fastq.aligner` | `bwa-mem3` | Built-in `bwa-mem3` or external `chromap`, `bwa-mem`, `bwa-mem2`. |
| `fastq.chunks` | automatic | Explicit chunk count; otherwise use `chunk_size_bytes`. |
| `fastq.chunk_size_bytes` | 2147483648 | Target size in compressed bytes for estimating chunk count. |
| `fastq.index_threads_per_mate` | 8 | Indexing threads for each of R1 and R2. |
| `fastq.index_search_mib` | 128 | Maximum R2 boundary search distance in MiB. |
| `fastq.jobs` | 4 | Maximum concurrently mapped chunks. |
| `fastq.threads_per_job` | 8 | Aligner threads per chunk. |
| `fastq.index_dir` | unset | Reuse a previously created FASTQ index. |
| `fastq.expected_pairs` | unset | Check total pair count after mapping. |

Resource limits are per stage, not a shared process-wide CPU limit. An explicit
setting in an existing YAML file overrides these defaults.

```sh
./rust-hic version
./rust-hic help run
./rust-hic run config.yaml
```

## Pipeline outputs

A successful `run` writes a grouped result directory. A `prefix` labels every
data and QC file. It defaults to the R1 filename stem with a trailing mate token
removed (`GM12878_100M_R1.fastq.gz` becomes `GM12878_100M`), or to the first
pairs input's stem for a pairs run. Set `output.prefix` to override it:

```yaml
output:
  prefix: GM12878_100M
```

```text
result/
├── matrix/<prefix>.mcool
├── pairs/<prefix>.filtered-nodups.pairs              # --keep-pairs
├── pairs/<prefix>.filtered-nodups.pairs.gz + .px2    # compressed_pairs
├── alignments/<prefix>.bam                           # --keep-bam
├── qc/<prefix>.report.html
├── qc/<prefix>.report.json
├── qc/<prefix>.contact-stats.json
├── qc/<prefix>.filter-stats.json                     # restriction filtering
├── logs/events.jsonl and stage logs
├── logs/slurm/                                       # Slurm mode
├── config.resolved.yaml
├── status.json
├── MANIFEST.tsv
├── work-dir.txt                                      # retained on failure
└── failure.json                                      # on failure
```

| File | Contents |
| --- | --- |
| `matrix/<prefix>.mcool` | Multi-resolution contact matrix. |
| `pairs/<prefix>.filtered-nodups.pairs` | Sorted, filtered and deduplicated contacts, when `--keep-pairs` or `keep_pairs: true` is set. |
| `pairs/<prefix>.filtered-nodups.pairs.gz` and `.px2` | Optional compressed pairs and pairix index when `compressed_pairs: true`. |
| `alignments/<prefix>.bam` | Aligner-emitted records in one unsorted BAM, when `--keep-bam` or `keep_bam: true` is set. |
| `qc/...` | Run report and machine-readable statistics. |
| `logs/` | Progress events and stage logs; `logs/slurm/` holds scheduler logs. |
| `status.json` | Completion status and elapsed time. `failure.json` records a failed run. |
| `MANIFEST.tsv` | Size and path of every file in the result, relative to its root. |

`output.prefix` must be a single filename component: ASCII letters, digits,
`.`, `_` and `-`, at most 128 characters, and it cannot start or end with `.`
or contain `..`. A derived prefix that fails these rules stops preflight with a
message asking for an explicit `output.prefix`.

Older results used a flat directory (`contacts.mcool`,
`filtered-nodups.pairs`, `report.json`, and stage logs at the top level).
`run report` still reads those results; new runs always use the grouped layout.

The mcool stores its chromosome, pixel, bin and index datasets with gzip compression.
By default, chromosomes are sorted in natural numeric order. Standard
chromosomes come first, followed by other contigs (for example, `chrUn_2`
before `chrUn_10`). All contigs in the chromosome sizes file are included.
Sorting uses chromosome names only: numbered chromosomes are followed by X,
Y, Z, W and M/MT. Names with or without `chr` are supported, and chromosome
names are not changed; additional contigs follow in natural order. Assembly
names do not affect the sorting.
Set `sort_chromosomes: false` to use the original chromosome sizes file order.

By default, neither BAM nor uncompressed pairs are retained. For example,
`rust-hic run config.yaml --keep-bam --keep-pairs` retains both. BAM contains
the records emitted by the aligner, before pairs filtering and deduplication;
it is not coordinate sorted and has no index. Chromap itself may omit reads
from its SAM output. Because Chromap emits either pairs or SAM, retaining its
BAM runs a second alignment pass and adds substantial time and I/O. With BAM
disabled, no BAM file or BAM output stream is created. `--keep-bam` requires
FASTQ input.

## Run report

Open `qc/<prefix>.report.html` in a web browser after a successful `run`. It is self-contained
and does not need an internet connection. It summarizes input read-pair counts
when available, contacts after each processing stage, duplicates, cis/trans and
distance distributions, strand and pair-type composition, matrix balance
status at each resolution, chunk details, and stage timings. Percentages name
their denominator. Counts after aligner selection are not necessarily a mapping
rate, and duplicates removed inside an aligner are not included in the pipeline
deduplication count. In local mode, alignment and chunk preparation overlap.
In Slurm mode, chunk preparation follows the alignment array and total elapsed
time includes scheduler waits.

`qc/<prefix>.report.json` contains the same statistics for downstream tools. Missing or
unavailable values are `null`, not zero. A resolution with no usable balance
weights is reported separately from converged or non-converged resolutions.
The report summarizes observable results; it does not assign a pass/fail
quality grade or establish biological validity.

To generate a report for an earlier successful run, choose an existing empty
directory for the two report files:

```sh
mkdir report-output
./rust-hic run report existing-result --output report-output
```

With `--output`, the two files are written directly into the chosen directory as
`report.html` and `report.json`. Without it, regeneration writes
`qc/<prefix>.report.html` and `qc/<prefix>.report.json` in the result directory.

Older runs did not record distance, strand and pair-type distributions; these
are marked unavailable. Add `--recompute-contact-stats` to read the entire
final uncompressed pairs file once (which must have been retained) and calculate those distributions without
modifying the original result. This can take time for large datasets. The
command does not overwrite an existing report. Report creation is part of a
successful `run`; an error leaves the output and `failure.json` for inspection.

## Map FASTQ chunks

To map an existing FASTQ index to one pairs file without running the rest of
the pipeline:

```sh
./rust-hic map fastq-index \
  -r reference.fa -x reference.fa \
  --jobs 4 --threads-per-job 8 -o mapped.pairs
```

`map` reads the paired FASTQ paths from the index manifest and checks that the
inputs still match. The reference FASTA and aligner index must already exist;
external aligners must also be installed. For BWA, `-x` is the prefix used to
build its index. Each selected chunk is first mapped to a temporary pairs
file. After all chunks succeed, their records are combined in chunk order
into the file specified by `-o`, with one header. The output is uncompressed
and unsorted, before filtering or deduplication. `clean` can sort and deduplicate
the combined output. Existing output files are never overwritten; failed
runs retain temporary files and logs for inspection.
Pass extra aligner options with e.g. `--aligner-arg "-k 19 -T 30"`;
the argument string is split using shell-style quoting without invoking a shell.
The option may be repeated, and options for program-managed inputs, outputs,
or threads are rejected. In a `run` configuration, use `fastq.aligner_args`
as a YAML list of individual arguments.
Use `--chunk N` to map only one chunk, numbered from 0 (for example,
`--chunk 0` selects the first chunk). Without it, all chunks are mapped.
Run `rust-hic map --help` for options.

## Clean pairs

Use `clean` to sort and merge pairs files, apply optional restriction-fragment
filtering, and remove duplicate contacts across all inputs:

```sh
./rust-hic clean sample.pairs \
  --chromsizes hg38.chrom.sizes --fragments hg38.fragments.bed \
  --assembly hg38 --threads 4 -o sample.clean.pairs
```

The output is a sorted, uncompressed pairs file. Without `--fragments`,
same-fragment filtering is skipped. The output must not already exist; stage
logs are retained in a temporary workspace only on failure. `--threads` is a
limit per stage, not a global CPU limit. Run `rust-hic clean --help` for
options.

Deduplication defaults to `--max-mismatch 3`, matching pairtools 1.1.3's
default `--method max`: both positions may differ by at most 3 bp, with matching
chromosomes and strands. Duplicate clusters are transitive within each batch.
Use `clean --max-mismatch 1` (or `dedup_max_mismatch: 1` in YAML) for runHiC's
explicit setting. Like pairtools' default scipy backend, batches contain 10,000
records and carry over up to 100 retained records, so unusually dense clusters
at batch boundaries can be missed. This change from the previous fixed value of
1 can remove more contacts from the same data.

Pairs inputs to `run`, `clean` and `matrix` must be uncompressed `.pairs`
files. `compressed_pairs` creates an additional compressed output for later use
by other tools; it does not enable compressed input.

Pairs processing supports two column layouts: the standard eight columns
`readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type`, optionally followed
by `rfrag1 rfrag_start1 rfrag_end1 rfrag2 rfrag_start2 rfrag_end2` (14 columns
in total). The `#columns:` header and every row must match, and all inputs to a
merge must use the same headers. Other layouts are rejected rather than
silently discarding columns. Existing fragment annotations are preserved,
including when a pair is flipped during sorting. Use `--fragments` only with
eight-column input; 14-column input already has fragment annotations and must
not be filtered again with this option. These rules also apply to `run` with
`input_sorted: true` and to `matrix` input.

## Build a matrix from cleaned pairs

After `clean`, create a multi-resolution `.mcool` from the resulting pairs file:

```sh
./rust-hic matrix sample.clean.pairs \
  --chromsizes hg38.chrom.sizes --assembly hg38 --threads 4 -o sample.mcool
```

`matrix` bins and balances the input records. It does not sort or deduplicate
pairs; use `clean` first. Use `--no-sort-chromosomes` to preserve the input
chromosome order.
Set `--resolutions` to a
comma-separated list of bin sizes in bases. Each must be a positive multiple of
the first, which must be the smallest. The complete example configuration lists
the default resolutions; run `rust-hic matrix --help` for options.

The output must not exist. Intermediate files are created next to the output
and removed after success. If processing fails, the workspace and stage logs
are retained and their location is printed.

## BWA alignment to pairs conversion

For BWA-MEM and BWA-MEM2, `map` converts the aligner's SAM output to pairs using
the `5unique` walk policy. These conversion options are not passed to BWA:

- `--min-mapq N`: MAPQ threshold for classifying an alignment as unique
  (default 1).
- `--max-molecule-size N`: molecule size threshold in bases for resolving walks
  (default 750).
- `--max-inter-align-gap N`: gap threshold in bases when grouping alignments
  (default 20).

For example, use `--min-mapq 10` to change the pairs conversion threshold; use
`--aligner-arg "-T 30"` to pass `-T 30` to BWA itself. These conversion options
do not apply to Chromap, which produces pairs directly; non-default SAM parsing
options are rejected with Chromap. In a `run` configuration,
set the corresponding `fastq.sam` fields `min_mapq`, `max_molecule_size`, and
`max_inter_align_gap`.

## Compressed pairs

With `compressed_pairs: true`, the pipeline also writes
`pairs/<prefix>.filtered-nodups.pairs.gz` and its `.px2` pairix index in Rust, using the
specified `compression_threads` (default 4). Compression reads the temporary
plain pairs after matrix generation, so its time is included in the pipeline
total and it requires space for both files.

The pairix writer accepts mapped, tab-delimited 8- or 14-column pairs sorted by
chromosome pair and first position, with coordinates from 1 through 2^30. It
rejects a repeated chromosome-pair block and unsupported input; failed runs
can leave partial files.

## Choosing an aligner

The `map` and `run` commands default to the built-in `bwa-mem3` aligner, which
needs no executable and shares one reference index across concurrent chunks.
Set `-x` (or `fastq.aligner_index` in YAML) to an existing BWA-MEM2 or BWA-MEM3
reference index prefix; a BWA-MEM2 index can be reused without rebuilding.
One thread reads each chunk, and
`--threads-per-job` threads align batches; structured alignment fields are
converted directly to pairs. The implementation uses 4096-pair batches and
BWA `-S -P` Hi-C settings. Differences from BWA-MEM2 may occur. The built-in
aligner accepts `--aligner-arg` options `-k`, `-w`, `-T`, `-A`, `-B` and `-c`;
other aligner arguments are rejected. On x86-64 it requires SSE4.2.

For example, map one chunk using the built-in aligner:

```sh
./rust-hic map fastq-index --chunk 0 --aligner bwa-mem3 \
  -r reference.fa -x reference.fa --threads-per-job 16 -o mapped.pairs
```

This creates `mapped.pairs`. The reference index must already exist;
the built-in aligner does not build it.

For external BWA, set `aligner` to `bwa-mem` or `bwa-mem2`, choose the corresponding
executable, and set `aligner_index` to the prefix used by `bwa index` or
`bwa-mem2 index`. All BWA paths use a `5unique`-style walk policy, retaining
`UU`, `UR` and `RU` pairs. Compared with pairtools 1.1.3, clipping before the
first alignment does not create a gap, all internal gaps are checked, and a
walk whose first alignment is unusable falls back to unique-alignment
selection rather than rescue. These cases can produce different pairs. The SAM
parser accepts BWA's M/I/D/S/H CIGAR operations; other operations fail
explicitly. The compatibility checks cover real
BWA-MEM 0.7.19 and BWA-MEM2 2.2.1 runs on synthetic reads, and randomized and
targeted SAM comparisons with pairtools 1.1.3. They do not prove equivalence
for every possible SAM or biological dataset.
The preflight checks the standard BWA-MEM or BWA-MEM2 2.2.1 index files
before the FASTQ scan; alternative index layouts are not currently supported.

`aligner_executable` may be omitted when the corresponding `chromap`, `bwa`
or `bwa-mem2` executable is on `PATH`; an explicit path takes precedence.

For each BWA version, a test using 100,000 real U937_WT read pairs and a
small synthetic reference built from 2,000 reads per mate produced 4,077
selected pairs identical to pairtools on the same chunk SAM output. Running
the aligner once on the unchunked input selected one additional pair in this
test. Chunking can therefore change alignment results; this experiment does
not establish biological accuracy or generalize the size of the difference.

## Restriction-fragment filtering

Restriction-fragment filtering is **off by default**, including in the example
configuration. This matches runHiC's default of leaving `--add-frag` off; it
does not imply that the rest of the pipelines have identical defaults or results.
When enabled, the workflow assigns each mapped end to a predicted restriction
fragment and removes pairs whose ends fall in the same fragment of the same
chromosome. This includes dangling and self-ligation pairs; other orientations
within a fragment are also removed. The remaining pairs are deduplicated and
used to generate the matrix. Enabling this filter therefore changes the output.
`filter-stats.json` records removed and retained counts.

Choose a fragment definition that matches the actual experiment. For the
two-enzyme Arima Hi-C scheme, enable in-silico digestion with:

```yaml
digest:
  enzyme: Arima
  reference: hg38.fa
```

The program scans the reference for MboI (`^GATC`) and HinfI (`G^ANTC`)
sites and generates a BED using runHiC's coordinate convention. `^` marks
the cleavage position and `N` is a variable base. The uncompressed FASTA
must contain every chromosome in `chromsizes` at the stated length. The BED is
written to the scratch workspace and removed after a successful run. On the hg38 reference used for local testing, the
generated BED matched runHiC's byte for byte. This does not establish
equivalence for every possible reference or Arima kit chemistry.

Alternatively, provide a precomputed fragment BED:

```yaml
fragments: hg38.fragments.bed
```

Use a BED from the correct reference assembly and enzyme or enzyme combination.
It must contain three columns (chromosome, start, end), with contiguous
intervals starting at zero for each chromosome. Automatic digestion currently
supports only the two-enzyme Arima scheme; use a precomputed BED for other
restriction enzymes. Specify at most one of `digest` and `fragments`.

For experiments without fixed restriction sites, such as Micro-C, leave both
settings unset; do not use an unrelated restriction-fragment BED. The rest of
the pipeline still runs, but Micro-C and other experiment-specific workflows
have not yet been validated on real datasets.

## Alternative inputs and resources

To start from pairs instead,
replace `fastq` with `"input_pairs": ["a.pairs", "b.pairs"]` and set
`"input_sorted": true` only for already sorted inputs. The output directory
must not exist; the scratch directory must exist and have sufficient space.
Scratch intermediates are removed after success; on failure their path is
recorded in `work-dir.txt` for inspection. `status.json` is written only
on success; logs and partial output may
remain on failure. The default ten resolutions and 1,000 ICE iterations can
be changed with `resolutions` (integer multiples of the first value) and
`balance_max_iterations`. A converged balance is not guaranteed for sparse
data; inspect each resolution's balancing attributes.

`merge_scratch_dir` is optional and defaults to `scratch_dir`. Its temporary
space needs can approach the size of the uncompressed merged pairs. `/dev/shm`
uses memory and may count against a cluster job's memory limit.

To confine merge, filtering, deduplication, statistics, matrix creation and
optional pair compression to specified Linux CPU IDs, set e.g.
`"postprocess_cpus": [0, 2, 4, 6, 8, 10]` inside `resources`. The IDs must
be within the process's current affinity, including any Slurm cpuset; an empty,
duplicated or unavailable list is rejected. Indexing, alignment and chunk
sorting keep their own resource settings. CPU IDs are logical CPUs; choose
physical-core siblings deliberately if that matters for a benchmark.

The base-Cooler converter uses a bounded reader thread to supply batches while
`matrix_threads` Rayon workers parse and bin records in parallel. The same
workers sort the resulting runs. The reader is an additional thread; explicit
`postprocess_cpus` confines it to the chosen CPU set as well.

## FASTQ indexing and reuse

Indexing scans the compressed input to build gzip recovery points but validates
only boundary read IDs; it does not validate all interior pairs.
Boundary search assumes correctly paired inputs with unique read IDs and uses
the first matching R2 ID; it does not detect repeated IDs in the search window.
Index construction uses rapidgzip; the two mates are indexed concurrently.
The former experimental `fastq.index_backend` setting has been removed; delete
it from older YAML configurations.
The generated indexes retain checkpoints at approximately 64 MiB decoded intervals. During mapping,
all aligners read and validate mate records in pairs, including their read IDs.
For external BWA, one worker thread per chunk reads both mates with rapidgzip,
interleaves the records, and streams them to BWA in paired-input mode. A
separate thread converts BWA's SAM output to pairs. Built-in BWA-MEM3 reads
paired records and converts its structured alignment fields directly to pairs.
For Chromap, one worker feeds both mate pipes with bounded
memory, and Chromap writes pairs directly. If Chromap consumes one mate substantially
ahead of the other, the lagging input is temporarily buffered on disk in the
mapping workspace. Allow enough free disk for up to one uncompressed mate per
concurrent chunk. Older
sparse indexes with gaps over 128 MiB are rejected;
rebuild them before reuse to avoid excessive memory use. Actual checkpoint
gaps depend on the gzip block boundaries. `jobs × threads_per_job`
sets aligner threads alone; FASTQ streaming and SAM conversion also consume CPU,
so this is not a global CPU limit. Sorting jobs, merge workers and
matrix threads are separate limits rather than a shared process-tree limit.

To build an index independently and reuse it on later runs:

```sh
./rust-hic index \
  data/R1.fastq.gz data/R2.fastq.gz -o new-index-directory \
  --chunk-size 2GiB
# Add index_dir: new-index-directory under fastq before run.
./rust-hic run config.yaml
```

`index` takes its settings from command-line options; it does not read a
configuration file and needs no aligner or reference. Paths are relative to the
current working directory. The output directory may already exist, but index
output files must not already exist. Run
`rust-hic index --help` for all options. Indexing progress is shown on the
terminal and saved to `index.log` in the output directory. During rapidgzip
scanning, progress reports an estimated percentage, decoded bytes, elapsed time,
and average speed for each mate. The percentage uses the compression ratio of
an initial sample and may be inaccurate if the ratio changes; it can be unknown
if the sample cannot be used. By
default, chunk count is calculated from the larger compressed FASTQ file and a
2 GiB target. Boundaries are chosen by uncompressed offsets and aligned on read IDs,
so actual compressed chunk sizes may vary.
Use `--chunk-size` with bytes or units such as `GB` and `GiB`, or use `--chunks`
to specify a count instead. The options are mutually exclusive. The main pipeline
uses the same 2 GiB default, configurable with `fastq.chunk_size_bytes`; if
`fastq.chunks` is specified, it takes precedence. When reusing an index, omit
`fastq.chunks` to use its recorded count, or set it to that count. The run checks the index's input mate paths,
chunk count, input fingerprints and boundary IDs. Fingerprints cover file
metadata and file edges, not the complete compressed content. If
`expected_pairs` is specified, the total is checked after alignment; a newly
built index does not provide a full pair count. Reusing an index avoids
rebuilding recovery points; interior paired IDs are validated when each chunk is mapped.

## Errors and retained files

`run --check` reports the failing setting or input path and exits nonzero.
Fix the reported problem and run it again; it does not reserve space or prevent
files from changing after the check. A regular `run` performs the same
preflight automatically. If the pipeline fails after creating its output
directory, inspect `failure.json` if present, `logs/events.jsonl`, stage logs, and the scratch
path in `work-dir.txt`. Scratch is retained on failure.
Choose a new output directory for a retry; existing outputs are not overwritten.

Standalone `map`, `clean`, and `matrix` retain their temporary workspace and
print its path when processing fails. An existing output file is never replaced.
The program reports invalid options and missing files before starting work,
but it cannot predict all later filesystem, aligner, data, or memory failures.
