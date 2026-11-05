# Slurm execution

`rust-hic run --mode slurm` runs the existing FASTQ-to-mcool workflow across Slurm
CPU nodes. It submits an index job, a job array with one alignment task per
chunk, and a final job that sorts, filters, merges and builds the matrix on one
node. Submission returns after the jobs are accepted; the submitting shell
need not stay open. No Nextflow installation or scheduler daemon is required.

The Slurm commands `sbatch`, `scontrol`, `squeue`, and `scancel` must be on `PATH` on the
submission host. Slurm must support job arrays, dependencies, and
`--kill-on-invalid-dep`; the submitting user must be allowed to update and
release their own array tasks. Each node must be able to access the input FASTQs,
reference, reference index, `scratch_dir` and output directory through the
same paths. An external aligner executable must also be accessible at the same
path on each node. The shared scratch filesystem must allow execution of the
copied rust-hic executable; node-local `/tmp` is not suitable for `scratch_dir`.
Keep inputs and reference indexes unchanged until the workflow has finished.

## Submit

Start from the [minimal](example.minimal.config.yaml) or
[complete](example.config.yaml) configuration. For example, to request 16
chunks, up to 16 concurrent alignment tasks, and 16 alignment threads per task:

```yaml
fastq:
  read1: /shared/sample_R1.fastq.gz
  read2: /shared/sample_R2.fastq.gz
  reference: /shared/hg38.fa
  aligner_index: /shared/hg38.fa
  chunks: 16
  jobs: 16
  threads_per_job: 16
  index_threads_per_mate: 8
```

Fill in the remaining required settings, including a new `output_dir` and an
existing shared `scratch_dir`, then submit:

```sh
rust-hic run config.yaml --mode slurm --account ACCOUNT --partition PARTITION \
  --time 04:00:00 --index-mem 16G --map-mem 32G --finish-mem 64G
```

Account and partition can be omitted if your Slurm site supplies defaults.
The time and memory values are **per job**, not totals for the workflow.
`--index-time`, `--map-time` and `--finish-time` can override `--time` for
individual stages. `--max-parallel` overrides `fastq.jobs` for the array.
The same `--keep-bam` and `--keep-pairs` options as local mode are available.
Add `--check` to validate the configuration and supported submission options
without submitting jobs or creating the output directory. It checks that
`sbatch` and `scontrol` are available, but does not ask the scheduler to
validate the account, partition, time or memory requests. It cannot verify
that paths are shared with compute nodes or that `postprocess_cpus` is allowed
on the final node; a reused FASTQ index is fully verified by the index job.
Slurm validates scheduler settings when jobs are submitted.

By default the index job requests `2 × index_threads_per_mate + 1` CPUs, each
alignment task requests `threads_per_job + 1` CPUs for the aligner and reader,
and the final job requests the maximum of `(sort_threads + 1) × sort_jobs`,
`merge_workers + dedup_workers + 3`, and `matrix_threads + 1`. These are
conservative requests, and actual CPU use varies by phase.
`--index-cpus`, `--map-cpus` and `--finish-cpus` override those requests;
they do not change the workers in the YAML configuration. Adjust the workers
as well when changing CPU requests. Slurm CPU units and actual placement depend
on the site. Check that the requests fit the partition.

Each BWA-MEM3 alignment task loads its own reference index. Memory use can
therefore grow with the number of tasks on a node; `--map-mem` is reserved for
**each** task. Shared-storage bandwidth, scheduler wait and the slowest chunk
also affect total time. The array permits execution on multiple nodes but does
not force a separate node for every chunk. More tasks or threads do not
necessarily make the workflow faster.

Site-specific options can be passed to all three jobs, for example
`--sbatch-arg=--qos=normal`. Repeat this option as needed, using the
`--name=value` form. Resource, array, dependency, log and other options managed
by rust-hic cannot be overridden this way. `--exclusive` and `--no-requeue`
are also accepted as individual `--sbatch-arg` values.

## Monitor, cancel and troubleshoot

```sh
rust-hic run status /shared/new-result
rust-hic run cancel /shared/new-result
```

Submission prints the job IDs and saves them in `slurm-jobs.json` inside the
result directory, along with an eight-character random group ID. The index and
final processing jobs are named `rust-hic-<group-id>-main`; each alignment task
is named `rust-hic-<group-id>-mapping-<chunk>`, with chunk numbers starting at
zero. `status` shows the group ID and workflow state, and queries `squeue`;
while the shared workspace exists, it also shows completed chunks. Slurm's
`sacct` can provide finished job states and resource usage. Scheduler
logs are in the result's `slurm/` directory. On failure, inspect `failure.json`
if present, those logs, and the shared workspace named in `work-dir.txt`. `cancel` only
cancels the recorded active jobs; it leaves output and scratch files for
inspection.

If `squeue` truncates job names, use
`squeue -u "$USER" -o "%.18i %.40j %.2t %M %R"` to show the full group and
chunk IDs.

The array starts only after the index job succeeds. The final job runs after
all array tasks have ended and refuses to create a matrix unless every chunk
has a valid completion marker and its expected files. Each task writes into a
separate attempt directory and publishes its completed chunk atomically; a
restarted index or alignment task can reuse a completed result. The final job
is submitted with `--no-requeue`: automatic recovery of partially completed
sorting or matrix work is not yet implemented. A failed run is retained for
inspection and does not overwrite an existing output directory.

Unlike local `run`, sorting starts after all distributed alignment tasks have
finished. The report's alignment interval runs from the first task's start to
the last task's finish and can include waiting for later array tasks. Total
elapsed time includes scheduler queue waits; consult Slurm accounting and
per-task logs when comparing compute time. The report lists final chunk
preparation separately; its alignment-through-preparation interval also
includes any wait for the final job and the optional restriction digest.
Successful workflows remove their
owned scratch workspace after publishing the result; a cleanup failure leaves
it in place and is reported in the final job log.
