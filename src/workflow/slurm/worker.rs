//! Execute Slurm stages and publish complete chunk artifacts atomically.
use super::state::{
    clear_failure, hash_file, now_ms, read_json, record_failure, write_json, ChunkDone, IndexDone,
    Plan, SCHEMA,
};
use crate::{
    mapping::align,
    workflow::{
        chunks, config::Config, finish_stages, layout::Layout, process::affinity_mask, RunClock,
    },
};
use anyhow::{bail, Context, Result};
use std::{
    env,
    fs::{self, File, OpenOptions},
    io::{self, Write},
    path::{Path, PathBuf},
    sync::atomic::AtomicBool,
    time::Duration,
};

fn load_plan(work: &Path) -> Result<(Plan, Config)> {
    let work = fs::canonicalize(work)
        .with_context(|| format!("open shared workspace {}", work.display()))?;
    let plan: Plan = read_json(&work.join("plan.json"))?;
    if plan.schema != SCHEMA || plan.work != work || !(1..=256).contains(&plan.chunks) {
        bail!("invalid or incompatible Slurm workflow plan")
    }
    let config_path = plan.output.join("config.resolved.yaml");
    if hash_file(&config_path)? != plan.config_sha256 {
        bail!("resolved configuration changed after submission")
    }
    if hash_file(&env::current_exe()?)? != plan.executable_sha256 {
        bail!("workflow executable changed after submission")
    }
    let c: Config = serde_yaml::from_reader(File::open(&config_path)?)?;
    if c.output_dir != plan.output || c.fastq.as_ref().and_then(|f| f.chunks) != Some(plan.chunks) {
        bail!("workflow configuration does not match its plan")
    }
    Ok((plan, c))
}

fn index_path(plan: &Plan, c: &Config) -> Result<PathBuf> {
    let fastq = c.fastq.as_ref().context("missing FASTQ configuration")?;
    Ok(fastq
        .index_dir
        .clone()
        .unwrap_or_else(|| plan.work.join("fastq-index")))
}

fn index_done(plan: &Plan, c: &Config) -> Result<IndexDone> {
    let state: IndexDone = read_json(&plan.work.join("index-complete.json"))
        .context("index job did not publish a completion marker")?;
    if state.schema != SCHEMA
        || state.run_id != plan.run_id
        || state.chunks != plan.chunks
        || state.end_ms < state.start_ms
    {
        bail!("index completion marker does not match the workflow")
    }
    let path = index_path(plan, c)?;
    if hash_file(&path.join("manifest.json"))? != state.manifest_sha256 {
        bail!("FASTQ index manifest changed after the index job")
    }
    Ok(state)
}

pub(crate) fn task(phase: &str, work: &Path, chunk: Option<usize>) -> Result<()> {
    let (plan, mut c) = load_plan(work)?;
    let layout = Layout::from_config(&plan.output, &c)?;
    if matches!((phase, chunk), ("index", None)) {
        clear_failure(&plan.output, "index")?;
    }
    let result = match (phase, chunk) {
        ("index", None) => run_index(&plan, &layout, &mut c),
        ("map", Some(i)) => run_chunk(&plan, &c, i),
        ("finish", None) => run_finish(&plan, &layout, &mut c),
        _ => bail!("invalid internal Slurm task"),
    };
    if let Err(error) = &result {
        if phase != "map" {
            record_failure(&plan.output, phase, error);
        }
    }
    result
}

fn run_index(plan: &Plan, layout: &Layout, c: &mut Config) -> Result<()> {
    if plan.work.join("index-complete.json").is_file() {
        index_done(plan, c)?;
        return Ok(());
    }
    let start_ms = now_ms()?;
    let fastq = c.fastq.as_mut().context("missing FASTQ configuration")?;
    let index = if fastq.index_dir.is_some() {
        align::index(
            fastq,
            &env::current_exe()?,
            &plan.work.join("unused-index"),
            &layout.log("fastq-index"),
            false,
            true,
        )?
    } else {
        let final_index = plan.work.join("fastq-index");
        if final_index.exists() {
            // A completed rename may have survived a preemption before its marker.
            align::validate_reusable_index(fastq, &final_index)?;
            crate::fastq_index::verify(&final_index, 1)?;
            final_index
        } else {
            let attempt = tempfile::Builder::new()
                .prefix("index-")
                .tempdir_in(plan.work.join("attempts"))?;
            let candidate = attempt.path().join("fastq-index");
            let result = align::index(
                fastq,
                &env::current_exe()?,
                &candidate,
                &layout.log("fastq-index"),
                false,
                false,
            );
            if let Err(e) = result {
                let retained = attempt.keep();
                return Err(e)
                    .with_context(|| format!("index attempt retained in {}", retained.display()));
            }
            fs::rename(&candidate, &final_index).context("publish completed FASTQ index")?;
            final_index
        }
    };
    if fastq.chunks != Some(plan.chunks) {
        bail!("index chunk count differs from submitted array")
    }
    let state = IndexDone {
        schema: SCHEMA,
        run_id: plan.run_id.clone(),
        chunks: plan.chunks,
        start_ms,
        end_ms: now_ms()?,
        manifest_sha256: hash_file(&index.join("manifest.json"))?,
    };
    write_json(&plan.work.join("index-complete.json"), &state)?;
    eprintln!(
        "event=slurm_index_complete chunks={} index={}",
        plan.chunks,
        index.display()
    );
    Ok(())
}

fn chunk_dir(plan: &Plan, chunk: usize) -> PathBuf {
    plan.work.join("chunks").join(format!("{chunk:04}"))
}

fn validate_chunk(plan: &Plan, c: &Config, chunk: usize, index: &IndexDone) -> Result<ChunkDone> {
    let dir = chunk_dir(plan, chunk);
    let done: ChunkDone = read_json(&dir.join("complete.json"))
        .with_context(|| format!("alignment chunk {chunk} is incomplete"))?;
    if done.schema != SCHEMA
        || done.run_id != plan.run_id
        || done.chunk != chunk
        || done.manifest_sha256 != index.manifest_sha256
        || done.end_ms < done.start_ms
    {
        bail!("invalid completion marker for alignment chunk {chunk}")
    }
    let files = dir.join(format!("align-{chunk:04}"));
    if fs::metadata(files.join("alignments.pairs"))?.len() != done.pairs_bytes
        || done.pairs_bytes == 0
    {
        bail!("alignment chunk {chunk} pairs file is missing or changed")
    }
    match (c.keep_bam, done.bam_bytes) {
        (true, Some(size))
            if size > 0 && fs::metadata(files.join("alignments.bam"))?.len() == size => {}
        (false, None) => {}
        _ => bail!("alignment chunk {chunk} BAM file is missing or changed"),
    }
    let stats: serde_json::Value = read_json(&dir.join(format!("align-{chunk:04}-stats.json")))?;
    if stats.get("chunk").and_then(|v| v.as_u64()) != Some(chunk as u64)
        || stats.get("input_pairs").and_then(|v| v.as_u64()) != Some(done.input_pairs)
    {
        bail!("alignment chunk {chunk} statistics disagree with its completion marker")
    }
    Ok(done)
}

fn run_chunk(plan: &Plan, c: &Config, chunk: usize) -> Result<()> {
    if chunk >= plan.chunks {
        bail!("chunk {chunk} is out of range 0..{}", plan.chunks)
    }
    let index = index_done(plan, c)?;
    let final_dir = chunk_dir(plan, chunk);
    if final_dir.exists() {
        validate_chunk(plan, c, chunk, &index)?;
        eprintln!("event=slurm_chunk_reused chunk={chunk}");
        return Ok(());
    }
    let start_ms = now_ms()?;
    let attempt = tempfile::Builder::new()
        .prefix(&format!("chunk-{chunk:04}-"))
        .tempdir_in(plan.work.join("attempts"))?;
    let fastq = c.fastq.as_ref().context("missing FASTQ configuration")?;
    let cancelled = AtomicBool::new(false);
    let result = align::map_with_completed(
        fastq,
        &index_path(plan, c)?,
        attempt.path(),
        attempt.path(),
        Some(chunk),
        &cancelled,
        c.keep_bam,
        |_, _| Ok(()),
    );
    if let Err(e) = result {
        let retained = attempt.keep();
        return Err(e).with_context(|| {
            format!(
                "alignment chunk {chunk} attempt retained in {}",
                retained.display()
            )
        });
    }
    // Refuse to publish if the index was replaced while this chunk ran.
    index_done(plan, c)?;
    let files = attempt.path().join(format!("align-{chunk:04}"));
    let pairs_bytes = fs::metadata(files.join("alignments.pairs"))?.len();
    let bam_bytes = if c.keep_bam {
        Some(fs::metadata(files.join("alignments.bam"))?.len())
    } else {
        None
    };
    let stats: serde_json::Value =
        read_json(&attempt.path().join(format!("align-{chunk:04}-stats.json")))?;
    let input_pairs = stats
        .get("input_pairs")
        .and_then(|v| v.as_u64())
        .context("missing chunk input count")?;
    let done = ChunkDone {
        schema: SCHEMA,
        run_id: plan.run_id.clone(),
        chunk,
        start_ms,
        end_ms: now_ms()?,
        manifest_sha256: index.manifest_sha256.clone(),
        input_pairs,
        pairs_bytes,
        bam_bytes,
    };
    write_json(&attempt.path().join("complete.json"), &done)?;
    let attempt = attempt.keep();
    if let Err(e) = fs::rename(&attempt, &final_dir) {
        if final_dir.exists() {
            validate_chunk(plan, c, chunk, &index)?;
            eprintln!(
                "event=slurm_chunk_already_complete chunk={chunk} unused_attempt={}",
                attempt.display()
            );
            return Ok(());
        }
        return Err(e).with_context(|| {
            format!(
                "publish chunk {chunk}; attempt retained in {}",
                attempt.display()
            )
        });
    }
    eprintln!("event=slurm_chunk_complete chunk={chunk} pairs={input_pairs}");
    Ok(())
}

fn write_initial_events(
    plan: &Plan,
    layout: &Layout,
    index: &IndexDone,
    chunks: &[ChunkDone],
) -> Result<()> {
    let first = chunks
        .iter()
        .map(|c| c.start_ms)
        .min()
        .context("no alignment chunks")?;
    let last = chunks
        .iter()
        .map(|c| c.end_ms)
        .max()
        .context("no alignment chunks")?;
    // Nodes normally have synchronized clocks. Clamp small clock skew so stage
    // durations never become negative in the report.
    let index_start = index.start_ms.max(plan.submitted_ms);
    let index_end = index.end_ms.max(index_start);
    let map_start = first.max(index_end);
    let map_end = last.max(map_start);
    let path = layout.events();
    let mut f = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)
        .with_context(|| format!("create workflow event log {}", path.display()))?;
    for (event, time) in [
        ("pipeline_start", plan.submitted_ms),
        ("fastq_index_start", index_start),
        ("fastq_index_done", index_end),
        ("alignment_start", map_start),
        ("alignment_done", map_end),
    ] {
        let elapsed = time.saturating_sub(plan.submitted_ms) as f64 / 1000.0;
        writeln!(
            f,
            "{}",
            serde_json::json!({"event":event,"elapsed_seconds":elapsed})
        )?;
    }
    f.flush()?;
    Ok(())
}

fn copy_chunk_logs(plan: &Plan, layout: &Layout, chunk: usize) -> Result<()> {
    let dir = chunk_dir(plan, chunk);
    let prefix = format!("align-{chunk:04}");
    for entry in fs::read_dir(&dir)? {
        let entry = entry?;
        if !entry.file_type()?.is_file() {
            continue;
        }
        let name = entry.file_name();
        let Some(name_text) = name.to_str() else {
            continue;
        };
        if name_text.starts_with(&prefix)
            && (name_text.ends_with(".log") || name_text.ends_with(".json"))
        {
            fs::copy(entry.path(), layout.logs().join(name))?;
        }
    }
    Ok(())
}

fn run_finish(plan: &Plan, layout: &Layout, c: &mut Config) -> Result<()> {
    if layout.status().exists() {
        bail!("workflow output already has a completion status")
    }
    let index = index_done(plan, c)?;
    let mut done = Vec::with_capacity(plan.chunks);
    for chunk in 0..plan.chunks {
        done.push(validate_chunk(plan, c, chunk, &index)?);
    }
    let total = done
        .iter()
        .try_fold(0u64, |sum, d| sum.checked_add(d.input_pairs))
        .context("alignment pair count overflow")?;
    let fastq = c.fastq.as_ref().context("missing FASTQ configuration")?;
    let manifest: serde_json::Value = read_json(&index_path(plan, c)?.join("manifest.json"))?;
    for expected in [
        fastq.expected_pairs,
        manifest.get("pairs").and_then(|x| x.as_u64()),
    ]
    .into_iter()
    .flatten()
    {
        if total != expected {
            bail!("aligner read {total} pairs but expected {expected}")
        }
    }
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(plan.work.join("finish.started"))
        .context("final processing was already started; inspect the retained workspace")?;
    let postprocess_cpus = affinity_mask(c.resources.postprocess_cpus.as_deref())?;
    for chunk in 0..plan.chunks {
        copy_chunk_logs(plan, layout, chunk)?;
    }
    write_initial_events(plan, layout, &index, &done)?;
    let offset = now_ms()?.saturating_sub(plan.submitted_ms);
    let started = RunClock::with_offset(Duration::from_millis(offset));
    crate::workflow::digest_stage(c, layout, &plan.work, started)?;
    let mapped: Vec<_> = (0..plan.chunks)
        .map(|i| chunk_dir(plan, i).join(format!("align-{i:04}/alignments.pairs")))
        .collect();
    let sorted = chunks::prepare_mapped(
        c,
        &mapped,
        &plan.work,
        layout,
        started,
        postprocess_cpus.as_ref(),
    )?;
    if c.keep_bam {
        let parts: Vec<_> = (0..plan.chunks)
            .map(|i| chunk_dir(plan, i).join(format!("align-{i:04}/alignments.bam")))
            .collect();
        crate::mapping::bam::merge(&parts, &layout.bam())?;
        crate::workflow::log(&layout.logs(), "bam_merge_done", started)?;
    }
    finish_stages(
        c,
        &plan.work,
        layout,
        &sorted,
        c.fragments.is_some(),
        false,
        started,
        &postprocess_cpus,
    )?;
    // A previous scheduler attempt may have failed before being requeued.
    match fs::remove_file(layout.failure()) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::NotFound => {}
        Err(e) => return Err(e).context("remove obsolete failure record"),
    }
    crate::workflow::log(&layout.logs(), "pipeline_done", started)?;
    write_json(
        &layout.status(),
        &serde_json::json!({
            "status":"complete", "elapsed_seconds":started.elapsed().as_secs_f64()
        }),
    )?;
    // Publish the completed result before deleting scratch. If a node is lost
    // during cleanup, the matrix and its completion status remain usable.
    match fs::remove_dir_all(&plan.work) {
        Ok(()) => {
            if let Err(e) = fs::remove_file(layout.work_dir_file()) {
                eprintln!("warning: unable to remove work-dir.txt: {e}");
            }
        }
        Err(e) => eprintln!(
            "warning: completed scratch workspace retained at {}: {e}",
            plan.work.display()
        ),
    }
    eprintln!("report: {}", layout.report_html().display());
    if let Err(e) = layout.write_manifest() {
        eprintln!(
            "warning: unable to write {}: {e:#}",
            layout.manifest().display()
        );
    }
    Ok(())
}
