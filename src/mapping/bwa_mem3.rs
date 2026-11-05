//! In-process bwa-mem3 alignment and direct structured-field to pairs conversion.
use super::{align::FastqConfig, sam};
use crate::fastq_index::{FastqPair, PairedFastqReader};
use anyhow::{bail, Context, Result};
use bwa_mem3_rs::{
    AlignScratch, AlignedFields, AlignedFieldsSink, BwaIndex, IdBases, Mate, MemOpts, MemPeStat,
    ReadPair, RecordOrigin, RecordSink, ResidentCohort,
};
use crossbeam_channel::{bounded, SendTimeoutError};
use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc,
    },
    thread,
    time::{Duration, Instant},
};

const BATCH_PAIRS: usize = 4096;

pub fn options(c: &FastqConfig) -> Result<MemOpts> {
    let mut opts = MemOpts::new()?;
    opts.set_pe(true).set_hic_no_pairing_no_rescue();
    let mut args = c.aligner_args.iter();
    while let Some(flag) = args.next() {
        if !["-k", "-w", "-T", "-A", "-B", "-c"].contains(&flag.as_str()) {
            bail!(
                "unsupported built-in bwa-mem3 argument {flag}; supported: -k, -w, -T, -A, -B, -c"
            )
        }
        let value = args
            .next()
            .with_context(|| format!("missing value for bwa-mem3 option {flag}"))?;
        let v: i32 = value
            .parse()
            .with_context(|| format!("invalid value for bwa-mem3 option {flag}"))?;
        match flag.as_str(){
            "-k" if v>0 => {opts.set_min_seed_len(v);},
            "-w" if v>0 => {opts.set_band_width(v);},
            "-T" if v>=0 => {opts.set_minimum_score(v);},
            "-A" if v>0 => {opts.set_match_score(v);},
            "-B" if v>=0 => {opts.set_mismatch_penalty(v);},
            "-c" if v>0 => {opts.set_max_occurrences(v);},
            _ => bail!("unsupported or invalid built-in bwa-mem3 argument {flag}; supported: -k, -w, -T, -A, -B, -c"),
        }
    }
    Ok(opts)
}

struct Fields<'a> {
    index: &'a BwaIndex,
    min_mapq: u8,
    records: Vec<(Vec<sam::Alignment>, Vec<sam::Alignment>)>,
    error: Option<anyhow::Error>,
}
impl AlignedFieldsSink for Fields<'_> {
    fn emit(
        &mut self,
        origin: RecordOrigin,
        mate: Mate,
        _primary: bool,
        fields: AlignedFields<'_>,
    ) {
        if self.error.is_some() {
            return;
        }
        let RecordOrigin::Pair(i) = origin else {
            self.error = Some(anyhow::anyhow!("unexpected single-end alignment"));
            return;
        };
        let Some(pair) = self.records.get_mut(i) else {
            self.error = Some(anyhow::anyhow!("invalid alignment pair index"));
            return;
        };
        match sam::from_fields(&fields, self.index, self.min_mapq) {
            Ok(a) => match mate {
                Mate::R1 => pair.0.push(a),
                Mate::R2 => pair.1.push(a),
            },
            Err(e) => self.error = Some(e),
        }
    }
}
struct Batch {
    first: u64,
    pairs: Vec<FastqPair>,
}
struct Done {
    data: Vec<u8>,
    counts: sam::Counts,
    bam_data: Option<Vec<u8>>,
}

struct Packed<'a> {
    fields: Fields<'a>,
    data: Vec<u8>,
}
impl RecordSink for Packed<'_> {
    fn emit(&mut self, origin: RecordOrigin, body: &[u8]) {
        if self.fields.error.is_some() {
            return;
        }
        let result = (|| -> Result<()> {
            if body.len() < 32 {
                bail!("truncated BAM core from bwa-mem3")
            }
            let tid = i32::from_le_bytes(body[0..4].try_into()?);
            let pos = i32::from_le_bytes(body[4..8].try_into()?);
            let name_len = body[8] as usize;
            let mapq = body[9];
            let cigar_len = u16::from_le_bytes(body[12..14].try_into()?) as usize;
            let flag = u16::from_le_bytes(body[14..16].try_into()?);
            let start = 32usize
                .checked_add(name_len)
                .context("invalid BAM name length")?;
            let end = start
                .checked_add(cigar_len * 4)
                .context("invalid BAM cigar length")?;
            let cigar_bytes = body.get(start..end).context("truncated BAM cigar")?;
            let cigar: Vec<u32> = cigar_bytes
                .chunks_exact(4)
                .map(|x| u32::from_le_bytes(x.try_into().unwrap()))
                .collect();
            let fields = AlignedFields {
                tid,
                pos,
                next_tid: 0,
                next_pos: 0,
                tlen: 0,
                flag,
                bin: 0,
                mapq,
                query_start: 0,
                query_end: 0,
                cigar: &cigar,
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
            let mate = match flag & 0xc0 {
                0x40 => Mate::R1,
                0x80 => Mate::R2,
                _ => bail!("invalid BAM mate flag"),
            };
            self.fields.emit(origin, mate, false, fields);
            self.data
                .extend_from_slice(&u32::try_from(body.len())?.to_le_bytes());
            self.data.extend_from_slice(body);
            Ok(())
        })();
        if let Err(e) = result {
            self.fields.error = Some(e)
        }
    }
}
fn align_batch(
    index: &BwaIndex,
    opts: &MemOpts,
    scratch: &mut AlignScratch,
    batch: Batch,
    c: &FastqConfig,
    keep_bam: bool,
) -> Result<Done> {
    let cohort = ResidentCohort::new(false)?;
    let mut range = cohort.reserve_pairs(batch.pairs.len())?;
    let pairs: Vec<_> = batch
        .pairs
        .iter()
        .map(|p| ReadPair {
            name_r1: &p.name,
            seq_r1: &p.seq1,
            qual_r1: Some(&p.qual1),
            name_r2: &p.name,
            seq_r2: &p.seq2,
            qual_r2: Some(&p.qual2),
        })
        .collect();
    cohort.write_pairs(&mut range, &pairs)?;
    cohort.seed_extend(index, opts, scratch, &mut range)?;
    // Both pairing and mate rescue are disabled (-S -P), so no insert-size
    // model is consulted. The API still requires a model for paired reads.
    let pestat = MemPeStat::zero()?;
    let mut sink = Fields {
        index,
        min_mapq: c.sam.min_mapq,
        records: (0..pairs.len()).map(|_| (Vec::new(), Vec::new())).collect(),
        error: None,
    };
    let ids = IdBases {
        first_single_id: 0,
        first_pair_id: batch.first,
    };
    let bam_data = if keep_bam {
        let mut packed = Packed {
            fields: sink,
            data: Vec::new(),
        };
        cohort.pair_emit(
            index,
            opts,
            scratch,
            &mut range,
            Some(&pestat),
            ids,
            0,
            &mut packed,
        )?;
        sink = packed.fields;
        Some(packed.data)
    } else {
        cohort.pair_emit_fields(
            index,
            opts,
            scratch,
            &mut range,
            Some(&pestat),
            ids,
            0,
            &mut sink,
        )?;
        None
    };
    if let Some(e) = sink.error {
        return Err(e);
    }
    let mut out = Vec::new();
    let mut counts = sam::Counts::default();
    for (p, (mut a, mut b)) in batch.pairs.iter().zip(sink.records) {
        if a.is_empty() || b.is_empty() {
            bail!("bwa-mem3 emitted incomplete paired alignment")
        }
        let name = std::str::from_utf8(&p.name).context("FASTQ read name is not UTF-8")?;
        sam::write_pair(
            name,
            &mut a,
            &mut b,
            &mut out,
            &mut counts,
            c.sam.max_molecule_size,
            c.sam.max_inter_align_gap,
        )?;
    }
    Ok(Done {
        data: out,
        counts,
        bam_data,
    })
}

#[allow(clippy::too_many_arguments)] // Shared chunk context and preloaded aligner state.
pub fn align_one(
    c: &FastqConfig,
    index_dir: &Path,
    work: &Path,
    logs: &Path,
    i: usize,
    output: PathBuf,
    index: &Arc<BwaIndex>,
    opts: &Arc<MemOpts>,
    cancelled: &AtomicBool,
    keep_bam: bool,
) -> Result<(PathBuf, u64)> {
    let started = Instant::now();
    let dir = work.join(format!("align-{i:04}"));
    fs::create_dir(&dir)?;
    let mut bam_writer = if keep_bam {
        Some(super::bam::create_index_writer(
            &dir.join("alignments.bam"),
            index,
        )?)
    } else {
        None
    };
    let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(&output)?);
    writeln!(writer,"## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type")?;
    eprintln!("event=chunk_start chunk={i} output={}", output.display());
    let n = c.threads_per_job;
    let cap = n.saturating_mul(2).clamp(1, 64);
    let (tx, rx) = bounded::<Batch>(cap);
    let (done_tx, done_rx) = bounded::<Result<Done>>(cap);
    let cancel = AtomicBool::new(false);
    let (input, selected) = thread::scope(|scope| -> Result<(u64, u64)> {
        let decoder_cancel = &cancel;
        let decoder = thread::Builder::new()
            .name(format!("fq-{i:04}"))
            .spawn_scoped(scope, move || -> Result<u64> {
                let mut reader = PairedFastqReader::open(index_dir, i)?;
                let mut first = 0u64;
                loop {
                    if decoder_cancel.load(Ordering::Relaxed) || cancelled.load(Ordering::Acquire) {
                        bail!("alignment cancelled")
                    }
                    let mut pairs = Vec::with_capacity(BATCH_PAIRS);
                    while pairs.len() < BATCH_PAIRS {
                        match reader.next_pair()? {
                            Some(p) => pairs.push(p),
                            None => break,
                        }
                    }
                    if pairs.is_empty() {
                        break;
                    }
                    let count = pairs.len() as u64;
                    let mut batch = Batch { first, pairs };
                    loop {
                        match tx.send_timeout(batch, Duration::from_millis(100)) {
                            Ok(()) => break,
                            Err(SendTimeoutError::Timeout(b)) => {
                                batch = b;
                                if decoder_cancel.load(Ordering::Relaxed)
                                    || cancelled.load(Ordering::Acquire)
                                {
                                    bail!("alignment cancelled")
                                }
                            }
                            Err(SendTimeoutError::Disconnected(_)) => {
                                bail!("alignment workers stopped")
                            }
                        }
                    }
                    first = first
                        .checked_add(count)
                        .context("FASTQ pair count overflow")?;
                }
                Ok(first)
            })?;
        let mut workers = Vec::new();
        for j in 0..n {
            let rx = rx.clone();
            let send = done_tx.clone();
            let idx = Arc::clone(index);
            let opts = Arc::clone(opts);
            let cancel = &cancel;
            workers.push(
                thread::Builder::new()
                    .name(format!("bwa-{i:03}-{j:02}"))
                    .spawn_scoped(scope, move || -> Result<()> {
                        let mut scratch = AlignScratch::new()?;
                        while !cancel.load(Ordering::Relaxed) && !cancelled.load(Ordering::Acquire)
                        {
                            let batch = match rx.recv_timeout(Duration::from_millis(100)) {
                                Ok(b) => b,
                                Err(crossbeam_channel::RecvTimeoutError::Timeout) => continue,
                                Err(crossbeam_channel::RecvTimeoutError::Disconnected) => break,
                            };
                            let result = align_batch(&idx, &opts, &mut scratch, batch, c, keep_bam);
                            let failed = result.is_err();
                            if send.send(result).is_err() {
                                break;
                            }
                            if failed {
                                cancel.store(true, Ordering::Relaxed);
                                break;
                            }
                        }
                        Ok(())
                    })?,
            );
        }
        drop(rx);
        drop(done_tx);
        let (mut input, mut selected) = (0u64, 0u64);
        let mut failure = None;
        while let Ok(result) = done_rx.recv() {
            match result {
                Ok(done) => {
                    if let Some(data) = done.bam_data {
                        if let Err(e) = bam_writer
                            .as_mut()
                            .expect("BAM writer")
                            .get_mut()
                            .write_all(&data)
                        {
                            failure = Some(e.into());
                            cancel.store(true, Ordering::Relaxed);
                            break;
                        }
                    }
                    if let Err(e) = writer.write_all(&done.data) {
                        failure = Some(e.into());
                        cancel.store(true, Ordering::Relaxed);
                        break;
                    }
                    input += done.counts.input_pairs;
                    selected += done.counts.selected_pairs;
                }
                Err(e) => {
                    failure = Some(e);
                    cancel.store(true, Ordering::Relaxed);
                    break;
                }
            }
        }
        drop(done_rx);
        let decoded = decoder
            .join()
            .map_err(|_| anyhow::anyhow!("FASTQ decoder panicked"))?;
        for worker in workers {
            let res = worker
                .join()
                .map_err(|_| anyhow::anyhow!("bwa-mem3 worker panicked"))?;
            if failure.is_none() {
                if let Err(e) = res {
                    failure = Some(e)
                }
            }
        }
        if let Some(e) = failure {
            return Err(e);
        }
        if cancelled.load(Ordering::Acquire) {
            bail!("alignment cancelled")
        }
        let decoded = decoded?;
        if decoded != input {
            bail!("bwa-mem3 processed {input} pairs but decoded {decoded}")
        }
        Ok((input, selected))
    })?;
    writer.flush()?;
    if let Some(mut writer) = bam_writer {
        super::bam::finish(&mut writer).context("finish BAM chunk")?;
    }
    serde_json::to_writer_pretty(
        File::create(logs.join(format!("align-{i:04}-stats.json")))?,
        &serde_json::json!({"chunk":i,"aligner":"bwa-mem3","input_pairs":input,"selected_pairs":selected,"elapsed_seconds":started.elapsed().as_secs_f64()}),
    )?;
    eprintln!(
        "event=chunk_done chunk={i} input_pairs={input} selected_pairs={selected} output={}",
        output.display()
    );
    Ok((output, input))
}
