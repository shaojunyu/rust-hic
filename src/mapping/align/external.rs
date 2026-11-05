//! External aligner processes and FIFO transport.
use super::*;

fn fifo(path: &Path) -> Result<()> {
    let name = CString::new(path.as_os_str().as_bytes())?;
    if unsafe { libc::mkfifo(name.as_ptr(), 0o600) } != 0 {
        return Err(std::io::Error::last_os_error()).context("create FIFO");
    };
    Ok(())
}
fn open_writer(path: &Path, cancel: &AtomicBool) -> Result<File> {
    loop {
        if cancel.load(Ordering::Relaxed) {
            bail!("alignment stopped before opening FASTQ FIFO")
        };
        match OpenOptions::new()
            .write(true)
            .custom_flags(libc::O_NONBLOCK)
            .open(path)
        {
            Ok(f) => {
                let flags =
                    unsafe { libc::fcntl(std::os::fd::AsRawFd::as_raw_fd(&f), libc::F_GETFL) };
                if flags < 0
                    || unsafe {
                        libc::fcntl(
                            std::os::fd::AsRawFd::as_raw_fd(&f),
                            libc::F_SETFL,
                            flags & !libc::O_NONBLOCK,
                        )
                    } < 0
                {
                    return Err(std::io::Error::last_os_error().into());
                }
                return Ok(f);
            }
            Err(e) if e.raw_os_error() == Some(libc::ENXIO) => {
                thread::sleep(Duration::from_millis(25))
            }
            Err(e) => return Err(e.into()),
        }
    }
}
fn terminate(p: &mut Child) {
    if p.try_wait().ok().flatten().is_none() {
        let _ = p.kill();
    }
    let _ = p.wait();
}

// Chromap emits either pairs or SAM, so retaining its alignments requires a
// second pass. Keeping the first pass unchanged preserves its Hi-C semantics.
fn chromap_bam(
    c: &FastqConfig,
    index: &Path,
    dir: &Path,
    logs: &Path,
    i: usize,
    cancelled: &AtomicBool,
) -> Result<()> {
    let f1 = dir.join("bam-read1.fifo");
    let f2 = dir.join("bam-read2.fifo");
    let sam_output = dir.join("bam-output.sam");
    fifo(&f1)?;
    fifo(&f2)?;
    // A chunk-specific stdout symlink also gives Chromap's low-memory
    // temporary files unique names without changing its working directory.
    std::os::unix::fs::symlink("/dev/stdout", &sam_output)?;
    let mut child = Command::new(&c.aligner_executable)
        .args(["--preset", "hic", "--low-mem", "--SAM", "-x"])
        .arg(&c.aligner_index)
        .arg("-r")
        .arg(&c.reference)
        .arg("-t")
        .arg(c.threads_per_job.to_string())
        .arg("-1")
        .arg(&f1)
        .arg("-2")
        .arg(&f2)
        .arg("-o")
        .arg(&sam_output)
        .args(&c.aligner_args)
        .stdout(Stdio::piped())
        .stderr(Stdio::from(File::create(
            logs.join(format!("align-{i:04}-bam.log")),
        )?))
        .spawn()
        .context("start Chromap BAM pass")?;
    let stdout = child.stdout.take().context("capture Chromap SAM")?;
    let cancel = AtomicBool::new(false);
    let result = thread::scope(|scope| -> Result<()> {
        let (tx, rx) = mpsc::channel();
        let worker = tx.clone();
        let cancel_ref = &cancel;
        let first = &f1;
        let second = &f2;
        scope.spawn(move || {
            let result =
                crate::fastq_index::stream_split(index, i, [first, second], cancel_ref).map(|_| ());
            let _ = worker.send(result);
        });
        let bam_path = dir.join("alignments.bam");
        scope.spawn(move || {
            let result = super::super::bam::sam_to_bam(BufReader::new(stdout), &bam_path);
            let _ = tx.send(result);
        });
        let mut finished = 0;
        let mut status = None;
        loop {
            if cancelled.load(Ordering::Acquire) {
                cancel.store(true, Ordering::Release);
                terminate(&mut child);
                bail!("alignment cancelled")
            }
            while let Ok(result) = rx.try_recv() {
                if let Err(e) = result {
                    cancel.store(true, Ordering::Release);
                    terminate(&mut child);
                    return Err(e).context("Chromap BAM pass");
                }
                finished += 1;
            }
            if status.is_none() {
                status = child.try_wait()?;
            }
            if let Some(s) = status {
                if !s.success() {
                    cancel.store(true, Ordering::Release);
                    bail!("Chromap BAM pass failed: {s}")
                }
                if finished == 2 {
                    break;
                }
            }
            thread::sleep(Duration::from_millis(50));
        }
        Ok(())
    });
    if result.is_err() {
        cancel.store(true, Ordering::Release);
        terminate(&mut child);
    }
    result?;
    fs::remove_file(f1)?;
    fs::remove_file(f2)?;
    fs::remove_file(sam_output)?;
    Ok(())
}
#[allow(clippy::too_many_arguments)] // Chunk context and shared aligner state.
pub(super) fn align_one(
    c: &FastqConfig,
    index: &Path,
    work: &Path,
    logs: &Path,
    i: usize,
    mem3: Option<(&Arc<bwa_mem3_rs::BwaIndex>, &Arc<bwa_mem3_rs::MemOpts>)>,
    cancelled: &AtomicBool,
    keep_bam: bool,
) -> Result<(PathBuf, u64)> {
    if cancelled.load(Ordering::Acquire) {
        bail!("alignment cancelled")
    }
    if let Some((reference, opts)) = mem3 {
        let output = work.join(format!("align-{i:04}")).join("alignments.pairs");
        return super::super::bwa_mem3::align_one(
            c, index, work, logs, i, output, reference, opts, cancelled, keep_bam,
        );
    }
    let started = Instant::now();
    let dir = work.join(format!("align-{i:04}"));
    fs::create_dir(&dir)?;
    let f1 = dir.join("read1.fifo");
    let f2 = dir.join("read2.fifo");
    let is_chromap = c.aligner == "chromap";
    fifo(&f1)?;
    if is_chromap {
        fifo(&f2)?;
    }
    let output = dir.join("alignments.pairs");
    eprintln!("event=chunk_start chunk={i} output={}", output.display());
    let logpath = logs.join(format!("align-{i:04}.log"));
    let log = File::create(&logpath)?;
    let mut command = Command::new(&c.aligner_executable);
    if is_chromap {
        command
            .args(["--preset", "hic", "--low-mem", "-x"])
            .arg(&c.aligner_index)
            .arg("-r")
            .arg(&c.reference)
            .arg("-t")
            .arg(c.threads_per_job.to_string())
            .arg("-1")
            .arg(&f1)
            .arg("-2")
            .arg(&f2)
            .arg("-o")
            .arg(&output)
            .args(&c.aligner_args)
            .stdout(Stdio::null());
    } else {
        command
            .args(["mem", "-SP", "-p", "-t"])
            .arg(c.threads_per_job.to_string())
            .args(&c.aligner_args)
            .arg(&c.aligner_index)
            .arg(&f1)
            .stdout(Stdio::piped());
    }
    let mut aligner = command
        .stderr(Stdio::from(log))
        .spawn()
        .context("start aligner")?;
    enum Finished {
        Stream(Result<()>),
        Parser(Result<super::super::sam::Counts>),
        Bam(Result<()>),
    }
    let cancel = Arc::new(AtomicBool::new(false));
    let (tx, rx) = mpsc::channel();
    let count = thread::scope(|scope| -> Result<Option<super::super::sam::Counts>> {
        let worker_cancel = Arc::clone(&cancel);
        let worker_tx = tx.clone();
        let first = f1.clone();
        let second = f2.clone();
        let errpath = logs.join(format!("align-{i:04}-read.log"));
        let worker = thread::Builder::new()
            .name(format!("fq-{i:04}"))
            .spawn_scoped(scope, move || {
                let result = (|| -> Result<()> {
                    let mut log = File::create(errpath)?;
                    if is_chromap {
                        let bytes = crate::fastq_index::stream_split(
                            index,
                            i,
                            [&first, &second],
                            &worker_cancel,
                        )?;
                        writeln!(
                            log,
                            "paired_stream chunk={i} R1_bytes={} R2_bytes={}",
                            bytes[0], bytes[1]
                        )?;
                    } else {
                        let output = open_writer(&first, &worker_cancel)?;
                        let pairs =
                            crate::fastq_index::stream_paired(index, i, output, &worker_cancel)?;
                        writeln!(log, "paired_stream chunk={i} pairs={pairs}")?;
                    }
                    Ok(())
                })();
                let _ = worker_tx.send(Finished::Stream(result));
            });
        if let Err(error) = worker {
            cancel.store(true, Ordering::Relaxed);
            terminate(&mut aligner);
            return Err(error).context("start paired FASTQ worker");
        }
        if let Some(stdout) = aligner.stdout.take() {
            let tx = tx.clone();
            let outpath = output.clone();
            let reader: Box<dyn std::io::Read + Send> = if keep_bam {
                let (read, write) = std::os::unix::net::UnixStream::pair()?;
                let bam_tx = tx.clone();
                let bam_path = dir.join("alignments.bam");
                scope.spawn(move || {
                    let result = super::super::bam::sam_to_bam(BufReader::new(read), &bam_path);
                    let _ = bam_tx.send(Finished::Bam(result));
                });
                Box::new(super::super::bam::Tee {
                    input: stdout,
                    copy: write,
                })
            } else {
                Box::new(stdout)
            };
            let parser = thread::Builder::new()
                .name(format!("pairs-{i:04}"))
                .spawn_scoped(scope, move || {
                    let result = (|| -> Result<_> {
                        let output = BufWriter::with_capacity(1024 * 1024, File::create(outpath)?);
                        super::super::sam::convert(
                            BufReader::with_capacity(1024 * 1024, reader),
                            output,
                            c.sam.min_mapq,
                            c.sam.max_molecule_size,
                            c.sam.max_inter_align_gap,
                        )
                    })();
                    let _ = tx.send(Finished::Parser(result));
                });
            if let Err(error) = parser {
                cancel.store(true, Ordering::Relaxed);
                terminate(&mut aligner);
                return Err(error).context("start SAM-to-pairs worker");
            }
        }
        drop(tx);
        let (mut done, mut parsed, mut result) = (0, None, None);
        let mut bam_done = !keep_bam || is_chromap;
        loop {
            if cancelled.load(Ordering::Acquire) {
                result = Some(anyhow::anyhow!("alignment cancelled"));
                break;
            }
            while let Ok(r) = rx.try_recv() {
                match r {
                    Finished::Stream(r) => {
                        done += 1;
                        if let Err(e) = r {
                            result = Some(e);
                            break;
                        }
                    }
                    Finished::Parser(r) => match r {
                        Ok(n) => parsed = Some(n),
                        Err(e) => {
                            result = Some(e);
                            break;
                        }
                    },
                    Finished::Bam(r) => {
                        bam_done = true;
                        if let Err(e) = r {
                            result = Some(e);
                            break;
                        }
                    }
                }
            }
            if result.is_some() {
                break;
            }
            let status = match aligner.try_wait() {
                Ok(s) => s,
                Err(e) => {
                    result = Some(e.into());
                    break;
                }
            };
            if let Some(s) = status {
                cancel.store(true, Ordering::Relaxed);
                if !s.success() {
                    result = Some(anyhow::anyhow!("{} chunk {i} failed: {s}", c.aligner));
                    break;
                }
                if done == 1 && (is_chromap || parsed.is_some()) && bam_done {
                    break;
                }
            }
            thread::sleep(Duration::from_millis(50));
        }
        if let Some(e) = result {
            cancel.store(true, Ordering::Relaxed);
            terminate(&mut aligner);
            return Err(e);
        };
        Ok(parsed)
    })
    .and_then(|parsed| -> Result<u64> {
        if let Some(n) = parsed {
            eprintln!(
                "event=sam_parsed chunk={i} input_pairs={} selected_pairs={}",
                n.input_pairs, n.selected_pairs
            );
            serde_json::to_writer(
                File::create(logs.join(format!("align-{i:04}-sam-stats.json")))?,
                &n,
            )?;
            return Ok(n.input_pairs);
        }
        let text = fs::read_to_string(&logpath)?;
        let count = text
            .lines()
            .filter_map(|line| line.trim().strip_prefix("Number of reads:"))
            .map(|s| s.trim().trim_end_matches('.').parse::<u64>())
            .collect::<std::result::Result<Vec<_>, _>>()?;
        if count.len() != 1 || count[0] % 2 != 0 {
            bail!("chromap chunk {i} did not report a valid read count")
        };
        Ok(count[0] / 2)
    })?;
    fs::remove_file(f1)?;
    if is_chromap {
        fs::remove_file(f2)?;
    }
    if keep_bam && is_chromap {
        chromap_bam(c, index, &dir, logs, i, cancelled)?;
    }
    serde_json::to_writer_pretty(
        File::create(logs.join(format!("align-{i:04}-stats.json")))?,
        &serde_json::json!({"chunk":i,"aligner":c.aligner,"input_pairs":count,"elapsed_seconds":started.elapsed().as_secs_f64()}),
    )?;
    eprintln!("event=chunk_done chunk={i} output={}", output.display());
    Ok((output, count))
}
