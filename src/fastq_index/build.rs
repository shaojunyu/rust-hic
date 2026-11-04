//! Build and compact FASTQ recovery indexes.
use super::*;

pub(super) fn compact_existing(out: &Path, threads: u32) -> Result<()> {
    verify(out, threads)?;
    let m: Manifest = serde_json::from_reader(File::open(out.join("manifest.json"))?)?;
    compact_index(
        &out.join(&m.inputs[0].index),
        &m.boundaries
            .iter()
            .map(|b| b.r1_uncompressed_offset)
            .collect::<Vec<_>>(),
    )?;
    compact_index(
        &out.join(&m.inputs[1].index),
        &m.boundaries
            .iter()
            .map(|b| b.r2_uncompressed_offset)
            .collect::<Vec<_>>(),
    )?;
    verify(out, threads)
}
fn export_sparse(
    reader: &mut rapidgzip::RapidgzipReader,
    path: &Path,
    offsets: &[u64],
) -> Result<()> {
    let mut writer = sparse::SparseWriter::new(path, offsets)?;
    reader.export_to(&mut writer)?;
    writer.finish()?;
    Ok(())
}
/// Estimate compression ratio from the difference between two prefixes.
/// Both prefixes are truncated gzip streams, so the last read may fail; bytes
/// returned before that failure remain usable for this approximate estimate.
fn estimate_decoded_size(input: &Path) -> Option<u64> {
    const SAMPLE_STEP: u64 = 8 * 1024 * 1024;
    const MAX_DECODED: u64 = 256 * 1024 * 1024;
    fn sample(input: &Path, bytes: u64) -> Option<u64> {
        let mut source = File::open(input).ok()?.take(bytes);
        let mut prefix = tempfile::NamedTempFile::new().ok()?;
        if io::copy(&mut source, &mut prefix).ok()? != bytes {
            return None;
        }
        let mut reader = rapidgzip::RapidgzipReader::open(prefix.path(), 1).ok()?;
        let mut decoded = 0u64;
        let mut buffer = vec![0u8; 1024 * 1024];
        while let Ok(n) = reader.read(&mut buffer) {
            decoded = decoded.checked_add(n as u64)?;
            if n == 0 {
                break;
            }
            if decoded >= MAX_DECODED {
                return None;
            }
        }
        Some(decoded)
    }
    let compressed_size = input.metadata().ok()?.len();
    if compressed_size <= 2 * SAMPLE_STEP {
        return None;
    }
    let first = sample(input, SAMPLE_STEP)?;
    let second = sample(input, 2 * SAMPLE_STEP)?;
    let decoded_step = second.checked_sub(first).filter(|&n| n > 0)?;
    u64::try_from((decoded_step as u128 * compressed_size as u128) / SAMPLE_STEP as u128).ok()
}
/// Collect gzip recovery points with rapidgzip without parsing FASTQ records.
fn build_dense(
    input: &Path,
    threads: u32,
    mate: &str,
) -> Result<(rapidgzip::RapidgzipReader, u64)> {
    let estimate = estimate_decoded_size(input);
    let mut r = rapidgzip::RapidgzipReader::open(input, threads)?;
    let started = Instant::now();
    let mut last_report = started;
    let mut n = 0u64;
    let mut buffer = vec![0u8; 1024 * 1024];
    loop {
        let count = r.read(&mut buffer)?;
        n = n
            .checked_add(count as u64)
            .context("uncompressed size overflow")?;
        let elapsed = started.elapsed();
        if count == 0 || last_report.elapsed().as_secs() >= 3 {
            let seconds = elapsed.as_secs_f64().max(0.001);
            let progress = if count == 0 {
                "100.0%".to_owned()
            } else if let Some(total) = estimate.filter(|&total| total > 0) {
                format!("{:.1}%", (n as f64 * 100.0 / total as f64).min(99.0))
            } else {
                "unknown".to_owned()
            };
            eprintln!("index_progress mate={mate} estimated={progress} decoded={:.2} MiB elapsed={seconds:.1}s rate={:.1} MiB/s{}", n as f64 / 1048576.0, n as f64 / 1048576.0 / seconds, if count == 0 { " scan_complete" } else { "" });
            last_report = Instant::now();
        }
        if count == 0 {
            break;
        }
    }
    Ok((r, n))
}
/// Find a four-line FASTQ record near an arbitrary uncompressed offset.
/// Inputs are assumed to contain valid records with unique read IDs.
fn find_record(
    input: &mut (impl Read + Seek + ?Sized),
    start: u64,
    end: u64,
    wanted: Option<&[u8]>,
) -> Result<Option<(u64, Vec<u8>)>> {
    use std::collections::VecDeque;
    // Begin one byte before the requested range so that a record header
    // exactly at `start` is not discarded when skipping a partial line.
    let seek = start.saturating_sub(1);
    input.seek(SeekFrom::Start(seek))?;
    let mut r = BufReader::with_capacity(1024 * 1024, input);
    let mut pos = seek;
    let mut line = Vec::new();
    if start > 0 {
        let n = r.read_until(b'\n', &mut line)?;
        pos += n as u64;
    }
    let mut pending: VecDeque<(u64, Vec<u8>)> = VecDeque::new();
    loop {
        while pending.len() < 4 {
            let mut line = Vec::new();
            let off = pos;
            let n = r.read_until(b'\n', &mut line)?;
            if n == 0 {
                break;
            };
            pos += n as u64;
            pending.push_back((off, line));
        }
        if pending.len() < 4 || pending[0].0 >= end {
            break;
        }
        let h = &pending[0].1;
        let valid = h.starts_with(b"@") && h.len() > 1 && pending[2].1.starts_with(b"+");
        if valid {
            let id = h.split(|b| b.is_ascii_whitespace()).next().unwrap();
            let id = id
                .strip_suffix(b"/1")
                .or_else(|| id.strip_suffix(b"/2"))
                .unwrap_or(id);
            let off = pending[0].0;
            if id.len() > 1 && wanted.is_none_or(|w| w == id) {
                return Ok(Some((off, id.to_vec())));
            }
            for _ in 0..4 {
                pending.pop_front();
            }
        } else {
            pending.pop_front();
        }
    }
    Ok(None)
}

/// Search outward from the predicted offset, retaining the entire original
/// search window. Each probe reads forwards because gzip cannot be decoded
/// backwards; expanding the probes limits repeated seeks when the estimate is
/// inaccurate.
fn find_record_near(
    input: &mut (impl Read + Seek + ?Sized),
    start: u64,
    end: u64,
    predicted: u64,
    wanted: &[u8],
    initial_span: u64,
) -> Result<Option<(u64, Vec<u8>)>> {
    if start > end || initial_span == 0 {
        bail!("invalid boundary search range")
    }
    let center = predicted.clamp(start, end);
    let (mut left, mut right, mut span) = (center, center, initial_span);
    while left > start || right < end {
        let next_right = right.saturating_add(span).min(end);
        if next_right > right {
            if let Some(result) = find_record(input, right, next_right, Some(wanted))? {
                return Ok(Some(result));
            }
            right = next_right;
        }
        let next_left = left.saturating_sub(span).max(start);
        if next_left < left {
            if let Some(result) = find_record(input, next_left, left, Some(wanted))? {
                return Ok(Some(result));
            }
            left = next_left;
        }
        span = span.saturating_mul(2);
    }
    Ok(None)
}
pub(super) fn build_index(
    r1: &Path,
    r2: &Path,
    out: &Path,
    chunk_spec: &str,
    threads: u32,
    window: u64,
) -> Result<()> {
    if threads < 1 || window < 1 {
        bail!("threads and search window must be positive")
    }
    fs::create_dir_all(out)?;
    if out.join("manifest.json").exists() {
        bail!("manifest already exists; choose a new output directory")
    }
    let start = Instant::now();
    let before1 = fingerprint(r1, "r1.gzidx".into(), 0)?;
    let before2 = fingerprint(r2, "r2.gzidx".into(), 0)?;
    let chunks = if let Some(size) = chunk_spec.strip_prefix("size:") {
        let size: u64 = size.parse().context("invalid chunk size")?;
        if size == 0 {
            bail!("chunk size must be positive")
        }
        usize::try_from(before1.size.max(before2.size).div_ceil(size)).context("too many chunks")?
    } else {
        chunk_spec.parse::<usize>().context("invalid chunk count")?
    };
    if !(1..=256).contains(&chunks) {
        bail!("chunk count {chunks} is outside the supported range 1..=256; increase chunk size")
    }
    eprintln!("building gzip recovery points backend=rapidgzip without scanning FASTQ IDs");
    let scan_started = Instant::now();
    let ((mut ra, oa), (mut rb, ob)) = std::thread::scope(|scope| -> Result<_> {
        let a = scope.spawn(|| build_dense(r1, threads, "R1"));
        let b = scope.spawn(|| build_dense(r2, threads, "R2"));
        Ok((
            a.join()
                .map_err(|_| anyhow::anyhow!("R1 worker panicked"))??,
            b.join()
                .map_err(|_| anyhow::anyhow!("R2 worker panicked"))??,
        ))
    })?;
    let scan_seconds = scan_started.elapsed().as_secs_f64();
    let after1 = fingerprint(r1, "r1.gzidx".into(), oa)?;
    let after2 = fingerprint(r2, "r2.gzidx".into(), ob)?;
    if (before1.size, before1.modified_ns, &before1.edge_sha256)
        != (after1.size, after1.modified_ns, &after1.edge_sha256)
        || (before2.size, before2.modified_ns, &before2.edge_sha256)
            != (after2.size, after2.modified_ns, &after2.edge_sha256)
    {
        bail!("input changed during indexing")
    }
    if oa == 0 || ob == 0 {
        bail!("empty FASTQ")
    }
    eprintln!("chunk_count={chunks} R1_decoded_bytes={oa} R2_decoded_bytes={ob}");
    let (mut boundaries, mut last1, mut last2) = (Vec::new(), None, None);
    let (mut r1_search_seconds, mut r2_search_seconds) = (0.0, 0.0);
    for i in 0..chunks {
        let a = ((oa as u128 * i as u128) / chunks as u128) as u64;
        let search_started = Instant::now();
        let (off1, id) = find_record(&mut ra, a, oa.min(a.saturating_add(window)), None)?
            .context("no valid R1 record in boundary window")?;
        r1_search_seconds += search_started.elapsed().as_secs_f64();
        let predicted = ((ob as u128 * off1 as u128) / oa as u128) as u64;
        let start2 = predicted.saturating_sub(window);
        let end2 = ob.min(predicted.saturating_add(window));
        let search_started = Instant::now();
        let (off2, _) = find_record_near(&mut rb, start2, end2, predicted, &id, 4 * 1024 * 1024)?
            .with_context(|| {
            format!(
                "R2 ID not found within ±{window} bytes near {predicted}; increase search window"
            )
        })?;
        r2_search_seconds += search_started.elapsed().as_secs_f64();
        if last1.is_some_and(|x| off1 <= x) || last2.is_some_and(|x| off2 <= x) {
            bail!("boundary offsets not strictly increasing")
        }
        last1 = Some(off1);
        last2 = Some(off2);
        boundaries.push(Boundary {
            pair: None,
            r1_uncompressed_offset: off1,
            r2_uncompressed_offset: off2,
            read_id: String::from_utf8(id)?,
            r1_checkpoint: None,
            r2_checkpoint: None,
        });
    }
    if boundaries[0].r1_uncompressed_offset != 0 || boundaries[0].r2_uncompressed_offset != 0 {
        bail!("first read pair must match at start of both files")
    }
    boundaries.push(Boundary {
        pair: None,
        r1_uncompressed_offset: oa,
        r2_uncompressed_offset: ob,
        read_id: String::new(),
        r1_checkpoint: None,
        r2_checkpoint: None,
    });
    let export_started = Instant::now();
    for (reader, path, offsets) in [
        (
            &mut ra,
            out.join("r1.gzidx"),
            boundaries
                .iter()
                .map(|b| b.r1_uncompressed_offset)
                .collect::<Vec<_>>(),
        ),
        (
            &mut rb,
            out.join("r2.gzidx"),
            boundaries
                .iter()
                .map(|b| b.r2_uncompressed_offset)
                .collect::<Vec<_>>(),
        ),
    ] {
        let temp = path.with_extension("gzidx.pending");
        export_sparse(reader, &temp, &offsets)?;
        fs::rename(temp, path)?;
    }
    let c1 = checkpoints(&out.join("r1.gzidx"))?;
    let c2 = checkpoints(&out.join("r2.gzidx"))?;
    let export_seconds = export_started.elapsed().as_secs_f64();
    for b in &mut boundaries {
        b.r1_checkpoint = nearest(&c1, b.r1_uncompressed_offset);
        b.r2_checkpoint = nearest(&c2, b.r2_uncompressed_offset);
    }
    let manifest=Manifest{version:1,description:"Gzip recovery points plus approximate uncompressed-byte-proportional boundaries; only boundary IDs matched.".into(),inputs:vec![after1,after2],pairs:None,chunks,boundaries};
    let tmp = out.join("manifest.json.tmp");
    fs::write(&tmp, serde_json::to_vec_pretty(&manifest)?)?;
    fs::rename(tmp, out.join("manifest.json"))?;
    eprintln!("index completed in {:.2}s", start.elapsed().as_secs_f64());
    let verify_started = Instant::now();
    verify(out, threads)?;
    let verify_seconds = verify_started.elapsed().as_secs_f64();
    eprintln!(
        "event=index_timing scan_seconds={scan_seconds:.3} r1_boundary_seconds={r1_search_seconds:.3} r2_boundary_seconds={r2_search_seconds:.3} export_seconds={export_seconds:.3} verify_seconds={verify_seconds:.3} total_seconds={:.3}",
        start.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{find_record, find_record_near};
    use std::io::Cursor;

    #[test]
    fn boundary_search_stops_at_first_match() {
        let mut data = b"@target/2 comment\nX?\n+\nI\n@target/2\nA\n+\nI\n".to_vec();
        data.resize(3 * 1024 * 1024, b'x');
        let mut input = Cursor::new(data);
        let result = find_record(&mut input, 0, 3 * 1024 * 1024, Some(b"@target"))
            .unwrap()
            .unwrap();
        assert_eq!(result, (0, b"@target".to_vec()));
        // The buffered reader may prefetch, but must not scan the whole window.
        assert!(input.position() <= 1024 * 1024);
    }

    #[test]
    fn boundary_search_handles_offsets_and_missing_ids() {
        let data = b"@first/1\nAC\n+\nII\n@second/1 comment\nTG\n+\nII\n";
        let offset = data.windows(7).position(|w| w == b"@second").unwrap() as u64;
        let mut input = Cursor::new(data.as_slice());
        assert_eq!(
            find_record(&mut input, 3, data.len() as u64, None).unwrap(),
            Some((offset, b"@second".to_vec()))
        );
        let mut input = Cursor::new(data.as_slice());
        assert_eq!(
            find_record(&mut input, offset, data.len() as u64, Some(b"@second")).unwrap(),
            Some((offset, b"@second".to_vec()))
        );
        let mut input = Cursor::new(data.as_slice());
        assert!(
            find_record(&mut input, 0, data.len() as u64, Some(b"@missing"))
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn outward_search_covers_both_sides_and_exact_probe_boundaries() {
        let mut data = Vec::new();
        let mut records = Vec::new();
        for i in 0..12 {
            let offset = data.len() as u64;
            let id = format!("@read{i:02}");
            data.extend_from_slice(format!("{id}/2 comment\nACGT\n+\nIIII\n").as_bytes());
            records.push((offset, id));
        }
        let end = data.len() as u64;
        for (offset, id) in &records {
            for predicted in [0, offset.saturating_sub(1), *offset, offset + 1, end] {
                for step in [1, 7, 31] {
                    let mut input = Cursor::new(data.as_slice());
                    assert_eq!(
                        find_record_near(&mut input, 0, end, predicted, id.as_bytes(), step)
                            .unwrap(),
                        Some((*offset, id.as_bytes().to_vec())),
                        "offset={offset} predicted={predicted} step={step}"
                    );
                }
            }
        }
        let (offset, id) = &records[5];
        let mut input = Cursor::new(data.as_slice());
        assert!(
            find_record_near(&mut input, 0, *offset, 0, id.as_bytes(), 7)
                .unwrap()
                .is_none()
        );
        let mut input = Cursor::new(data.as_slice());
        assert_eq!(
            find_record_near(&mut input, *offset, end, end, id.as_bytes(), 7).unwrap(),
            Some((*offset, id.as_bytes().to_vec()))
        );
        let mut input = Cursor::new(data.as_slice());
        assert!(
            find_record_near(&mut input, 0, end, end / 2, b"@missing", 7)
                .unwrap()
                .is_none()
        );
    }
}
