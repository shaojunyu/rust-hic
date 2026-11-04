use super::format::{read_header, Format};
use anyhow::{bail, Context, Result};
use std::{
    collections::{BTreeMap, HashMap, HashSet, VecDeque},
    env,
    io::{self, BufRead, BufReader, BufWriter, Write},
    ops::Range,
    sync::mpsc,
    thread,
    time::Instant,
};

type Group = (String, String, u8, u8);
#[derive(Clone)]
struct Pair {
    line: String,
    c1: Range<usize>,
    c2: Range<usize>,
    p1: i64,
    p2: i64,
    s1: u8,
    s2: u8,
    kind: Range<usize>,
}
impl Pair {
    fn parse(line: String, format: Format) -> Result<Self> {
        let f = format.fields(&line)?;
        let p1: i64 = f[2].parse().context("pos1")?;
        let p2: i64 = f[4].parse().context("pos2")?;
        if f[1] == "!"
            || f[3] == "!"
            || p1 <= 0
            || p2 <= 0
            || p1 == i64::MAX
            || p2 == i64::MAX
            || !["+", "-"].contains(&f[5])
            || !["+", "-"].contains(&f[6])
        {
            bail!("unsupported pair")
        }
        // The fields borrow from `line`; keep their byte ranges instead of
        // allocating separate chromosome and pair-type strings per record.
        let range = |field: &str| {
            let start = field.as_ptr() as usize - line.as_ptr() as usize;
            start..start + field.len()
        };
        Ok(Self {
            c1: range(f[1]),
            c2: range(f[3]),
            p1,
            p2,
            s1: f[5].as_bytes()[0],
            s2: f[6].as_bytes()[0],
            kind: range(f[7]),
            line,
        })
    }
    fn group(&self) -> Group {
        (
            self.line[self.c1.clone()].to_owned(),
            self.line[self.c2.clone()].to_owned(),
            self.s1,
            self.s2,
        )
    }
    fn key(&self) -> (&str, &str, i64, i64, &str) {
        (
            &self.line[self.c1.clone()],
            &self.line[self.c2.clone()],
            self.p1,
            self.p2,
            &self.line[self.kind.clone()],
        )
    }
}
fn root(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i]
    }
    i
}
fn union(parent: &mut [usize], a: usize, b: usize) {
    let x = root(parent, a);
    let y = root(parent, b);
    if x != y {
        parent[x] = y
    }
}
#[derive(Default)]
struct Coords {
    entries: HashMap<(i64, i64), usize>,
    range: Option<(i64, i64)>,
}
impl Coords {
    fn insert(&mut self, pos: (i64, i64), index: usize) {
        self.range = Some(match self.range {
            Some((min, max)) => (min.min(pos.0), max.max(pos.0)),
            None => (pos.0, pos.0),
        });
        self.entries.insert(pos, index);
    }
}
struct Batch {
    rows: Vec<Pair>,
    parent: Vec<usize>,
    groups: HashMap<Group, Coords>,
    local_output: Vec<u8>,
    local_carry: VecDeque<Pair>,
    local_duplicates: u64,
}
fn neighbors(coords: &Coords, p: &Pair, mismatch: i64) -> Vec<usize> {
    let Some((min, max)) = coords.range else {
        return Vec::new();
    };
    if p.p1 < min.saturating_sub(mismatch) || p.p1 > max.saturating_add(mismatch) {
        return Vec::new();
    }
    if mismatch > 16 {
        return coords
            .entries
            .iter()
            .filter_map(|(&(x, y), &i)| {
                (x.abs_diff(p.p1) <= mismatch as u64 && y.abs_diff(p.p2) <= mismatch as u64)
                    .then_some(i)
            })
            .collect();
    }
    let mut found = Vec::new();
    for dx in -mismatch..=mismatch {
        for dy in -mismatch..=mismatch {
            if let (Some(x), Some(y)) = (p.p1.checked_add(dx), p.p2.checked_add(dy)) {
                if let Some(&i) = coords.entries.get(&(x, y)) {
                    found.push(i);
                }
            }
        }
    }
    found
}
fn prepare(lines: Vec<String>, format: Format, mismatch: i64) -> Result<Batch> {
    let rows: Vec<Pair> = lines
        .into_iter()
        .map(|line| Pair::parse(line, format))
        .collect::<Result<_>>()?;
    for w in rows.windows(2) {
        if w[0].key() > w[1].key() {
            bail!("input not sorted")
        }
    }
    let mut parent: Vec<usize> = (0..rows.len()).collect();
    let mut groups: HashMap<Group, Coords> = HashMap::new();
    for (i, p) in rows.iter().enumerate() {
        let coords = groups.entry(p.group()).or_default();
        for j in neighbors(coords, p, mismatch) {
            union(&mut parent, i, j)
        }
        coords.insert((p.p1, p.p2), i);
    }
    let mut local_output = Vec::new();
    let mut local_carry = VecDeque::new();
    let mut seen = HashSet::new();
    let mut local_duplicates = 0;
    for (i, p) in rows.iter().enumerate() {
        if seen.insert(root(&mut parent, i)) {
            local_output.extend_from_slice(p.line.as_bytes());
            local_carry.push_back(i);
            if local_carry.len() > 100 {
                local_carry.pop_front();
            }
        } else {
            local_duplicates += 1;
        }
    }
    let local_carry = local_carry.into_iter().map(|i| rows[i].clone()).collect();
    Ok(Batch {
        rows,
        parent,
        groups,
        local_output,
        local_carry,
        local_duplicates,
    })
}
fn process(
    mut batch: Batch,
    carry: &mut VecDeque<Pair>,
    out: &mut impl Write,
    mismatch: i64,
) -> Result<(u64, u64, bool)> {
    let n = batch.rows.len();
    let mut touches = false;
    'outer: for p in carry.iter() {
        if let Some(coords) = batch.groups.get(&p.group()) {
            if !neighbors(coords, p, mismatch).is_empty() {
                touches = true;
                break 'outer;
            }
        }
    }
    if !touches {
        out.write_all(&batch.local_output)?;
        *carry = batch.local_carry;
        return Ok((n as u64, batch.local_duplicates, false));
    }
    for (k, p) in carry.iter().enumerate() {
        let i = n + k;
        batch.parent.push(i);
        let coords = batch.groups.entry(p.group()).or_default();
        for j in neighbors(coords, p, mismatch) {
            union(&mut batch.parent, i, j)
        }
        coords.insert((p.p1, p.p2), i);
    }
    let mut seen = HashSet::new();
    for i in n..batch.parent.len() {
        seen.insert(root(&mut batch.parent, i));
    }
    carry.clear();
    let mut duplicates = 0;
    for (i, p) in batch.rows.into_iter().enumerate() {
        if seen.insert(root(&mut batch.parent, i)) {
            out.write_all(p.line.as_bytes())?;
            carry.push_back(p);
            if carry.len() > 100 {
                carry.pop_front();
            }
        } else {
            duplicates += 1
        }
    }
    Ok((n as u64, duplicates, true))
}
fn run() -> Result<()> {
    let args: Vec<_> = env::args().collect();
    if args.len() != 3 {
        bail!("usage: dedup_pairs_pipeline THREADS MAX_MISMATCH < sorted.pairs > nodups.pairs")
    }
    let mismatch: i64 = args[2].parse().context("invalid max mismatch")?;
    if mismatch < 0 {
        bail!("max mismatch must be non-negative")
    }
    let threads: usize = args[1].parse()?;
    if !(1..=64).contains(&threads) {
        bail!("threads must be 1..64")
    }
    let start = Instant::now();
    let (job_tx, job_rx) = crossbeam_channel::bounded::<(usize, Vec<String>, Format)>(threads * 2);
    let (res_tx, res_rx) = mpsc::sync_channel::<(usize, Result<Batch>)>(threads * 2);
    let (permit_tx, permit_rx) = mpsc::sync_channel::<()>(threads * 4);
    for _ in 0..threads * 4 {
        permit_tx.send(())?;
    }
    let (head_tx, head_rx) = mpsc::sync_channel(1);
    let reader = thread::spawn(move || -> Result<(f64, f64, f64)> {
        let stdin = io::stdin();
        let mut input = BufReader::with_capacity(1024 * 1024, stdin.lock());
        let (headers, format) = read_header(&mut input)?;
        if !headers
            .iter()
            .any(|l| l.starts_with("#sorted: chr1-chr2-pos1-pos2"))
        {
            bail!("missing sorted pairs header")
        }
        head_tx.send(headers)?;
        let mut id = 0;
        let (mut permit_wait, mut read_time, mut send_wait) = (0.0, 0.0, 0.0);
        loop {
            let stamp = Instant::now();
            permit_rx.recv()?;
            permit_wait += stamp.elapsed().as_secs_f64();
            let stamp = Instant::now();
            let mut lines = Vec::with_capacity(10_000);
            for _ in 0..10_000 {
                let mut line = String::new();
                if input.read_line(&mut line)? == 0 {
                    break;
                }
                lines.push(line)
            }
            read_time += stamp.elapsed().as_secs_f64();
            if lines.is_empty() {
                break;
            }
            let stamp = Instant::now();
            job_tx.send((id, lines, format))?;
            send_wait += stamp.elapsed().as_secs_f64();
            id += 1
        }
        Ok((permit_wait, read_time, send_wait))
    });
    let mut workers = Vec::new();
    for _ in 0..threads {
        let rx = job_rx.clone();
        let tx = res_tx.clone();
        workers.push(thread::spawn(move || {
            let (mut wait, mut compute, mut send) = (0.0, 0.0, 0.0);
            loop {
                let stamp = Instant::now();
                let job = rx.recv();
                wait += stamp.elapsed().as_secs_f64();
                match job {
                    Ok((id, lines, format)) => {
                        let stamp = Instant::now();
                        let batch = prepare(lines, format, mismatch);
                        compute += stamp.elapsed().as_secs_f64();
                        let stamp = Instant::now();
                        let sent = tx.send((id, batch));
                        send += stamp.elapsed().as_secs_f64();
                        if sent.is_err() {
                            break;
                        }
                    }
                    Err(_) => break,
                }
            }
            (wait, compute, send)
        }));
    }
    drop(res_tx);
    let mut out = BufWriter::with_capacity(1024 * 1024, io::stdout().lock());
    for h in head_rx.recv().context("header reader failed")? {
        out.write_all(h.as_bytes())?
    }
    let mut pending = BTreeMap::new();
    let mut next = 0;
    let mut carry = VecDeque::new();
    let (mut total, mut dup) = (0, 0);
    let mut reconciled = 0u64;
    let mut previous: Option<Pair> = None;
    let (mut result_wait, mut ordered_time) = (0.0, 0.0);
    loop {
        let stamp = Instant::now();
        let message = res_rx.recv();
        result_wait += stamp.elapsed().as_secs_f64();
        let Ok((id, result)) = message else { break };
        pending.insert(id, result);
        while let Some(result) = pending.remove(&next) {
            let stamp = Instant::now();
            let batch = result?;
            if let (Some(prev), Some(first)) = (&previous, batch.rows.first()) {
                if prev.key() > first.key() {
                    bail!("input not sorted")
                }
            }
            previous = batch.rows.last().cloned();
            let (n, d, touched) = process(batch, &mut carry, &mut out, mismatch)?;
            reconciled += u64::from(touched);
            total += n;
            dup += d;
            next += 1;
            let _ = permit_tx.send(());
            ordered_time += stamp.elapsed().as_secs_f64();
            if next % 1000 == 0 {
                eprintln!(
                    "event=dedup_progress records={total} duplicates={dup} elapsed_seconds={:.3}",
                    start.elapsed().as_secs_f64()
                )
            }
        }
    }
    let (permit_wait, read_time, send_wait) = reader
        .join()
        .map_err(|_| anyhow::anyhow!("reader panicked"))??;
    let (mut worker_wait, mut worker_compute, mut worker_send) = (0.0, 0.0, 0.0);
    for w in workers {
        let (a, b, c) = w.join().map_err(|_| anyhow::anyhow!("worker panicked"))?;
        worker_wait += a;
        worker_compute += b;
        worker_send += c;
    }
    out.flush()?;
    eprintln!("event=dedup_done input_records={total} duplicates={dup} output_records={} elapsed_seconds={:.3}",total-dup,start.elapsed().as_secs_f64());
    eprintln!("event=timing reader_permit_wait={permit_wait:.3} reader_read={read_time:.3} reader_send_wait={send_wait:.3} worker_wait_sum={worker_wait:.3} worker_compute_sum={worker_compute:.3} worker_send_wait_sum={worker_send:.3} ordered_wait={result_wait:.3} ordered_process={ordered_time:.3}");
    eprintln!("event=boundary_stats batches={next} reconciled_batches={reconciled}");
    Ok(())
}
pub(crate) fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e:#}");
        std::process::exit(1)
    }
}
