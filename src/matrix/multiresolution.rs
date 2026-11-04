use super::attr::{int_attr, str_attr};
use anyhow::{bail, Context, Result};
use hdf5::{
    types::{EnumMember, EnumType, FixedAscii, IntSize, TypeDescriptor},
    File, Group,
};
use rayon::prelude::*;
use std::{
    env,
    fs::{self, File as StdFile},
    io::{BufRead, BufReader},
    path::Path,
    time::Instant,
};
#[derive(Clone)]
struct Chrom {
    name: String,
    len: u64,
}
struct Pixels {
    x: Vec<i64>,
    y: Vec<i64>,
    v: Vec<i32>,
}
fn dbl_attr(d: &hdf5::Dataset, k: &str, v: f64) -> Result<()> {
    d.new_attr::<f64>().shape(()).create(k)?.write_scalar(&v)?;
    Ok(())
}
fn bool_attr(d: &hdf5::Dataset, k: &str, v: bool) -> Result<()> {
    d.new_attr::<bool>().shape(()).create(k)?.write_scalar(&v)?;
    Ok(())
}
fn median(mut data: Vec<f64>) -> f64 {
    if data.is_empty() {
        return f64::NAN;
    }
    data.sort_unstable_by(f64::total_cmp);
    let n = data.len();
    if n % 2 == 1 {
        data[n / 2]
    } else {
        (data[n / 2 - 1] + data[n / 2]) / 2.0
    }
}
fn margins(
    p: &Pixels,
    bias: Option<&[f64]>,
    binary: bool,
    nbins: usize,
    threads: usize,
) -> Vec<f64> {
    let n = p.x.len();
    let part = n.div_ceil(threads);
    let pieces: Vec<Vec<f64>> = (0..threads)
        .into_par_iter()
        .map(|i| {
            let mut out = vec![0f64; nbins];
            let lo = (i * part).min(n);
            let hi = ((i + 1) * part).min(n);
            for j in lo..hi {
                let a = p.x[j] as usize;
                let b = p.y[j] as usize;
                if b - a < 2 {
                    continue;
                }
                let val = if binary { 1.0 } else { p.v[j] as f64 };
                let v = if let Some(w) = bias {
                    w[a] * w[b] * val
                } else {
                    val
                };
                out[a] += v;
                out[b] += v;
            }
            out
        })
        .collect();
    let mut result = vec![0f64; nbins];
    for piece in pieces {
        for (a, b) in result.iter_mut().zip(piece) {
            *a += b
        }
    }
    result
}
struct Balance {
    weights: Vec<f64>,
    scale: f64,
    var: f64,
    converged: bool,
    iters: usize,
}
fn balance(
    p: &Pixels,
    offsets: &[i64],
    threads: usize,
    start: &Instant,
    res: u64,
    max_iters: usize,
) -> Balance {
    let n = *offsets.last().unwrap() as usize;
    let nnz = margins(p, None, true, n, threads);
    let raw = margins(p, None, false, n, threads);
    let mut bias: Vec<f64> = nnz
        .iter()
        .map(|x| if *x < 10.0 { 0.0 } else { 1.0 })
        .collect();
    let mut normalized = raw.clone();
    for w in offsets.windows(2) {
        let (lo, hi) = (w[0] as usize, w[1] as usize);
        let m = median(raw[lo..hi].iter().copied().filter(|x| *x > 0.0).collect());
        for v in &mut normalized[lo..hi] {
            *v /= m
        }
    }
    let logs: Vec<f64> = normalized
        .iter()
        .copied()
        .filter(|x| *x > 0.0)
        .map(f64::ln)
        .collect();
    let med = median(logs.clone());
    let mad = median(logs.into_iter().map(|x| (x - med).abs()).collect());
    let cutoff = (med - 5.0 * mad).exp();
    for (i, v) in normalized.iter().enumerate() {
        if *v < cutoff {
            bias[i] = 0.0
        }
    }
    let (mut scale, mut var, mut iters) = (f64::NAN, 0.0, 0usize);
    for iter in 0..max_iters {
        let marg = margins(p, Some(&bias), false, n, threads);
        let vals: Vec<f64> = marg.iter().copied().filter(|x| *x != 0.0).collect();
        iters = iter + 1;
        if vals.is_empty() {
            bias.fill(f64::NAN);
            scale = f64::NAN;
            var = 0.0;
            break;
        }
        scale = vals.iter().sum::<f64>() / (vals.len() as f64);
        var = vals.iter().map(|x| (x - scale) * (x - scale)).sum::<f64>() / (vals.len() as f64);
        for (w, m) in bias.iter_mut().zip(marg) {
            if m != 0.0 {
                *w /= m / scale
            }
        }
        if iter % 10 == 0 {
            eprintln!("event=balance_progress resolution={res} iteration={} variance={var:.8} elapsed_seconds={:.3}",iter+1,start.elapsed().as_secs_f64())
        }
        if var < 1e-5 {
            break;
        }
    }
    for w in &mut bias {
        if *w == 0.0 {
            *w = f64::NAN
        } else {
            *w /= scale.sqrt()
        }
    }
    Balance {
        weights: bias,
        scale,
        var,
        converged: var < 1e-5,
        iters,
    }
}
#[allow(clippy::too_many_arguments)] // Group output combines data and balancing parameters.
fn write_group(
    g: &Group,
    chroms: &[Chrom],
    resolution: u64,
    p: &Pixels,
    threads: usize,
    start: &Instant,
    assembly: &str,
    max_iters: usize,
) -> Result<()> {
    let cg = g.create_group("chroms")?;
    let bg = g.create_group("bins")?;
    let pg = g.create_group("pixels")?;
    let ig = g.create_group("indexes")?;
    let names: Vec<FixedAscii<256>> = chroms
        .iter()
        .map(|c| FixedAscii::from_ascii(&c.name))
        .collect::<std::result::Result<_, _>>()?;
    let lens: Vec<i32> = chroms.iter().map(|c| c.len as i32).collect();
    cg.new_dataset_builder().with_data(&names).create("name")?;
    cg.new_dataset_builder().with_data(&lens).create("length")?;
    let (mut bc, mut bs, mut be, mut co) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for (i, c) in chroms.iter().enumerate() {
        co.push(bc.len() as i64);
        let mut s = 0;
        while s < c.len {
            bc.push(i as i32);
            bs.push(s as i32);
            be.push((s + resolution).min(c.len) as i32);
            s += resolution
        }
    }
    co.push(bc.len() as i64);
    let chromtype = TypeDescriptor::Enum(EnumType {
        size: IntSize::U4,
        signed: true,
        members: chroms
            .iter()
            .enumerate()
            .map(|(i, c)| EnumMember {
                name: c.name.clone(),
                value: i as u64,
            })
            .collect(),
    });
    let d = bg
        .new_dataset_builder()
        .empty_as(&chromtype)
        .shape(bc.len())
        .create("chrom")?;
    let dtype = d.dtype()?;
    let status = hdf5::sync::sync(|| unsafe {
        hdf5_sys::h5d::H5Dwrite(d.id(), dtype.id(), 0, 0, 0, bc.as_ptr().cast())
    });
    if status < 0 {
        bail!("enum write failed")
    }
    bg.new_dataset_builder().with_data(&bs).create("start")?;
    bg.new_dataset_builder().with_data(&be).create("end")?;
    ig.new_dataset_builder()
        .with_data(&co)
        .create("chrom_offset")?;
    pg.new_dataset_builder().with_data(&p.x).create("bin1_id")?;
    pg.new_dataset_builder().with_data(&p.y).create("bin2_id")?;
    pg.new_dataset_builder().with_data(&p.v).create("count")?;
    let mut offsets = vec![0i64; bc.len() + 1];
    let mut cur = 0;
    for (i, &x) in p.x.iter().enumerate() {
        while cur <= x as usize {
            offsets[cur] = i as i64;
            cur += 1
        }
    }
    while cur < offsets.len() {
        offsets[cur] = p.x.len() as i64;
        cur += 1
    }
    ig.new_dataset_builder()
        .with_data(&offsets)
        .create("bin1_offset")?;
    let sum = p.v.iter().map(|x| *x as i64).sum::<i64>();
    for (k, v) in [
        ("format", "HDF5::Cooler"),
        ("bin-type", "fixed"),
        ("storage-mode", "symmetric-upper"),
        ("genome-assembly", assembly),
        ("generated-by", "rust-hic Rust"),
        ("metadata", "{}"),
    ] {
        str_attr(g, k, v)?
    }
    for (k, v) in [
        ("format-version", 3),
        ("bin-size", resolution as i64),
        ("nchroms", chroms.len() as i64),
        ("nbins", bc.len() as i64),
        ("nnz", p.x.len() as i64),
        ("sum", sum),
    ] {
        int_attr(g, k, v)?
    }
    eprintln!(
        "event=balance_start resolution={resolution} pixels={} elapsed_seconds={:.3}",
        p.x.len(),
        start.elapsed().as_secs_f64()
    );
    let w = balance(p, &co, threads, start, resolution, max_iters);
    let ds = bg
        .new_dataset_builder()
        .with_data(&w.weights)
        .chunk((65536.min(bc.len().max(1)),))
        .deflate(6)
        .create("weight")?;
    dbl_attr(&ds, "tol", 1e-5)?;
    dbl_attr(&ds, "min_nnz", 10.0)?;
    dbl_attr(&ds, "min_count", 0.0)?;
    dbl_attr(&ds, "mad_max", 5.0)?;
    dbl_attr(&ds, "ignore_diags", 2.0)?;
    dbl_attr(&ds, "scale", w.scale)?;
    dbl_attr(&ds, "var", w.var)?;
    bool_attr(&ds, "cis_only", false)?;
    bool_attr(&ds, "converged", w.converged)?;
    bool_attr(&ds, "divisive_weights", false)?;
    dbl_attr(&ds, "iterations", w.iters as f64)?;
    dbl_attr(&ds, "max_iters", max_iters as f64)?;
    eprintln!("event=balance_done resolution={resolution} iterations={} converged={} variance={} elapsed_seconds={:.3}",w.iters,w.converged,w.var,start.elapsed().as_secs_f64());
    Ok(())
}
fn coarsen(base: &Pixels, chroms: &[Chrom], base_res: u64, res: u64) -> Result<Pixels> {
    let mut map = Vec::new();
    let mut next = 0u64;
    for c in chroms {
        let n = c.len.div_ceil(base_res);
        for i in 0..n {
            map.push(next + (i * base_res) / res)
        }
        next += c.len.div_ceil(res)
    }
    if next > u32::MAX as u64 {
        bail!("too many bins")
    }
    let mut rows: Vec<(u64, i32)> = (0..base.x.len())
        .into_par_iter()
        .map(|i| {
            (
                (map[base.x[i] as usize] << 32) | map[base.y[i] as usize],
                base.v[i],
            )
        })
        .collect();
    rows.par_sort_unstable_by_key(|x| x.0);
    let (mut x, mut y, mut v) = (Vec::new(), Vec::new(), Vec::new());
    let mut i = 0;
    while i < rows.len() {
        let key = rows[i].0;
        let mut count = 0i64;
        while i < rows.len() && rows[i].0 == key {
            count += rows[i].1 as i64;
            i += 1
        }
        x.push((key >> 32) as i64);
        y.push((key & 0xffff_ffff) as i64);
        v.push(i32::try_from(count).context("pixel count overflow")?)
    }
    Ok(Pixels { x, y, v })
}
fn run() -> Result<()> {
    let a: Vec<_> = env::args().collect();
    if a.len() != 7 {
        bail!(
            "usage: build_mcool BASE.cool CHROMSIZES OUTPUT.mcool ASSEMBLY RESOLUTIONS_CSV THREADS"
        )
    };
    let start = Instant::now();
    let threads: usize = a[6].parse()?;
    if !(1..=64).contains(&threads) {
        bail!("invalid threads")
    };
    let max_iters: usize = env::var(crate::matrix::BALANCE_MAX_ITERS_ENV)
        .unwrap_or_else(|_| "1000".into())
        .parse()?;
    if !(1..=10000).contains(&max_iters) {
        bail!("invalid RUST_HIC_BALANCE_MAX_ITERS")
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;
    let mut chroms = Vec::new();
    for l in BufReader::new(StdFile::open(&a[2])?).lines() {
        let l = l?;
        let f: Vec<_> = l.split_whitespace().collect();
        if f.len() != 2 {
            bail!("invalid chromsizes")
        };
        let len: u64 = f[1].parse()?;
        if len == 0 || len > i32::MAX as u64 {
            bail!("invalid chromosome length")
        };
        chroms.push(Chrom {
            name: f[0].into(),
            len,
        })
    }
    if chroms.is_empty() {
        bail!("empty chromsizes")
    }
    let input = File::open(&a[1])?;
    let base_res: i64 = input.attr("bin-size")?.read_scalar()?;
    if base_res <= 0 {
        bail!("invalid base binsize")
    };
    let base_res = base_res as u64;
    let orig_names: Vec<FixedAscii<256>> = input.dataset("chroms/name")?.read_raw()?;
    let orig_lens: Vec<i32> = input.dataset("chroms/length")?.read_raw()?;
    if orig_names.len() != chroms.len()
        || orig_lens.len() != chroms.len()
        || chroms
            .iter()
            .enumerate()
            .any(|(i, c)| orig_names[i].as_str() != c.name || orig_lens[i] as u64 != c.len)
    {
        bail!("base chromosome table differs from chromsizes")
    }
    let resolutions: Vec<u64> = a[5]
        .split(',')
        .map(|s| s.parse::<u64>())
        .collect::<std::result::Result<_, _>>()?;
    if resolutions.is_empty()
        || resolutions
            .iter()
            .any(|r| *r < base_res || r % base_res != 0)
    {
        bail!("resolutions must be positive multiples of base resolution")
    };
    let mut sorted = resolutions.clone();
    sorted.sort_unstable();
    sorted.dedup();
    if sorted.len() != resolutions.len() {
        bail!("duplicate resolutions")
    }
    let output = Path::new(&a[3]);
    if output.exists() {
        bail!("output exists")
    };
    let temp = output.with_extension(format!("mcool.partial.{}", std::process::id()));
    if temp.exists() {
        bail!("partial output exists")
    };
    let base = Pixels {
        x: input.dataset("pixels/bin1_id")?.read_raw()?,
        y: input.dataset("pixels/bin2_id")?.read_raw()?,
        v: input.dataset("pixels/count")?.read_raw()?,
    };
    if base.x.len() != base.y.len() || base.x.len() != base.v.len() {
        bail!("pixel lengths differ")
    };
    let nbins: u64 = chroms.iter().map(|c| c.len.div_ceil(base_res)).sum();
    let mut prior: Option<(i64, i64)> = None;
    for ((&x, &y), &count) in base.x.iter().zip(&base.y).zip(&base.v) {
        if x < 0 || y < x || y as u64 >= nbins || count <= 0 || prior.is_some_and(|p| p >= (x, y)) {
            bail!("base pixels must be positive, symmetric-upper, in range and strictly sorted")
        }
        prior = Some((x, y));
    }
    eprintln!(
        "event=base_loaded pixels={} elapsed_seconds={:.3}",
        base.x.len(),
        start.elapsed().as_secs_f64()
    );
    let file = File::create(&temp)?;
    let root = file.group("/")?;
    str_attr(&root, "format", "HDF5::MCOOL")?;
    int_attr(&root, "format-version", 2)?;
    let gs = root.create_group("resolutions")?;
    for res in sorted {
        let g = gs.create_group(&res.to_string())?;
        if res == base_res {
            write_group(&g, &chroms, res, &base, threads, &start, &a[4], max_iters)?
        } else {
            eprintln!(
                "event=coarsen_start resolution={res} elapsed_seconds={:.3}",
                start.elapsed().as_secs_f64()
            );
            let p = coarsen(&base, &chroms, base_res, res)?;
            write_group(&g, &chroms, res, &p, threads, &start, &a[4], max_iters)?
        }
        file.flush()?
    }
    drop(gs);
    drop(root);
    drop(file);
    fs::hard_link(&temp, output).context("publish without overwriting")?;
    fs::remove_file(temp)?;
    eprintln!(
        "event=mcool_done elapsed_seconds={:.3}",
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
pub(crate) fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e:#}");
        std::process::exit(1)
    }
}
