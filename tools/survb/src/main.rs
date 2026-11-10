//! survb: trimming (Lee 2009) selection bound per perturbation on a sparse single-cell cache.
//!
//! For each perturbation p with unit direction u_p (gene weights) and surviving fraction S_p:
//!   score of a cell c on u_p:  y = sum_j u_pj * ln(1 + CP10k_cj)
//!   m      = mean(y | cells of p) - mean(y | control cells)
//!   d_top  = mean of the top ceil(S n) control scores - mean of all control scores  (>= 0)
//!   d_bot  = mean of the bottom ceil(S n) control scores - mean of all control scores (<= 0)
//! Selection alone (keeping a fraction S of cells that looked like controls) can move the mean by at most d_top
//! upwards or |d_bot| downwards. Control scores for all perturbations are a sparse x dense product, computed in
//! blocks of perturbations, parallel over cells (rayon), with the dense block row-major so the inner loop is a
//! contiguous axpy that the compiler vectorizes with NEON.
//!
//! usage: survb pb <cache_dir> <out_prefix>   (pseudobulk means: all cells and two interleaved halves)
//!        survb <cache_dir> <dirs.tsv: pert gene_index weight> <pert.tsv: pert n S chronos> <out.tsv> [half]
//! half = 0 or 1 keeps only the perturbed cells at even or odd positions of the group member list (the same
//! interleaved halves as the pseudobulk tool), so a direction built on one half is measured on the other.
mod cache;
mod pb;
use anyhow::{bail, Result};
use cache::Cache;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::time::Instant;

const BLOCK: usize = 256;
const NONE: u32 = u32::MAX;

fn read_tsv(p: &str) -> Result<Vec<Vec<String>>> {
    Ok(std::fs::read_to_string(p)?
        .lines()
        .skip(1)
        .map(|l| l.split('\t').map(String::from).collect())
        .collect())
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    if a.len() == 4 && a[1] == "pb" {
        return pb::run(&a[2], &a[3]);
    }
    if a.len() != 5 && a.len() != 6 {
        bail!("usage: survb <cache_dir> <dirs.tsv> <pert.tsv> <out.tsv> [half]");
    }
    let half: Option<usize> = if a.len() == 6 { Some(a[5].parse()?) } else { None };
    let t0 = Instant::now();
    let c = Cache::open(&a[1])?;
    let ng = c.genes.len();
    // perturbations, in pert.tsv order
    let pert_rows = read_tsv(&a[3])?;
    let perts: Vec<String> = pert_rows.iter().map(|r| r[0].clone()).collect();
    let s_frac: Vec<f64> = pert_rows.iter().map(|r| r[2].parse().unwrap()).collect();
    let np = perts.len();
    let pidx: HashMap<&str, usize> = perts.iter().enumerate().map(|(i, p)| (p.as_str(), i)).collect();
    // dense U over the union of direction genes: u[row * np + p]
    let mut gmap = vec![NONE; ng];
    let mut nrow = 0usize;
    let dirs = read_tsv(&a[2])?;
    for r in &dirs {
        let g: usize = r[1].parse()?;
        if gmap[g] == NONE {
            gmap[g] = nrow as u32;
            nrow += 1;
        }
    }
    let mut u = vec![0f32; nrow * np];
    for r in &dirs {
        let p = pidx[r[0].as_str()];
        let g: usize = r[1].parse()?;
        u[gmap[g] as usize * np + p] = r[2].parse()?;
    }
    let gidx: HashMap<&str, usize> = c.groups.iter().enumerate().map(|(i, g)| (g.0.as_str(), i)).collect();
    let ctrl = c.cells(gidx["control"]);
    let nc = ctrl.len();
    eprintln!("survb: {np} perturbations, {nrow} direction genes, {nc} control cells ({:.1}s)", t0.elapsed().as_secs_f64());

    // ---- control scores, block by block
    let mut stats = vec![(0f64, 0f64, 0f64); np]; // mean, d_top, d_bot
    for b0 in (0..np).step_by(BLOCK) {
        let bw = BLOCK.min(np - b0);
        let mut sc = vec![0f32; nc * bw];
        sc.par_chunks_mut(bw).zip(ctrl.par_iter()).for_each(|(out, &cell)| {
            let cell = cell as usize;
            let s = c.scale[cell];
            let (ii, vv) = c.row(cell);
            for (&i, &v) in ii.iter().zip(vv) {
                let r = gmap[i as usize];
                if r == NONE {
                    continue;
                }
                let y = (v as f32 * s).ln_1p();
                let ur = &u[r as usize * np + b0..r as usize * np + b0 + bw];
                for (o, &w) in out.iter_mut().zip(ur) {
                    *o += y * w;
                }
            }
        });
        let st: Vec<(f64, f64, f64)> = (0..bw)
            .into_par_iter()
            .map(|j| {
                let mut col: Vec<f32> = (0..nc).map(|i| sc[i * bw + j]).collect();
                col.sort_unstable_by(|x, y| x.partial_cmp(y).unwrap());
                let mean = col.iter().map(|&x| x as f64).sum::<f64>() / nc as f64;
                let k = ((s_frac[b0 + j] * nc as f64).ceil() as usize).clamp(1, nc);
                let top = col[nc - k..].iter().map(|&x| x as f64).sum::<f64>() / k as f64;
                let bot = col[..k].iter().map(|&x| x as f64).sum::<f64>() / k as f64;
                (mean, top - mean, bot - mean)
            })
            .collect();
        stats[b0..b0 + bw].copy_from_slice(&st);
    }
    eprintln!("survb: control scores done ({:.1}s)", t0.elapsed().as_secs_f64());

    // ---- perturbed cells, each on its own direction
    let pres: Vec<(usize, f64, f64)> = (0..np)
        .into_par_iter()
        .map(|p| {
            let Some(&g) = gidx.get(perts[p].as_str()) else { return (0, f64::NAN, f64::NAN) };
            let cells = c.cells(g);
            let ys: Vec<f64> = cells
                .iter()
                .enumerate()
                .filter(|(k, _)| half.map_or(true, |h| k & 1 == h))
                .map(|(_, &cell)| {
                    let cell = cell as usize;
                    let s = c.scale[cell];
                    let (ii, vv) = c.row(cell);
                    let mut y = 0f64;
                    for (&i, &v) in ii.iter().zip(vv) {
                        let r = gmap[i as usize];
                        if r != NONE {
                            y += ((v as f32 * s).ln_1p() * u[r as usize * np + p]) as f64;
                        }
                    }
                    y
                })
                .collect();
            let n = ys.len();
            let m = ys.iter().sum::<f64>() / n.max(1) as f64;
            let v = ys.iter().map(|y| (y - m) * (y - m)).sum::<f64>() / (n.max(2) - 1) as f64;
            (n, m, v.sqrt())
        })
        .collect();
    let mut w = BufWriter::new(std::fs::File::create(&a[4])?);
    writeln!(w, "pert\tn_ctrl\tn_pert\tS\tmean_ctrl\tmean_pert\tsd_pert\tm\td_top\td_bot")?;
    for p in 0..np {
        let (mean, dt, db) = stats[p];
        let (n, mp, sd) = pres[p];
        writeln!(w, "{}\t{nc}\t{n}\t{:.4}\t{mean:.6}\t{mp:.6}\t{sd:.6}\t{:.6}\t{dt:.6}\t{db:.6}", perts[p], s_frac[p], mp - mean)?;
    }
    w.flush()?;
    eprintln!("survb: done {:.1}s", t0.elapsed().as_secs_f64());
    Ok(())
}
