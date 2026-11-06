#!/usr/bin/env python3
"""Phase 1 (SERVER): selection bound for every Tahoe-100M condition (sample x line).

Reads the per-cell PCA(50) memmap written by the sample-sufficiency pipeline (pass4a_stream.py):
    $TAHOE_PASS4/coords.f32   float32 (cells x 50)
    $TAHOE_PASS4/condid.i32   int32 condition id per cell
    $TAHOE_PASS4/ckpt.npz     cursor, cond_keys ("sample|line"), cond_drug {key: drug}
For each treated condition it computes
    S      surviving fraction = cells in the condition's well (sample) / median cells in DMSO wells
           (same plate when --sample-meta gives sample->plate, otherwise all DMSO wells; column S_scope says which)
    u, m   direction and length of (treated centroid - DMSO centroid of the same line)
    m_corr bias-corrected length, m^2 - tr(Sigma_t)/N - tr(Sigma_c)/N_c, floored at 0
    d_max  trimming bound from the DMSO cells of the same line (scripts/selection_bound.py)
    ratio  d_max / m_corr   (>= 1: selection alone can explain the whole shift)
Writes outputs/selection_bound_tahoe.csv and outputs/selection_bound_tahoe_summary.json.

    TAHOE_PASS4=data/tahoe_pass4 python3 scripts/phase1_bound_tahoe.py [--sample-meta FILE.csv]
Memory: about 0.5 GB for the DMSO cells plus one 8M-cell chunk; one streaming pass over the 18 GB memmap.
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from selection_bound import d_max_many

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--pass4", default=os.environ.get("TAHOE_PASS4", "data/tahoe_pass4"))
ap.add_argument("--nc", type=int, default=50)
ap.add_argument("--sample-meta", default=None, help="CSV with columns sample,plate (optional)")
ap.add_argument("--min-cells", type=int, default=50)
ap.add_argument("--out", default=os.path.join(ROOT, "outputs"))
ap.add_argument("--chunk", type=int, default=8_000_000)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
t0 = time.time()

z = np.load(os.path.join(a.pass4, "ckpt.npz"), allow_pickle=True)
cursor = int(z["cursor"]); keys = [str(k) for k in z["cond_keys"]]; cdrug = dict(z["cond_drug"].item())
K, NC = len(keys), a.nc
nrows = os.path.getsize(os.path.join(a.pass4, "coords.f32")) // (4 * NC)
coords = np.memmap(os.path.join(a.pass4, "coords.f32"), dtype=np.float32, mode="r", shape=(nrows, NC))[:cursor]
condid = np.memmap(os.path.join(a.pass4, "condid.i32"), dtype=np.int32, mode="r", shape=(cursor,))
meta = pd.DataFrame({"cid": range(K), "key": keys})
meta["sample"] = meta.key.str.rsplit("|", n=1).str[0]; meta["line"] = meta.key.str.rsplit("|", n=1).str[1]
meta["drug"] = meta.key.map(lambda k: cdrug.get(k, "?")); meta["dmso"] = meta.drug.str.startswith("DMSO")
print(f"[phase1] {cursor:,} cells, {K} conditions, {meta.dmso.sum()} DMSO", flush=True)

# pass 1: per-condition counts, sums, sums of squares (for centroids and trace of Sigma)
cnt = np.zeros(K); sm = np.zeros((K, NC)); ssq = np.zeros(K)
for s in range(0, cursor, a.chunk):
    e = min(s + a.chunk, cursor); cc = np.asarray(condid[s:e]); X = np.asarray(coords[s:e], dtype=np.float64)
    cnt += np.bincount(cc, minlength=K)
    for j in range(NC):
        sm[:, j] += np.bincount(cc, weights=X[:, j], minlength=K)
    ssq += np.bincount(cc, weights=(X * X).sum(1), minlength=K)
    print(f"  streamed {e:,}/{cursor:,}  ({time.time()-t0:.0f}s)", flush=True)
meta["N"] = cnt
cen = sm / np.maximum(cnt, 1)[:, None]
trS = ssq / np.maximum(cnt, 1) - (cen ** 2).sum(1)                 # tr(Sigma) per condition

# surviving fraction from cells per well
well = meta.groupby("sample").N.sum().rename("well_cells")
meta = meta.join(well, on="sample")
if a.sample_meta:
    sp = pd.read_csv(a.sample_meta)[["sample", "plate"]].astype(str)
    meta = meta.merge(sp, on="sample", how="left")
    dref = meta[meta.dmso].drop_duplicates("sample").groupby("plate").well_cells.median()
    meta["dmso_well_ref"] = meta.plate.map(dref); meta["S_scope"] = "plate"
    miss = meta.dmso_well_ref.isna()
    meta.loc[miss, "dmso_well_ref"] = meta[meta.dmso].drop_duplicates("sample").well_cells.median()
    meta.loc[miss, "S_scope"] = "global"
else:
    meta["dmso_well_ref"] = meta[meta.dmso].drop_duplicates("sample").well_cells.median(); meta["S_scope"] = "global"
meta["S"] = (meta.well_cells / meta.dmso_well_ref).clip(upper=1.0)

# DMSO reference per line: centroid, trace, and the cells themselves (for the bound)
dm_ids = meta[meta.dmso].groupby("line").cid.apply(list).to_dict()
ref_cen, ref_tr, ref_n = {}, {}, {}
for line, ids in dm_ids.items():
    n = cnt[ids].sum(); ref_n[line] = n
    ref_cen[line] = sm[ids].sum(0) / n
    ref_tr[line] = ssq[ids].sum() / n - (ref_cen[line] ** 2).sum()
want = np.isin(np.asarray(condid), meta[meta.dmso].cid.values)          # one boolean pass over cell ids
idx = np.flatnonzero(want)
dm_cells = np.asarray(coords[idx]); dm_cond = np.asarray(condid[idx])
line_of = meta.set_index("cid").line
dm_line = line_of.reindex(dm_cond).values
print(f"[phase1] gathered {len(idx):,} DMSO cells ({time.time()-t0:.0f}s)", flush=True)

rows = []
tr_meta = meta[(~meta.dmso) & (meta.N >= a.min_cells) & meta.line.isin(ref_cen.keys())]
for line, grp in tr_meta.groupby("line"):
    Xc = dm_cells[dm_line == line]
    v = cen[grp.cid.values] - ref_cen[line]
    m = np.linalg.norm(v, axis=1)
    m2c = m ** 2 - trS[grp.cid.values] / grp.N.values - ref_tr[line] / ref_n[line]
    mc = np.sqrt(np.clip(m2c, 0, None))
    U = v / np.maximum(m, 1e-12)[:, None]
    dmax = d_max_many(Xc, U, grp.S.values)
    out = grp[["sample", "line", "drug", "N", "well_cells", "S", "S_scope"]].copy()
    out["m"] = m; out["m_corr"] = mc; out["d_max"] = dmax
    out["ratio"] = np.where(mc > 0, dmax / np.maximum(mc, 1e-12), np.inf)
    rows.append(out)
res = pd.concat(rows, ignore_index=True)
res.to_csv(os.path.join(a.out, "selection_bound_tahoe.csv"), index=False)
bands = pd.cut(res.S, [0, .1, .25, .5, .75, 1.0001])
summ = dict(
    n_conditions=int(len(res)), S_scope=res.S_scope.value_counts().to_dict(),
    pct_ratio_ge_1=round(100 * float((res.ratio >= 1).mean()), 2),
    pct_ratio_ge_0_5=round(100 * float((res.ratio >= .5).mean()), 2),
    by_survival_band={str(k): dict(n=int(len(g)), median_ratio=round(float(g.ratio.median()), 3),
                                   pct_ratio_ge_1=round(100 * float((g.ratio >= 1).mean()), 1))
                      for k, g in res.groupby(bands, observed=True)},
    runtime_s=round(time.time() - t0))
json.dump(summ, open(os.path.join(a.out, "selection_bound_tahoe_summary.json"), "w"), indent=2)
print(json.dumps(summ, indent=2))
