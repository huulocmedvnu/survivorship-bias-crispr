"""Phase O2 final step for one Orion line: control cache for survb, directions and shifts from the per-target sums,
trimming bounds, share of the shift that selection could produce by Chronos bin.

usage: python3 o2_finalize.py HCT116|HEK293T
"""
import os, sys, json, glob, subprocess
import numpy as np, pandas as pd
from scipy.stats import spearmanr

LINE = sys.argv[1]
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..")
ST = os.path.join(H, "orion_o2", LINE); CA = os.path.join(ST, "cache"); os.makedirs(CA, exist_ok=True)
OUT = os.path.join(R, "outputs")
SURVB = os.path.join(R, "tools/survb/target/release/survb")
TOP, MINC = 200, 20

gm = pd.read_parquet(os.path.join(H, "meta/orion_gene_metadata.parquet")).set_index("gene_token_id").gene_name
tok = np.load(os.path.join(ST, "genes.npy")); names = gm.loc[tok].values; G = len(tok)

# ---- control cache in the survb layout (genes.txt, groups.tsv, members, indptr, lib, idx, cnt)
parts = [np.load(f) for f in sorted(glob.glob(os.path.join(ST, "ctrl", "*.npz")))]
nnz = np.array([len(p["idx"]) for p in parts]); ncell = np.array([len(p["lib"]) for p in parts])
indptr = np.zeros(ncell.sum() + 1, np.uint64)
pos, base = 1, 0
for p in parts:
    indptr[pos:pos + len(p["lib"])] = p["indptr"][1:] + base; pos += len(p["lib"]); base += len(p["idx"])
idx = np.concatenate([p["idx"] for p in parts]); cnt = np.concatenate([p["cnt"] for p in parts])
lib = np.concatenate([p["lib"] for p in parts]).astype(np.float32)
N = len(lib)
open(os.path.join(CA, "genes.txt"), "w").write("\n".join(names) + "\n")
open(os.path.join(CA, "groups.tsv"), "w").write(f"control\t0\t{N}\n")
np.arange(N, dtype=np.uint32).tofile(os.path.join(CA, "members.bin"))
indptr.tofile(os.path.join(CA, "indptr.bin")); lib.tofile(os.path.join(CA, "lib.bin"))
idx.tofile(os.path.join(CA, "idx.bin")); cnt.tofile(os.path.join(CA, "cnt.bin"))
# control means per gene: CP10k and log1p(CP10k), one batch file at a time to keep memory low
c_cp = np.zeros(G); c_lg = np.zeros(G)
for p in parts:
    rr = np.repeat(np.arange(len(p["lib"])), np.diff(p["indptr"]))
    v = p["cnt"] * 1e4 / p["lib"][rr]
    c_cp += np.bincount(p["idx"], weights=v, minlength=G); c_lg += np.bincount(p["idx"], weights=np.log1p(v), minlength=G)
c_cp /= N; c_lg /= N
print(f"{LINE}: control cache {N:,} cells, {G} genes, {len(idx):,} entries")

# ---- per-target directions and shifts
targets = json.load(open(os.path.join(ST, "targets.json")))
acc = {k: np.load(os.path.join(ST, f"acc_{k}.npy")) for k in ("cp", "lg", "n")}
ge = pd.read_csv("data/depmap/CRISPRGeneEffect.csv", index_col=0)
ge.columns = [c.split(" (")[0] for c in ge.columns]
chron = ge.loc["ACH-000971"] if LINE == "HCT116" else ge.median()
name_pos = pd.Series(range(G), index=names)
rows, dirs = [], []
for t, i in targets.items():
    n0, n1 = acc["n"][2 * i], acc["n"][2 * i + 1]
    if n0 < MINC or n1 < MINC:
        continue
    m0 = acc["cp"][2 * i] / n0; m1 = acc["cp"][2 * i + 1] / n1
    l0 = np.log2((m0 + 0.1) / (c_cp + 0.1)); l1 = np.log2((m1 + 0.1) / (c_cp + 0.1))
    mask = names != t
    sel = np.where(mask)[0][np.argsort(-np.abs(l0[mask]))[:TOP]]
    u = l0[sel] / np.linalg.norm(l0[sel])
    m = u @ (acc["lg"][2 * i + 1][sel] / n1 - c_lg[sel])               # measured on the other half
    rows.append(dict(pert=t, n=int(n0 + n1), n_h1=int(n1), m=m, rep=np.corrcoef(l0[sel], l1[sel])[0, 1],
                     chronos=chron.get(t, np.nan),
                     kd=np.log2((acc["cp"][2 * i, name_pos[t]] / n0 + acc["cp"][2 * i + 1, name_pos[t]] / n1) / 2 + 0.1)
                        - np.log2(c_cp[name_pos[t]] + 0.1) if t in name_pos.index else np.nan))
    dirs += [(t, int(j), float(w)) for j, w in zip(sel, u)]
d = pd.DataFrame(rows)
ref = d.n[d.chronos > -0.3].median()
d["S"] = (d.n / ref).clip(upper=1)
keep = set(d.pert[(d.rep > 0.3) & (d.m > 0) & d.chronos.notna()])      # only targets that enter the summary go to survb
dirs = [x for x in dirs if x[0] in keep]
fd, fp, fb = (os.path.join(ST, f"{k}.tsv") for k in ("dirs", "pert", "bounds"))
pd.DataFrame(dirs, columns=["pert", "gene_index", "weight"]).to_csv(fd, sep="\t", index=False)
d[d.pert.isin(keep)][["pert", "n", "S", "chronos"]].to_csv(fp, sep="\t", index=False)
subprocess.run([SURVB, CA, fd, fp, fb], check=True)
b = pd.read_csv(fb, sep="\t")[["pert", "d_top", "mean_ctrl"]]
d = d.merge(b, on="pert", how="left")
d.to_csv(os.path.join(OUT, f"orion_o2_{LINE}.tsv"), sep="\t", index=False)

s = d[(d.rep > 0.3) & (d.m > 0) & d.chronos.notna()].copy()
s["share"] = (s.d_top / s.m).clip(upper=1)
bins = pd.cut(s.chronos, [-9, -2, -1.5, -1, -0.5, -0.3, 9])
print(f"reference cells (Chronos > -0.3): {ref:.0f}; targets with a reproducible signature: {len(s)} of {len(d)}")
print(f"median share {s.share.median():.3f}, share >= 0.5 in {(s.share >= 0.5).mean():.1%}, whole shift {(s.d_top >= s.m).mean():.1%}")
r, p = spearmanr(s.share, s.chronos); print(f"Spearman share vs Chronos: rho {r:.3f}, p {p:.2g}")
print(s.groupby(bins, observed=True).agg(targets=("share", "size"), S=("S", "median"), share=("share", "median"),
      half_or_more=("share", lambda x: (x >= 0.5).mean()), kd=("kd", "median")).round(3).to_string())
