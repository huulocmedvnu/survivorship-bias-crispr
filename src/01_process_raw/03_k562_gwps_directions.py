"""Phase P step 2 input: one direction u per perturbation, taken from the day-6 essential screen (independent cells),
expressed on the day-8 genome-wide cache genes, plus the surviving fraction S at day 8.

u: top 200 genes by |lfc| of the day-6 all-cell pseudobulk against day-6 controls (genes with control mean >= 0.5
CP10k in both screens, target gene excluded), unit norm, signed.
S: day-8 cell count of the perturbation over the median day-8 count of genome-wide perturbations with Chronos
> -0.2 (non-essential in K562), capped at 1.
Output: data/intermediate/phaseP_dirs.tsv (pert, gene_index, weight) and source_data/phaseP_pert.tsv (pert, n8, S, chronos).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import DATA, SOURCE, SCORER  # noqa: E402
import os
import numpy as np
import pandas as pd

PB = f"{DATA}/pseudobulk/"
CACHE = f"{DATA}/dosecurve_cache/gwps/"
DEPMAP = f"{DATA}/depmap/CRISPRGeneEffect.csv"
OUT = str(SOURCE); INTER = f"{DATA}/intermediate"; os.makedirs(INTER, exist_ok=True)
EPS, MIN_CTRL, TOP = 0.1, 0.5, 200

g6 = open(PB + "k562ess.genes.txt").read().split()
gr6 = pd.read_csv(PB + "k562ess.groups.tsv", sep="\t", header=None, names=["g", "n"])
a6 = np.memmap(PB + "k562ess.f32", dtype=np.float32, mode="r").reshape(3, len(gr6), len(g6))
g8 = open(CACHE + "genes.txt").read().split()
gr8 = pd.read_csv(CACHE + "groups.tsv", sep="\t", header=None, names=["g", "off", "n"])
pb8 = np.memmap(PB + "gwps.f32", dtype=np.float32, mode="r").reshape(3, len(gr8), len(g8))
pos8 = pd.Series(range(len(g8)), index=g8)
pos6 = pd.Series(range(len(g6)), index=g6)
common = [g for g in g6 if g in pos8.index]
c6 = np.asarray(a6[0, 0, pos6[common].values], dtype=np.float64)
c8 = np.asarray(pb8[0, 0, pos8[common].values], dtype=np.float64)
assert gr6.g[0] == "control" and gr8.g[0] == "control"
keep = (c6 >= MIN_CTRL) & (c8 >= MIN_CTRL)
genes = np.array(common)[keep]
c6 = c6[keep]
i6, i8 = pos6[genes].values, pos8[genes].values

ge = pd.read_csv(DEPMAP, index_col=0)
ge.columns = [c.split(" (")[0] for c in ge.columns]
chronos = ge.loc["ACH-000551"]
n8 = gr8.set_index("g").n
ref = n8[[g for g in n8.index if g != "control" and chronos.get(g, np.nan) > -0.2]].median()
print(f"reference day-8 cells per non-essential perturbation: {ref:.0f}")

dirs, perts = [], []
r8 = set(gr8.g)
for r, p in enumerate(gr6.g):
    if p == "control" or p not in r8:
        continue
    l = np.log2((np.asarray(a6[0, r, i6], dtype=np.float64) + EPS) / (c6 + EPS))
    ok = genes != p
    idx = np.where(ok)[0][np.argsort(-np.abs(l[ok]))[:TOP]]
    u = l[idx] / np.linalg.norm(l[idx])
    dirs += [(p, int(i8[j]), float(w)) for j, w in zip(idx, u)]
    perts.append((p, int(n8[p]), min(1.0, n8[p] / ref), chronos.get(p, np.nan)))
pd.DataFrame(dirs, columns=["pert", "gene_index", "weight"]).to_csv(os.path.join(INTER, "phaseP_dirs.tsv"), sep="\t", index=False)
pp = pd.DataFrame(perts, columns=["pert", "n8", "S", "chronos"])
pp.to_csv(os.path.join(OUT, "phaseP_pert.tsv"), sep="\t", index=False)
print(f"{len(pp)} perturbations, {len(dirs)} direction entries, median S {pp.S.median():.2f}, "
      f"S < 0.5 in {(pp.S < 0.5).mean():.0%}")
