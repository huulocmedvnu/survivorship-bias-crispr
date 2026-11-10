"""Phase P: do day-8 (genome-wide) CRISPRi effects of
essential genes shrink relative to day 6 (essential library), as survivorship predicts?

Inputs (read-only): data/pseudobulk/{k562ess,gwps}.{f32,genes.txt,groups.tsv}
                    DepMap CRISPRGeneEffect.csv (Chronos), K562 = ACH-000551
Output: outputs/phaseP_killtest.tsv and a printed summary.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.formula.api as smf

PB = os.path.expanduser("data/pseudobulk/")
DEPMAP = os.path.expanduser("data/depmap/CRISPRGeneEffect.csv")
OUT = os.path.join(os.path.dirname(__file__), "..", "outputs")
EPS, MIN_CTRL, TOP, MIN_REP = 0.1, 0.5, 200, 0.3


def load(name):
    genes = open(PB + f"{name}.genes.txt").read().split()
    groups = pd.read_csv(PB + f"{name}.groups.tsv", sep="\t", header=None, names=["g", "n"])
    a = np.memmap(PB + f"{name}.f32", dtype=np.float32, mode="r").reshape(3, len(groups), len(genes))
    return genes, groups, a


g6, gr6, a6 = load("k562ess")
g8, gr8, a8 = load("gwps")
common = [g for g in g6 if g in set(g8)]
i6 = pd.Series(range(len(g6)), index=g6)[common].values
i8 = pd.Series(range(len(g8)), index=g8)[common].values
r6 = pd.Series(range(len(gr6)), index=gr6.g)
r8 = pd.Series(range(len(gr8)), index=gr8.g)
c6 = np.asarray(a6[0, r6["control"], i6], dtype=np.float64)
c8 = np.asarray(a8[0, r8["control"], i8], dtype=np.float64)
expr = (c6 >= MIN_CTRL) & (c8 >= MIN_CTRL)
genes = np.array(common)[expr]
c6, c8 = c6[expr], c8[expr]
i6, i8 = i6[expr], i8[expr]
gpos = pd.Series(range(len(genes)), index=genes)
print(f"genes expressed in both: {len(genes)}")

ge = pd.read_csv(DEPMAP, index_col=0)
ge.columns = [c.split(" (")[0] for c in ge.columns]
chronos = ge.loc["ACH-000551"]

lfc = lambda x, c: np.log2((x + EPS) / (c + EPS))
both = [p for p in gr6.g if p != "control" and p in r8.index]
rows = []
for p in both:
    x6 = [np.asarray(a6[k, r6[p], i6], dtype=np.float64) for k in range(3)]
    x8 = np.asarray(a8[0, r8[p], i8], dtype=np.float64)
    l6 = [lfc(x, c6) for x in x6]
    l8 = lfc(x8, c8)
    kd6 = l6[0][gpos[p]] if p in gpos.index else np.nan
    kd8 = l8[gpos[p]] if p in gpos.index else np.nan
    mask = np.ones(len(genes), bool)
    if p in gpos.index:
        mask[gpos[p]] = False
    res = {}
    for sel, mea in ((1, 2), (2, 1)):  # select genes in one day-6 half, measure in the other
        ls, lm = l6[sel], l6[mea]
        idx = np.where(mask)[0][np.argsort(-np.abs(ls[mask]))[:TOP]]
        u = ls[idx] / np.linalg.norm(ls[idx])
        res[sel] = dict(m6=lm[idx] @ u, m8=l8[idx] @ u, m6s=ls[idx] @ u,
                        rep=np.corrcoef(ls[idx], lm[idx])[0, 1])
    m6 = (res[1]["m6"] + res[2]["m6"]) / 2
    m8 = (res[1]["m8"] + res[2]["m8"]) / 2
    rows.append(dict(pert=p, n6=int(gr6.n[r6[p]]), n8=int(gr8.n[r8[p]]), kd6=kd6, kd8=kd8,
                     m6=m6, m8=m8, rep=(res[1]["rep"] + res[2]["rep"]) / 2,
                     placebo=res[1]["m6"] / res[2]["m6"] if res[2]["m6"] > 0 else np.nan,
                     chronos=chronos.get(p, np.nan)))
d = pd.DataFrame(rows)
d["R"] = d.m8 / d.m6
d["logR"] = np.log(d.R.where(d.R > 0))
d["dep8"] = np.log2(d.n8 / d.n8.median())
d["dep6"] = np.log2(d.n6 / d.n6.median())
os.makedirs(OUT, exist_ok=True)
d.to_csv(os.path.join(OUT, "phaseP_killtest.tsv"), sep="\t", index=False)

k = d[(d.rep > MIN_REP) & (d.m6 > 0) & d.chronos.notna()].copy()
print(f"perturbations in both screens: {len(d)}, reproducible day-6 signature with Chronos: {len(k)}")
print(f"median R = {k.R.median():.2f} (global day-8 / day-6 level, absorbs batch)")
for var in ["chronos", "dep8", "dep6"]:
    r, pv = spearmanr(k.R, k[var])
    print(f"Spearman R vs {var}: rho = {r:.3f}, p = {pv:.2g}")
# placebo: two held-out day-6 measurements (genes chosen in A measured in B, over chosen in B measured in A)
pl = k[np.isfinite(k.placebo) & (k.placebo > 0)]
print(f"placebo median = {pl.placebo.median():.2f}")
r, pv = spearmanr(pl.placebo, pl.chronos)
print(f"placebo (day-6 half vs half) vs chronos: rho = {r:.3f}, p = {pv:.2g}")
# adjusted model
kk = k[np.isfinite(k.logR) & k.kd6.notna()]
fit = smf.ols("logR ~ chronos + np.log(m6) + kd6", data=kk).fit()
print(f"OLS logR ~ chronos + log m6 + kd6 (n={len(kk)}): chronos beta = {fit.params['chronos']:.3f}, "
      f"p = {fit.pvalues['chronos']:.2g}")
# knockdown at day 8 vs day 6
kd = k.dropna(subset=["kd6", "kd8"])
r, pv = spearmanr(kd.kd8 - kd.kd6, kd.chronos)
print(f"knockdown change (kd8 - kd6, positive = weaker at day 8) vs chronos: rho = {r:.3f}, p = {pv:.2g}")
# by essentiality quartile
k["q"] = pd.qcut(k.chronos, 4, labels=["Q1 most essential", "Q2", "Q3", "Q4 least"])
print(k.groupby("q", observed=True).agg(n=("R", "size"), chronos=("chronos", "median"), R=("R", "median"),
                                        kd6=("kd6", "median"), kd8=("kd8", "median"),
                                        n8=("n8", "median"), n6=("n6", "median")).round(2).to_string())
