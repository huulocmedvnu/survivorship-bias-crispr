"""Summarise survb trimming bounds: share of the observed day-8 shift that selection alone could produce."""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr
O = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outputs/")
b = pd.read_csv(O + "phaseP_bounds.tsv", sep="\t")
p = pd.read_csv(O + "phaseP_pert.tsv", sep="\t")[["pert", "chronos"]]
k = pd.read_csv(O + "phaseP_killtest.tsv", sep="\t")[["pert", "R", "rep", "m6"]]
d = b.merge(p, on="pert").merge(k, on="pert")
d["se"] = d.sd_pert / np.sqrt(d.n_pert)
d = d[(d.rep > 0.3) & (d.m6 > 0) & (d.m > 2 * d.se) & d.chronos.notna()].copy()   # clear signal, same sign as day 6
d["share"] = (d.d_top / d.m).clip(upper=1)
d["explainable"] = d.d_top >= d.m
print(f"perturbations with a clear day-8 shift along the day-6 direction: {len(d)}")
print(f"S < 1 for {(d.S < 1).mean():.0%}, median share of shift selection could produce {d.share.median():.3f}")
print(f"selection alone could produce the whole shift: {d.explainable.sum()} ({d.explainable.mean():.1%})")
print(f"share >= 0.25: {(d.share >= 0.25).sum()}, share >= 0.5: {(d.share >= 0.5).sum()}")
r, pv = spearmanr(d.share, d.chronos); print(f"share vs Chronos: rho {r:.3f}, p {pv:.2g}")
r, pv = spearmanr(d.share, d.R); print(f"share vs day-8/day-6 ratio R: rho {r:.3f}, p {pv:.2g}")
d["q"] = pd.qcut(d.chronos, 4, labels=["Q1 most ess", "Q2", "Q3", "Q4 least"])
print(d.groupby("q", observed=True).agg(n=("share", "size"), S=("S", "median"), share=("share", "median"),
      share90=("share", lambda x: x.quantile(0.9)), explainable=("explainable", "mean"), R=("R", "median")).round(3).to_string())
print(d.sort_values("share", ascending=False)[["pert", "chronos", "S", "n_pert", "m", "d_top", "share", "R"]].head(15).round(3).to_string(index=False))
d.to_csv(O + "phaseP_bounds_summary.tsv", sep="\t", index=False)
