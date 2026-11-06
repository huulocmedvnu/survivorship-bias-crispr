#!/usr/bin/env python3
"""Survivorship, phase 0 (laptop, committed fixtures only): how much cell loss is there in Tahoe-100M, and does
the measured effect magnitude grow as fewer cells survive?

Tahoe wells pool ~44 cell lines that receive one drug at one dose. Two survival proxies:
  * well survival  = cells in the well / median cells in a DMSO well (killing that hits all lines alike)
  * relative survival = a line's share of its well / that line's median share in DMSO wells (differential killing)

    python3 scripts/phase0_cell_loss.py
"""
import json, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr

FX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "inputs")
g = pd.DataFrame(json.load(open(os.path.join(FX, "tahoe_direct_curves_full.json")))["per_group"])
q = pd.read_csv(os.path.join(FX, "tahoe_quota_per_condition.csv")).rename(columns={"cvcl": "line", "N0": "N"})
g = g.merge(q[["drug", "line", "N", "dose_uM", "m_corr", "regime"]], on=["drug", "line", "N"], how="left")
is_dmso = g.drug.str.startswith("DMSO")
g["share"] = g.N / g.groupby("sample").N.transform("sum")
g["S_rel"] = g.share / g.line.map(g[is_dmso].groupby("line").share.median())
g["S_well"] = g.groupby("sample").N.transform("sum") / g[is_dmso].groupby("sample").N.sum().median()
t = g[~is_dmso & g.dose_uM.notna()]
print(f"treated conditions with dose: {len(t)}")
print("relative survival 5/50/95%:", t.S_rel.quantile([.05, .5, .95]).round(2).tolist())
print(f"conditions from wells with <50% of DMSO cells: {100*(t.S_well<.5).mean():.1f}%  (<25%: {100*(t.S_well<.25).mean():.1f}%)")
print(t.groupby(pd.cut(t.S_well, [0, .1, .25, .5, .75, 1.5, 10]), observed=True)
       .agg(n=("m_corr", "size"), median_m=("m_corr", "median")).round(2).to_string())
print(f"Spearman(well survival, m) = {spearmanr(t.S_well, t.m_corr).statistic:.2f}")
