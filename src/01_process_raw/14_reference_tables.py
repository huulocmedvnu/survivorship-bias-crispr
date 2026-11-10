"""Write the three small reference tables that figures/make_figures.py reads, from the public raw files.

Inputs (see README, full reproduction):
  data/depmap/CRISPRGeneEffect.csv               DepMap Public 24Q4 Chronos gene effects
  data/replogle_supp/mmc1.xlsx, mmc2.xlsx         Tables S1 and S2 of Replogle et al. 2022 (per-perturbation summary statistics)
  data/dosecurve_cache/gwps/groups.tsv           cells per perturbation in the K562 genome-wide screen
Outputs:
  source_data/ref_k562_chronos.tsv                   gene, chronos (K562, ACH-000551)
  source_data/ref_replogle_knockdown.tsv             gene, d8, d6 (mean percent knockdown reported by the authors)
  source_data/ref_gwps_cells.tsv                     g, n
  source_data/phaseP_same_guides.tsv                 gene, same_guides (identical sgRNA pairs in both K562 libraries)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import DATA, SOURCE, SCORER  # noqa: E402
import os
import pandas as pd

D, O = f"{DATA}/", f"{SOURCE}/"

ge = pd.read_csv(D + "depmap/CRISPRGeneEffect.csv", index_col=0)
ge.columns = [c.split(" (")[0] for c in ge.columns]
ge.loc["ACH-000551"].dropna().rename("chronos").rename_axis("gene").to_csv(O + "ref_k562_chronos.tsv", sep="\t")

x = pd.ExcelFile(D + "replogle_supp/mmc2.xlsx")
a8, a6 = x.parse("TabA_K562_day8_summary_stat"), x.parse("TabB_K562_day6_summary_stat")
g = lambda s: s.astype(str).str.split("_").str[1]
a8["gene"], a6["gene"] = g(a8["genetic perturbation"]), g(a6["genetic perturbation"])
kd = a8.groupby("gene")["percent knockdown"].mean().to_frame("d8").join(
    a6.groupby("gene")["percent knockdown"].mean().rename("d6"), how="inner")
kd.to_csv(O + "ref_replogle_knockdown.tsv", sep="\t")

pd.read_csv(D + "dosecurve_cache/gwps/groups.tsv", sep="\t", header=None, names=["g", "o", "n"])[["g", "n"]].to_csv(
    O + "ref_gwps_cells.tsv", sep="\t", index=False)

# identical guide pairs in the two K562 libraries (Table S1 of Replogle et al. 2022)
x = pd.ExcelFile(D + "replogle_supp/mmc1.xlsx")
pairs = lambda t: t.groupby("gene").apply(
    lambda g: frozenset(zip(g["targeting sequence A"], g["targeting sequence B"])), include_groups=False)
l8, l6 = pairs(x.parse("TabA_K562_day8_library")), pairs(x.parse("TabB_K562_day6_library"))
both = l8.index.intersection(l6.index)
pd.Series({g: l8[g] == l6[g] for g in both}, name="same_guides").to_csv(O + "phaseP_same_guides.tsv", sep="\t")
