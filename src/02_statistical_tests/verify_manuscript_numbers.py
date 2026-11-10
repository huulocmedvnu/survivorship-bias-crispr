"""Recompute every number quoted in the Abstract and Results from source_data/ and compare it with the manuscript.

Each line prints the value stated in the manuscript, the value recomputed here (rounded the same way) and
MATCH or MISMATCH. The report is also written to output/tables/manuscript_claims_verification.txt, and all
estimates to output/tables/numbers.json. Exit status 1 if any value differs.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd  # noqa: E402
from analysis import read, write_numbers  # noqa: E402
from paths import TABLES  # noqa: E402

N = write_numbers()
Q1, Q2, Q3, Q4 = "Q1 (most essential)", "Q2", "Q3", "Q4 (least essential)"
kp, ak, sh, dr = N["k562_pairs"], N["k562_authors_knockdown"], N["selection_share"], N["drepano"]
dl = {(r["method"].split(" ")[0], r["outcome"].split(" ")[0]): r for r in N["dose_labels"]}


def f(x, d=2):
    return f"{x:.{d}f}"


def pv(p):
    """P value as written in the manuscript: two significant digits above 0.001, else a.b x 10^e."""
    if p >= 0.001:
        return f"{p:.2g}" if p < 0.01 else f"{p:.2f}"
    e = math.floor(math.log10(p)); m = p / 10 ** e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    return f"{m:.1f}e{e}"


def ci(lo, hi, d=2):
    return f"{f(lo, d)} to {f(hi, d)}"


hepg2 = read("phaseR_hepg2_pert.tsv")   # S = n / reference, so the reference median is n / S where S < 1
gb = sh["K562 genome-wide, day 8"]["by_bin"]["< -2"]

# (result, figure, quantity, value stated in the manuscript, value recomputed from source_data)
CLAIMS = [
    ("R1", "Results", "identical guide pairs in both K562 libraries", "2,062 of 2,063",
     f"{N['k562_same_guides']['identical']:,} of {N['k562_same_guides']['total']:,}"),
    ("R1", "Fig. 2c", "perturbations with a reproducible day 6 response", "1,235", f"{kp['n']:,}"),
    ("R1", "Fig. 2c", "Q1 median Chronos", "-1.89", f(kp["chronos_median"][Q1])),
    ("R1", "Fig. 2c", "Q1 response ratio day 8 / day 6", "0.68", f(kp["R_by_quartile"][Q1][0])),
    ("R1", "Fig. 2c", "Q1 response ratio, 95% CI", "0.64 to 0.74", ci(*kp["R_by_quartile"][Q1][1:])),
    ("R1", "Fig. 2c", "Q2 response ratio", "0.85", f(kp["R_by_quartile"][Q2][0])),
    ("R1", "Fig. 2c", "Q3 response ratio", "0.91", f(kp["R_by_quartile"][Q3][0])),
    ("R1", "Fig. 2c", "Q4 response ratio", "0.91", f(kp["R_by_quartile"][Q4][0])),
    ("R1", "Fig. 2c", "Q4 median Chronos", "-0.38", f(kp["chronos_median"][Q4])),
    ("R1", "Fig. 2c", "Spearman rho, ratio vs Chronos", "0.25", f(kp["rho"])),
    ("R1", "Fig. 2c", "P, ratio vs Chronos", "1.0e-18", pv(kp["p"])),
    ("R1", "Fig. 2c", "adjusted change in log ratio per unit Chronos", "0.13", f(kp["adjusted"]["beta"])),
    ("R1", "Fig. 2c", "adjusted, 95% CI", "0.08 to 0.18", ci(*kp["adjusted"]["ci"])),
    ("R1", "Fig. 2c", "adjusted, P", "1.4e-7", pv(kp["adjusted"]["p"])),
    ("R1", "Fig. 2c", "placebo median ratio (day 6 half / half)", "1.00", f(kp["placebo"]["median"])),
    ("R1", "Fig. 2c", "placebo Spearman rho", "0.03", f(kp["placebo"]["rho"])),
    ("R1", "Fig. 2c", "placebo P", "0.31", pv(kp["placebo"]["p"])),
    ("R1", "Fig. 2a", "median cells, day 8, Q1 and Q4", "125 and 228", f"{kp['n8_median'][Q1]:.0f} and {kp['n8_median'][Q4]:.0f}"),
    ("R1", "Fig. 2a", "median cells, day 6, Q1 and Q4", "153 and 142", f"{kp['n6_median'][Q1]:.0f} and {kp['n6_median'][Q4]:.0f}"),
    ("R2", "Fig. 2b", "authors' knockdown Q1, day 6 to day 8 (%)", "76 to 63", f"{ak['d6'][Q1]:.0f} to {ak['d8'][Q1]:.0f}"),
    ("R2", "Fig. 2b", "authors' knockdown Q4, day 6 to day 8 (%)", "86 to 85", f"{ak['d6'][Q4]:.0f} to {ak['d8'][Q4]:.0f}"),
    ("R2", "Fig. 2b", "Spearman rho, change vs Chronos", "-0.30", f(ak["rho"])),
    ("R2", "Fig. 2b", "P, change vs Chronos", "7.8e-40", pv(ak["p"])),
    ("R2", "Fig. 2b", "genes", "1,900", f"{ak['n']:,}"),
    ("R3", "Fig. 3a", "essential guides, relative cells at day 4, CRISPRi", "0.80", f(dr["CRISPRi, RPL9 and DCAF13"]["rel_n"][4])),
    ("R3", "Fig. 3a", "essential guides, relative cells at day 4, CRISPRko", "0.68", f(dr["CRISPRko, RPL9 and DCAF13"]["rel_n"][4])),
    ("R3", "Fig. 3a", "essential guides, relative cells at day 7, CRISPRi", "0.00", f(dr["CRISPRi, RPL9 and DCAF13"]["rel_n"][7])),
    ("R3", "Fig. 3a", "essential guides, relative cells at day 14, CRISPRko", "0.00", f(dr["CRISPRko, RPL9 and DCAF13"]["rel_n"][14])),
    ("R3", "Fig. 3b", "CRISPRi knockdown, essential, day 7 (guides)", "-1.85 (3)",
     f"{f(dr['CRISPRi, RPL9 and DCAF13']['kd'][7])} ({dr['CRISPRi, RPL9 and DCAF13']['guides'][7]})"),
    ("R3", "Fig. 3b", "CRISPRi knockdown, essential, day 14 (guides)", "-0.66 (2)",
     f"{f(dr['CRISPRi, RPL9 and DCAF13']['kd'][14])} ({dr['CRISPRi, RPL9 and DCAF13']['guides'][14]})"),
    ("R3", "Fig. 3b", "CRISPRi knockdown, non-essential, days 7 to 14 (range)", "-2.08 to -1.93",
     ci(min(dr["CRISPRi, non-essential targets"]["kd"][d] for d in (7, 10, 14)),
        max(dr["CRISPRi, non-essential targets"]["kd"][d] for d in (7, 10, 14)))),
    ("R4", "Fig. S2", "relative cells, Chronos < -2, HCT116", "0.36", f(N["orion_cells"]["HCT116 genome-wide"]["< -2"])),
    ("R4", "Fig. S2", "relative cells, Chronos < -2, HEK293T", "0.22", f(N["orion_cells"]["HEK293T genome-wide"]["< -2"])),
    ("R4", "Fig. S2", "relative cells, Chronos < -2, K562 genome-wide", "0.65", f(N["orion_cells"]["K562 genome-wide, day 8"]["< -2"])),
]
for lab, short, val, cint in [("K562 genome-wide, day 8", "K562 genome-wide", "0.24", "0.19 to 0.28"),
                              ("RPE1 essential, day 7", "RPE1", "0.38", "0.32 to 0.46"),
                              ("Jurkat essential", "Jurkat", "0.13", "0.05 to 0.27"),
                              ("HCT116 genome-wide", "HCT116", "0.38", "0.31 to 0.47"),
                              ("HEK293T genome-wide", "HEK293T", "0.99", "0.73 to 1.00")]:
    b = sh[lab]["by_bin"]["< -2"]
    CLAIMS += [("R4", "Fig. 4", f"explainable fraction, Chronos < -2, {short}", val, f(b["median"])),
               ("R4", "Fig. 4", f"explainable fraction, 95% CI, {short}", cint, ci(*b["ci"]))]
CLAIMS += [
    ("R4", "Fig. 4", "explainable fraction, Chronos > -0.3, all seven screens", "0 in 7 of 7",
     f"0 in {sum(v['by_bin']['> -0.3']['median'] == 0 for v in sh.values())} of {len(sh)}"),
    ("R4", "Fig. 4", "screens with a negative gradient (rho < 0)", "5 of 7", f"{sum(v['rho'] < 0 for v in sh.values())} of {len(sh)}"),
    ("R4", "Fig. 4", "Spearman rho, RPE1", "-0.26", f(sh["RPE1 essential, day 7"]["rho"])),
    ("R4", "Fig. 4", "P, RPE1", "6.1e-31", pv(sh["RPE1 essential, day 7"]["p"])),
    ("R4", "Fig. 4", "Spearman rho, HEK293T", "-0.11", f(sh["HEK293T genome-wide"]["rho"])),
    ("R4", "Fig. 4", "P, HEK293T", "8.8e-14", pv(sh["HEK293T genome-wide"]["p"])),
    ("R4", "Results", "HepG2 median cells of non-essential perturbations", "45", f"{(hepg2.n / hepg2.S)[hepg2.S < 1].median():.0f}"),
]
for key, short, val, cint, p in [("joint", "joint EM", "0.11", "0.04 to 0.18", "0.002"),
                                 ("Song", "Song", "0.085", "0.02 to 0.15", "0.007"),
                                 ("Mixscale", "Mixscale", "0.081", "0.02 to 0.15", "0.015")]:
    r = dl[(key, "decrease")]
    d = 2 if key == "joint" else 3
    CLAIMS += [("R5", "Fig. 5a", f"buffered minus other, knockdown decrease, {short}", val, f(r["beta"], d)),
               ("R5", "Fig. 5a", f"95% CI, {short}", cint, ci(r["lo"], r["hi"])),
               ("R5", "Fig. 5a", f"P, {short}", p, f"{r['p']:.3f}")]
r = dl[("joint", "effect")]
CLAIMS += [
    ("R5", "Fig. 5b", "buffered minus other, response ratio, joint EM", "-0.11", f(r["beta"])),
    ("R5", "Fig. 5b", "95% CI, joint EM", "-0.18 to -0.03", ci(r["lo"], r["hi"])),
    ("R5", "Fig. 5b", "P, joint EM", "0.005", f"{r['p']:.3f}"),
    ("R6", "Fig. S1a", "Spearman rho, Tahoe share vs PRISM viability", "0.005", f(N["tahoe"]["rho_rel"], 3)),
    ("R6", "Fig. S1a", "P", "0.63", pv(N["tahoe"]["p_rel"])),
    ("R6", "Fig. S1a", "well and line pairs, drugs", "8,581, 126", f"{N['tahoe']['pairs']:,}, {N['tahoe']['drugs']}"),
    ("R6", "Fig. S1b", "cells per well vs plate control at 0.05, 0.5, 5 uM", "0.82, 0.74, 1.06",
     ", ".join(f(v) for v in N["tahoe"]["well_S_by_dose"].values())),
    ("R6", "Fig. S1b", "largest ratio between the control wells of one plate", "2.6", f(N["tahoe"]["dmso_well_ratio_max"], 1)),
    ("Abstract", "Abstract", "explainable fraction in five of seven screens, median range", "13% to 99%",
     f"{100 * sh['Jurkat essential']['by_bin']['< -2']['median']:.0f}% to {100 * sh['HEK293T genome-wide']['by_bin']['< -2']['median']:.0f}%"),
]

lines, bad = [], 0
w = max(len(c[2]) for c in CLAIMS)
for res, fig, what, stated, got in CLAIMS:
    ok = stated == got
    bad += not ok
    lines.append(f"[{res:<8}] {fig:<8} {what:<{w}}  manuscript {stated:<16} recomputed {got:<16} {'MATCH' if ok else 'MISMATCH'}")
lines.append(f"\n{len(CLAIMS) - bad} of {len(CLAIMS)} values match the manuscript.")
report = "\n".join(lines)
print(report)
TABLES.mkdir(parents=True, exist_ok=True)
(TABLES / "manuscript_claims_verification.txt").write_text(report + "\n")
sys.exit(1 if bad else 0)
