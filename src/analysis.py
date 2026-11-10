"""Every estimate reported in the manuscript, computed from source_data/ only.

compute_all() returns the reported numbers (key "N", also written as output/tables/numbers.json) and the
summaries that the figure scripts plot. Bootstrap confidence intervals use 1,000 resamples from a single
generator, numpy.random.default_rng(1), drawn in a fixed order, so every run returns identical values.
"""
import json
import warnings
from functools import lru_cache

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr

from paths import SOURCE, TABLES

warnings.filterwarnings("ignore", category=FutureWarning)

QL = ["Q1 (most essential)", "Q2", "Q3", "Q4 (least essential)"]
BINS = [-9, -2, -1.5, -1, -0.5, -0.3, 9]
BL = ["< -2", "-2 to -1.5", "-1.5 to -1", "-1 to -0.5", "-0.5 to -0.3", "> -0.3"]


def read(name, **kw):
    return pd.read_csv(SOURCE / name, sep="\t", **kw)


@lru_cache(maxsize=1)
def compute_all():
    rng = np.random.default_rng(1)

    def boot_median(x, b=1000):
        x = np.asarray(x, float); x = x[np.isfinite(x)]
        if len(x) < 5:
            return np.nan, np.nan, np.nan
        m = np.median(rng.choice(x, (b, len(x))), axis=1)
        return np.median(x), np.quantile(m, 0.025), np.quantile(m, 0.975)

    N, P = {}, {}
    K562 = read("ref_k562_chronos.tsv", index_col=0).chronos

    # ---- Fig. 2: K562 day 6 against day 8
    k = read("phaseP_killtest.tsv")
    k = k[(k.rep > 0.3) & (k.m6 > 0) & k.chronos.notna()].copy()
    k["q"] = pd.qcut(k.chronos, 4, labels=QL)
    r, p = spearmanr(k.R, k.chronos)
    N["k562_pairs"] = dict(n=len(k), rho=round(r, 3), p=float(p),
                           R_by_quartile={q: [round(v, 3) for v in boot_median(k.R[k.q == q])] for q in QL},
                           chronos_median={q: round(k.chronos[k.q == q].median(), 2) for q in QL},
                           n8_median={q: float(k.n8[k.q == q].median()) for q in QL},
                           n6_median={q: float(k.n6[k.q == q].median()) for q in QL})
    kk = k[np.isfinite(k.logR) & k.kd6.notna()]
    f = smf.ols("logR ~ chronos + np.log(m6) + kd6", data=kk).fit()
    N["k562_pairs"]["adjusted"] = dict(n=int(f.nobs), beta=round(f.params["chronos"], 3),
                                       ci=[round(v, 3) for v in f.conf_int().loc["chronos"]], p=float(f.pvalues["chronos"]))
    pl = k[np.isfinite(k.placebo) & (k.placebo > 0)]
    r2, p2 = spearmanr(pl.placebo, pl.chronos)
    N["k562_pairs"]["placebo"] = dict(median=round(pl.placebo.median(), 3), rho=round(r2, 3), p=float(p2))
    kd = read("ref_replogle_knockdown.tsv", index_col=0)
    kd["chronos"] = kd.index.map(K562); kd = kd.dropna()
    kd["q"] = pd.qcut(kd.chronos, 4, labels=QL)
    r3, p3 = spearmanr(kd.d8 - kd.d6, kd.chronos)
    N["k562_authors_knockdown"] = dict(n=len(kd), rho=round(r3, 3), p=float(p3),
                                       d6={q: round(-100 * kd.d6[kd.q == q].median(), 1) for q in QL},
                                       d8={q: round(-100 * kd.d8[kd.q == q].median(), 1) for q in QL})
    lib = read("phaseP_same_guides.tsv", index_col=0)
    N["k562_same_guides"] = dict(identical=int(lib.iloc[:, 0].sum()), total=len(lib))
    P["fig2"] = dict(
        cells={d: [k[d][k.q == q].median() for q in QL] for d in ("n6", "n8")},
        knockdown={d: [-100 * kd[d][kd.q == q].median() for q in QL] for d in ("d6", "d8")},
        ratio=[boot_median(k.R[k.q == q]) for q in QL],
        placebo=[boot_median(pl.placebo[pl.q == q]) for q in QL])

    # ---- Fig. 3: Drepanos timecourse
    dq = read("phaseQ_drepano.tsv")
    N["drepano"], P["fig3"] = {}, []
    for mod, ess in [("CRISPRi", True), ("CRISPRi", False), ("CRISPRko", True), ("CRISPRko", False)]:
        s = dq[(dq["mod"] == mod) & (dq.essential == ess)]
        lab = f"{mod}, {'RPL9 and DCAF13' if ess else 'non-essential targets'}"
        rel = s.groupby("day").rel_n.median()
        sk = s[s.n >= 10].groupby("day").agg(kd=("kd", "median"), guides=("guide", "nunique"))
        N["drepano"][lab] = dict(rel_n={int(d): round(v, 2) for d, v in rel.items()},
                                 kd={int(d): round(v, 2) for d, v in sk.kd.items()},
                                 guides={int(d): int(v) for d, v in sk.guides.items()})
        P["fig3"].append(dict(mod=mod, essential=ess, label=lab, rel=rel, kd=sk))

    # ---- Fig. 4: fraction explainable by selective survival, by Chronos bin
    screens = {}
    s = read("phaseP_bounds_summary.tsv"); screens["K562 genome-wide, day 8"] = s[["chronos", "share"]]
    for line, lab in [("k562ess", "K562 essential, day 6"), ("rpe1", "RPE1 essential, day 7"),
                      ("jurkat", "Jurkat essential"), ("hepg2", "HepG2 essential")]:
        b = read(f"phaseR_{line}_bounds.tsv").merge(read(f"phaseR_{line}_pert.tsv")[["pert", "chronos"]], on="pert")
        b = b[(b.m > 2 * b.sd_pert / np.sqrt(b.n_pert)) & b.chronos.notna()]
        b["share"] = (b.d_top / b.m).clip(upper=1); screens[lab] = b[["chronos", "share"]]
    for line in ["HCT116", "HEK293T"]:
        b = read(f"orion_o2_{line}.tsv"); b = b[(b.rep > 0.3) & (b.m > 0) & b.chronos.notna() & b.d_top.notna()]
        b["share"] = (b.d_top / b.m).clip(upper=1); screens[f"{line} genome-wide"] = b[["chronos", "share"]]
    N["selection_share"], P["fig4"] = {}, {}
    for lab, d in screens.items():
        d = d.assign(bin=pd.cut(d.chronos, BINS, labels=BL))
        st = [boot_median(d.share[d.bin == bl]) + (int((d.bin == bl).sum()),) for bl in BL]
        r, p = spearmanr(d.share, d.chronos)
        N["selection_share"][lab] = dict(
            n=len(d), rho=round(r, 3), p=float(p), median_overall=round(d.share.median(), 3),
            by_bin={bl: dict(median=round(v[0], 3) if np.isfinite(v[0]) else None,
                             ci=[round(v[1], 3), round(v[2], 3)] if np.isfinite(v[1]) else None, n=v[3])
                    for bl, v in zip(BL, st)})
        P["fig4"][lab] = st

    # ---- Fig. 5: dose-response labels
    dl = read("phaseP_dose_labels.tsv")
    dl = dl[dl.m6 > 0].copy(); dl["lm6"] = np.log(dl.m6); dl["lnc"] = np.log(dl.n_cells); dl["lnde"] = np.log1p(dl.n_de_lfsr)
    rows = []
    for meth, lab in [("em", "joint EM (dosecurve)"), ("song", "Song et al. 2025"), ("mix", "Mixscale")]:
        x = dl.copy(); x["buf"] = (x[meth] == "buffered").astype(int)
        if meth == "em":
            x = x[x.em.isin(["buffered", "indeterminate"])]
        for y, ylab in [("dkd", "decrease in knockdown, day 6 to day 8"), ("R", "effect ratio, day 8 / day 6")]:
            z = x.dropna(subset=[y, "chronos", "kd6"])
            fm = smf.ols(f"{y} ~ buf + chronos + lm6 + lnc + lnde + kd6", data=z).fit()
            ci = fm.conf_int().loc["buf"]
            rows.append(dict(method=lab, outcome=ylab, beta=fm.params["buf"], lo=ci.iloc[0], hi=ci.iloc[1],
                             p=fm.pvalues["buf"], n=int(fm.nobs), n_buf=int(z.buf.sum())))
    t5 = pd.DataFrame(rows); N["dose_labels"] = t5.round(4).to_dict("records"); P["fig5"] = t5

    # ---- Fig. S1: Tahoe-100M cell numbers against PRISM viability
    tp = read("tahoe_t1_prism.tsv.gz")
    tp["vb"] = pd.cut(tp.v, [-0.01, 0.25, 0.5, 0.75, 0.95, 1.3], labels=["< 0.25", "0.25 to 0.5", "0.5 to 0.75", "0.75 to 0.95", "> 0.95"])
    ts = read("tahoe_t1_survival.tsv.gz")
    w = ts.groupby(["sample", "plate", "drug", "conc"], as_index=False).n.sum()
    dm = w[w.drug.str.contains("DMSO")].groupby("plate").n.mean().rename("dm"); w = w.join(dm, on="plate"); w["S"] = w.n / w.dm
    wt = w[~w.drug.str.contains("DMSO")]
    r1, p1 = spearmanr(tp.rel, tp.v); r2, p2 = spearmanr(tp.S, tp.v)
    N["tahoe"] = dict(pairs=len(tp), drugs=int(tp.key.nunique()), rho_rel=round(r1, 3), p_rel=float(p1), rho_S=round(r2, 3), p_S=float(p2),
                      well_S_by_dose={str(c): round(wt.S[wt.conc == c].median(), 2) for c in sorted(wt.conc.dropna().unique())},
                      dmso_well_ratio_max=round(w[w.drug.str.contains("DMSO")].groupby("plate").n.agg(lambda v: v.max() / v.min()).max(), 2))
    P["figS1"] = dict(prism=tp, wells=wt)

    # ---- Fig. S2: cells per perturbation by essentiality in genome-wide screens
    N["orion_cells"], P["figS2"] = {}, {}
    g8 = read("ref_gwps_cells.tsv")
    g8["chronos"] = g8.g.map(K562); g8 = g8.dropna(); g8["rel"] = g8.n / g8.n[g8.chronos > -0.3].median()
    sets = {"K562 genome-wide, day 8": g8}
    for line in ["HCT116", "HEK293T"]:
        sets[f"{line} genome-wide"] = read(f"orion_o1_{line}.tsv")
    for lab, d in sets.items():
        b = pd.cut(d.chronos, BINS, labels=BL); m = d.groupby(b, observed=True).rel.median()
        N["orion_cells"][lab] = {k_: round(v, 2) for k_, v in m.items()}
        P["figS2"][lab] = m
    return N, P


def write_numbers():
    N, _ = compute_all()
    TABLES.mkdir(parents=True, exist_ok=True)
    json.dump(N, open(TABLES / "numbers.json", "w"), indent=1, default=float)
    return N
