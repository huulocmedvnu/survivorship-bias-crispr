"""Numbers and figures for the manuscript, from the committed outputs only.
Writes figures/numbers.json (single source for every number quoted in the text) and figures/*.{png,pdf}.
Style: plotly_white, default Plotly colours, Arial, black text, no gridlines, lowercase panel letters."""
import os, json
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import spearmanr

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
O = os.path.join(R, "outputs/"); F = os.path.join(R, "figures/"); os.makedirs(F, exist_ok=True)
DEP = os.path.join(R, "data/depmap/")
SUPP = os.path.join(R, "data/replogle_supp/")
C = ["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3", "#FF6692"]
BLACK = "#111111"
rng = np.random.default_rng(1)
N = {}

ge = pd.read_csv(DEP + "CRISPRGeneEffect.csv", index_col=0); ge.columns = [c.split(" (")[0] for c in ge.columns]
K562 = ge.loc["ACH-000551"]


def boot_median(x, b=1000):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 5:
        return np.nan, np.nan, np.nan
    m = np.median(rng.choice(x, (b, len(x))), axis=1)
    return np.median(x), np.quantile(m, 0.025), np.quantile(m, 0.975)


def style(fig, w, h):
    fig.update_layout(template="plotly_white", width=w, height=h, font=dict(family="Arial", size=13, color=BLACK),
                      margin=dict(l=60, r=20, t=40, b=60), legend=dict(font=dict(color=BLACK)))
    fig.update_xaxes(showgrid=False, showline=True, linecolor=BLACK, ticks="outside", tickfont=dict(color=BLACK),
                     title_font=dict(color=BLACK))
    fig.update_yaxes(showgrid=False, showline=True, linecolor=BLACK, ticks="outside", tickfont=dict(color=BLACK),
                     title_font=dict(color=BLACK))
    for a in fig.layout.annotations:
        a.font = dict(family="Arial", size=14, color=BLACK)
    return fig


def save(fig, name):
    fig.write_image(F + name + ".png", scale=4); fig.write_image(F + name + ".pdf")


QL = ["Q1 (most essential)", "Q2", "Q3", "Q4 (least essential)"]

# ---------------------------------------------------------------- Figure 1: K562 day 6 vs day 8
k = pd.read_csv(O + "phaseP_killtest.tsv", sep="\t")
k = k[(k.rep > 0.3) & (k.m6 > 0) & k.chronos.notna()].copy()
k["q"] = pd.qcut(k.chronos, 4, labels=QL)
r, p = spearmanr(k.R, k.chronos)
N["k562_pairs"] = dict(n=len(k), rho=round(r, 3), p=float(p),
                       R_by_quartile={q: [round(v, 3) for v in boot_median(k.R[k.q == q])] for q in QL},
                       chronos_median={q: round(k.chronos[k.q == q].median(), 2) for q in QL},
                       n8_median={q: float(k.n8[k.q == q].median()) for q in QL}, n6_median={q: float(k.n6[k.q == q].median()) for q in QL})
kk = k[np.isfinite(k.logR) & k.kd6.notna()]
f = smf.ols("logR ~ chronos + np.log(m6) + kd6", data=kk).fit()
N["k562_pairs"]["adjusted"] = dict(n=int(f.nobs), beta=round(f.params["chronos"], 3), ci=[round(v, 3) for v in f.conf_int().loc["chronos"]], p=float(f.pvalues["chronos"]))
pl = k[np.isfinite(k.placebo) & (k.placebo > 0)]
r2, p2 = spearmanr(pl.placebo, pl.chronos)
N["k562_pairs"]["placebo"] = dict(median=round(pl.placebo.median(), 3), rho=round(r2, 3), p=float(p2))
# authors' knockdown (Table S2)
x = pd.ExcelFile(SUPP + "mmc2.xlsx")
a8 = x.parse("TabA_K562_day8_summary_stat"); a6 = x.parse("TabB_K562_day6_summary_stat")
g = lambda s: s.astype(str).str.split("_").str[1]
a8["gene"], a6["gene"] = g(a8["genetic perturbation"]), g(a6["genetic perturbation"])
kd = a8.groupby("gene")["percent knockdown"].mean().to_frame("d8").join(a6.groupby("gene")["percent knockdown"].mean().rename("d6"), how="inner")
kd["chronos"] = kd.index.map(K562); kd = kd.dropna()
kd["q"] = pd.qcut(kd.chronos, 4, labels=QL)
r3, p3 = spearmanr(kd.d8 - kd.d6, kd.chronos)
N["k562_authors_knockdown"] = dict(n=len(kd), rho=round(r3, 3), p=float(p3),
                                   d6={q: round(-100 * kd.d6[kd.q == q].median(), 1) for q in QL},
                                   d8={q: round(-100 * kd.d8[kd.q == q].median(), 1) for q in QL})
lib = pd.read_csv(O + "phaseP_same_guides.tsv", sep="\t", index_col=0)
N["k562_same_guides"] = dict(identical=int(lib.iloc[:, 0].sum()), total=len(lib))

fig = make_subplots(rows=1, cols=3, subplot_titles=("a  Cells per perturbation", "b  On-target knockdown (authors)", "c  Downstream response, day 8 / day 6"),
                    horizontal_spacing=0.09)
for j, (day, col) in enumerate([("n6", C[0]), ("n8", C[1])]):
    fig.add_bar(x=QL, y=[k[day][k.q == q].median() for q in QL], name="day 6" if j == 0 else "day 8", marker_color=col, row=1, col=1)
for j, (day, col) in enumerate([("d6", C[0]), ("d8", C[1])]):
    fig.add_bar(x=QL, y=[-100 * kd[day][kd.q == q].median() for q in QL], marker_color=col, showlegend=False, row=1, col=2)
med = [boot_median(k.R[k.q == q]) for q in QL]
fig.add_scatter(x=QL, y=[m[0] for m in med], mode="markers", marker=dict(symbol="diamond", size=11, color=C[1]), name="day 8 / day 6",
                error_y=dict(type="data", symmetric=False, array=[m[2] - m[0] for m in med], arrayminus=[m[0] - m[1] for m in med], color=C[1], thickness=1.2),
                row=1, col=3)
pm = [boot_median(pl.placebo[pl.q == q]) for q in QL]
fig.add_scatter(x=QL, y=[m[0] for m in pm], mode="markers", marker=dict(symbol="circle", size=9, color="#7F7F7F"), name="day 6 half / half",
                error_y=dict(type="data", symmetric=False, array=[m[2] - m[0] for m in pm], arrayminus=[m[0] - m[1] for m in pm], color="#7F7F7F", thickness=1.2),
                row=1, col=3)
fig.add_hline(y=1, line=dict(color="#7F7F7F", dash="dot", width=1), row=1, col=3)
fig.update_yaxes(title_text="median cells", row=1, col=1); fig.update_yaxes(title_text="median on-target knockdown (%)", range=[0, 100], row=1, col=2)
fig.update_yaxes(title_text="median ratio", range=[0, 1.3], row=1, col=3)
fig.update_xaxes(title_text="Chronos quartile in K562", tickangle=-30)
fig.update_layout(barmode="group")
save(style(fig, 1300, 460), "Figure2")

# ---------------------------------------------------------------- Figure 2: Drepanos timecourse
dq = pd.read_csv(O + "phaseQ_drepano.tsv", sep="\t")
fig = make_subplots(rows=1, cols=3, subplot_titles=("a  Cells per guide, relative", "b  On-target knockdown, CRISPR interference", "c  On-target knockdown, CRISPR knockout"),
                    horizontal_spacing=0.09)
N["drepano"] = {}
for mod, ess, col, dash in [("CRISPRi", True, C[1], "solid"), ("CRISPRi", False, C[0], "solid"), ("CRISPRko", True, C[1], "dash"), ("CRISPRko", False, C[0], "dash")]:
    s = dq[(dq["mod"] == mod) & (dq.essential == ess)]
    lab = f"{mod}, {'RPL9 and DCAF13' if ess else 'non-essential targets'}"
    rel = s.groupby("day").rel_n.median()
    sk = s[s.n >= 10].groupby("day").agg(kd=("kd", "median"), guides=("guide", "nunique"))
    fig.add_scatter(x=rel.index, y=rel.values, mode="lines+markers", line=dict(color=col, dash=dash), name=lab, row=1, col=1)
    fig.add_scatter(x=sk.index, y=sk.kd.values, mode="lines+markers+text", text=[str(v) for v in sk.guides],
                    textposition="top center" if ess else "bottom center", textfont=dict(size=11, color=BLACK),
                    line=dict(color=col, dash=dash), showlegend=False, row=1, col=2 if mod == "CRISPRi" else 3)
    N["drepano"][lab] = dict(rel_n={int(d): round(v, 2) for d, v in rel.items()}, kd={int(d): round(v, 2) for d, v in sk.kd.items()}, guides={int(d): int(v) for d, v in sk.guides.items()})
fig.update_xaxes(title_text="days after transduction", tickvals=[4, 7, 10, 14], range=[3, 15])
fig.update_yaxes(title_text="cells relative to non-essential guides", range=[-0.05, 1.1], row=1, col=1)
fig.update_yaxes(title_text="target log2 fold change", range=[-2.4, -0.3], row=1, col=2)
fig.update_yaxes(range=[-1.3, -0.3], row=1, col=3)
fig.update_layout(legend=dict(orientation="h", y=-0.28, x=0))
save(style(fig, 1300, 500), "Figure3")

# ---------------------------------------------------------------- Figure 3: selection share by Chronos bin across screens
BINS = [-9, -2, -1.5, -1, -0.5, -0.3, 9]; BL = ["< -2", "-2 to -1.5", "-1.5 to -1", "-1 to -0.5", "-0.5 to -0.3", "> -0.3"]
screens = {}
s = pd.read_csv(O + "phaseP_bounds_summary.tsv", sep="\t"); screens["K562 genome-wide, day 8"] = s[["chronos", "share"]]
for line, lab in [("k562ess", "K562 essential, day 6"), ("rpe1", "RPE1 essential, day 7"), ("jurkat", "Jurkat essential"), ("hepg2", "HepG2 essential")]:
    b = pd.read_csv(O + f"phaseR_{line}_bounds.tsv", sep="\t").merge(pd.read_csv(O + f"phaseR_{line}_pert.tsv", sep="\t")[["pert", "chronos"]], on="pert")
    b = b[(b.m > 2 * b.sd_pert / np.sqrt(b.n_pert)) & b.chronos.notna()]
    b["share"] = (b.d_top / b.m).clip(upper=1); screens[lab] = b[["chronos", "share"]]
for line in ["HCT116", "HEK293T"]:
    b = pd.read_csv(O + f"orion_o2_{line}.tsv", sep="\t"); b = b[(b.rep > 0.3) & (b.m > 0) & b.chronos.notna() & b.d_top.notna()]
    b["share"] = (b.d_top / b.m).clip(upper=1); screens[f"{line} genome-wide"] = b[["chronos", "share"]]
fig = go.Figure(); N["selection_share"] = {}
for i, (lab, d) in enumerate(screens.items()):
    d = d.assign(bin=pd.cut(d.chronos, BINS, labels=BL))
    st = [boot_median(d.share[d.bin == bl]) + (int((d.bin == bl).sum()),) for bl in BL]
    r, p = spearmanr(d.share, d.chronos)
    N["selection_share"][lab] = dict(n=len(d), rho=round(r, 3), p=float(p), median_overall=round(d.share.median(), 3),
                                     by_bin={bl: dict(median=round(v[0], 3) if np.isfinite(v[0]) else None, ci=[round(v[1], 3), round(v[2], 3)] if np.isfinite(v[1]) else None, n=v[3]) for bl, v in zip(BL, st)})
    ok = [j for j, v in enumerate(st) if v[3] >= 15]
    fig.add_scatter(x=[BL[j] for j in ok], y=[st[j][0] for j in ok], mode="lines+markers", name=lab, line=dict(color=C[i % 7], width=1.5),
                    marker=dict(symbol="diamond", size=9),
                    error_y=dict(type="data", symmetric=False, array=[st[j][2] - st[j][0] for j in ok], arrayminus=[st[j][0] - st[j][1] for j in ok], thickness=1, width=3))
fig.update_xaxes(title_text="Chronos gene effect of the target in the screened line (bins with at least 15 targets)")
fig.update_yaxes(title_text="median fraction explainable by selective survival", range=[-0.05, 1.05])
save(style(fig, 1000, 520), "Figure4")

# ---------------------------------------------------------------- Figure 4: dose-response labels and knockdown loss
dl = pd.read_csv(O + "phaseP_dose_labels.tsv", sep="\t")
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
        rows.append(dict(method=lab, outcome=ylab, beta=fm.params["buf"], lo=ci[0], hi=ci[1], p=fm.pvalues["buf"], n=int(fm.nobs), n_buf=int(z.buf.sum())))
t4 = pd.DataFrame(rows); N["dose_labels"] = t4.round(4).to_dict("records")
fig = make_subplots(rows=1, cols=2, subplot_titles=("a  Decrease in on-target knockdown, day 6 to day 8", "b  Downstream response ratio, day 8 / day 6"), horizontal_spacing=0.25)
for j, oc in enumerate(t4.outcome.unique()):
    s = t4[t4.outcome == oc]
    fig.add_scatter(x=s.beta, y=s.method, mode="markers", marker=dict(symbol="diamond", size=11, color=C[j]),
                    error_x=dict(type="data", symmetric=False, array=s.hi - s.beta, arrayminus=s.beta - s.lo, color=C[j], thickness=1.2), showlegend=False, row=1, col=j + 1)
    fig.add_vline(x=0, line=dict(color="#7F7F7F", dash="dot", width=1), row=1, col=j + 1)
fig.update_xaxes(title_text="buffered minus other genes (95% CI)")
save(style(fig, 1100, 380), "Figure5")

# ---------------------------------------------------------------- Supplementary figures: Tahoe and Orion cell numbers
tp = pd.read_csv(O + "tahoe_t1_prism.tsv.gz", sep="\t")
tp["vb"] = pd.cut(tp.v, [-0.01, 0.25, 0.5, 0.75, 0.95, 1.3], labels=["< 0.25", "0.25 to 0.5", "0.5 to 0.75", "0.75 to 0.95", "> 0.95"])
ts = pd.read_csv(O + "tahoe_t1_survival.tsv.gz", sep="\t")
w = ts.groupby(["sample", "plate", "drug", "conc"], as_index=False).n.sum()
dm = w[w.drug.str.contains("DMSO")].groupby("plate").n.mean().rename("dm"); w = w.join(dm, on="plate"); w["S"] = w.n / w.dm
wt = w[~w.drug.str.contains("DMSO")]
r1, p1 = spearmanr(tp.rel, tp.v); r2, p2 = spearmanr(tp.S, tp.v)
N["tahoe"] = dict(pairs=len(tp), drugs=int(tp.key.nunique()), rho_rel=round(r1, 3), p_rel=float(p1), rho_S=round(r2, 3), p_S=float(p2),
                  well_S_by_dose={str(c): round(wt.S[wt.conc == c].median(), 2) for c in sorted(wt.conc.dropna().unique())},
                  dmso_well_ratio_max=round(w[w.drug.str.contains("DMSO")].groupby("plate").n.agg(lambda v: v.max() / v.min()).max(), 2))
fig = make_subplots(rows=1, cols=2, subplot_titles=("a  Tahoe-100M, cell share against PRISM viability", "b  Tahoe-100M, cells per well against dose"), horizontal_spacing=0.12)
for j, (col, lab) in enumerate([("rel", "within-well share / DMSO share"), ("S", "cells / plate DMSO mean")]):
    fig.add_box(x=tp.vb, y=tp[col], name=lab, marker_color=C[j], boxpoints=False, row=1, col=1)
fig.add_box(x=wt.conc.map(lambda c: f"{c:g} µM"), y=wt.S, marker_color=C[2], boxpoints=False, showlegend=False, row=1, col=2)
fig.update_layout(boxmode="group"); fig.update_yaxes(range=[0, 2.5])
fig.update_xaxes(title_text="PRISM viability at the Tahoe dose", categoryorder="array", categoryarray=["> 0.95", "0.75 to 0.95", "0.5 to 0.75", "0.25 to 0.5", "< 0.25"], row=1, col=1); fig.update_xaxes(title_text="dose", categoryorder="array", categoryarray=[f"{c:g} µM" for c in sorted(wt.conc.dropna().unique())], row=1, col=2)
fig.update_yaxes(title_text="relative cells", row=1, col=1)
save(style(fig, 1100, 450), "FigureS1")

fig = go.Figure(); N["orion_cells"] = {}
g8 = pd.read_csv(os.path.expanduser("data/dosecurve_cache/gwps/groups.tsv"), sep="\t", header=None, names=["g", "o", "n"])
g8["chronos"] = g8.g.map(K562); g8 = g8.dropna(); g8["rel"] = g8.n / g8.n[g8.chronos > -0.3].median()
sets = {"K562 genome-wide, day 8": g8}
for line in ["HCT116", "HEK293T"]:
    sets[f"{line} genome-wide"] = pd.read_csv(O + f"orion_o1_{line}.tsv", sep="\t")
for i, (lab, d) in enumerate(sets.items()):
    b = pd.cut(d.chronos, BINS, labels=BL); m = d.groupby(b, observed=True).rel.median()
    N["orion_cells"][lab] = {k_: round(v, 2) for k_, v in m.items()}
    fig.add_scatter(x=list(m.index), y=m.values, mode="lines+markers", name=lab, marker=dict(symbol="diamond", size=9, color=C[[0, 4, 5][i]]), line=dict(color=C[[0, 4, 5][i]]))
fig.add_hline(y=1, line=dict(color="#7F7F7F", dash="dot", width=1))
fig.update_xaxes(title_text="Chronos gene effect of the target"); fig.update_yaxes(title_text="cells relative to non-essential targets", range=[0, 1.15])
save(style(fig, 900, 450), "FigureS2")

json.dump(N, open(os.path.join(R, "figures/numbers.json"), "w"), indent=1, default=float)
print(json.dumps({k_: v for k_, v in N.items() if k_ in ("k562_pairs", "k562_authors_knockdown", "tahoe", "orion_cells")}, indent=1, default=float)[:3000])
for lab, v in N["selection_share"].items():
    print(lab, v["n"], v["rho"], {b: v["by_bin"][b]["median"] for b in BL})
print(pd.DataFrame(N["dose_labels"]).to_string())
