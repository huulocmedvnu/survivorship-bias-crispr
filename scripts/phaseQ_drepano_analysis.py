"""Phase Q (Drepanos 2026 K562 timecourse, days 4-14 processed together): per guide and day, relative cell number,
on-target knockdown and effect magnitude along the day-4 direction."""
import os
import numpy as np
import pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
PB = os.path.join(D, "ext/drepano/pb/")
OUT = os.path.join(D, "outputs/")
DEPMAP = os.path.expanduser("data/depmap/CRISPRGeneEffect.csv")
DAYS, EPS, MIN_CTRL, TOP = [4, 7, 10, 14], 0.1, 0.5, 200

ge = pd.read_csv(DEPMAP, index_col=0)
ge.columns = [c.split(" (")[0] for c in ge.columns]
chronos = ge.loc["ACH-000551"]

pb = {d: {p: pd.read_csv(PB + f"day{d}_{p}.csv.gz", index_col=0) for p in ["all", "h0", "h1"]} for d in DAYS}
ncell = {d: pd.read_csv(PB + f"day{d}_n.csv").set_index("group").n for d in DAYS}
genes = sorted(set.intersection(*[set(pb[d]["all"].index) for d in DAYS]))


def ctrl_col(d, mod):
    cols = [c for c in pb[d]["all"].columns if c.startswith("CONTROL_")]
    return cols[0] if len(cols) == 1 else [c for c in cols if c.endswith("_" + mod)][0]


# genes expressed in every day's controls
expr = np.ones(len(genes), bool)
for d in DAYS:
    for mod in ["CRISPRi", "CRISPRko"]:
        expr &= pb[d]["all"].loc[genes, ctrl_col(d, mod)].values >= MIN_CTRL
genes = np.array(genes)[expr]
print(f"genes expressed in all controls: {len(genes)}")
lfc = lambda x, c: np.log2((x + EPS) / (c + EPS))

guides = [c for c in pb[4]["all"].columns if not c.startswith("CONTROL_")]
rows = []
for gcol in guides:
    mod, rest = gcol.split("|", 1)
    gene = rest.split("_")[0]
    c4 = pb[4]["all"].loc[genes, ctrl_col(4, mod)].values
    if gcol not in pb[4]["h0"] or gcol not in pb[4]["h1"]:
        continue
    l0 = lfc(pb[4]["h0"].loc[genes, gcol].values, c4)
    l1 = lfc(pb[4]["h1"].loc[genes, gcol].values, c4)
    mask = genes != gene
    dirs = []
    for ls, lm in ((l0, l1), (l1, l0)):
        idx = np.where(mask)[0][np.argsort(-np.abs(ls[mask]))[:TOP]]
        dirs.append((idx, ls[idx] / np.linalg.norm(ls[idx]), lm[idx] @ (ls[idx] / np.linalg.norm(ls[idx]))))
    m4 = np.mean([x[2] for x in dirs])
    rep = np.mean([np.corrcoef(l0[x[0]], l1[x[0]])[0, 1] for x in dirs])
    for d in DAYS:
        if gcol not in pb[d]["all"]:
            rows.append(dict(guide=gcol, gene=gene, mod=mod, day=d, n=0))
            continue
        c = pb[d]["all"].loc[genes, ctrl_col(d, mod)].values
        l = lfc(pb[d]["all"].loc[genes, gcol].values, c)
        md = np.mean([l[x[0]] @ x[1] for x in dirs])
        kd = lfc(pb[d]["all"].loc[[gene], gcol].values, pb[d]["all"].loc[[gene], ctrl_col(d, mod)].values)[0] \
            if gene in pb[d]["all"].index else np.nan
        rows.append(dict(guide=gcol, gene=gene, mod=mod, day=d, n=int(ncell[d].get(gcol, 0)), kd=kd,
                         m=md, m4=m4, rep4=rep))
t = pd.DataFrame(rows)
t["chronos"] = t.gene.map(chronos)
t["essential"] = t.gene.isin(["RPL9", "DCAF13"])
ref = t[~t.essential].groupby(["day", "mod"]).n.median().rename("n_ref")
t = t.join(ref, on=["day", "mod"])
t["rel_n"] = t.n / t.n_ref
t["ratio"] = t.m / t.m4
t.to_csv(OUT + "phaseQ_drepano.tsv", sep="\t", index=False)

print("Chronos of library genes (K562):", t.drop_duplicates("gene").set_index("gene").chronos.round(2).sort_values().head(6).to_dict())
for mod in ["CRISPRi", "CRISPRko"]:
    s = t[t["mod"] == mod]
    print(f"\n== {mod}: median per day, essential (RPL9, DCAF13) vs non-essential guides")
    agg = s.groupby(["essential", "day"]).agg(guides=("guide", "nunique"), rel_n=("rel_n", "median"),
                                             kd=("kd", "median"), effect_ratio=("ratio", "median")).round(2)
    print(agg.to_string())
    ne = s[(~s.essential) & (s.rep4 > 0.3)]
    q = ne.groupby("day").ratio.quantile([0.05, 0.95]).unstack().round(2)
    print("non-essential guides with reproducible day-4 signature, 5th-95th percentile of effect ratio:")
    print(q.to_string())
    print("essential guides, each:")
    print(s[s.essential].pivot_table(index="guide", columns="day", values=["rel_n", "ratio"]).round(2).to_string())
