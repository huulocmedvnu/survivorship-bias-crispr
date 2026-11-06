"""Phase Q, Dixit 2016 (GSE90063): K562 CRISPRko TF library, same library at day 7 and day 13.
Per target gene: cells (relative to the median of non-essential targets), on-target mRNA, and day-13 / day-7
effect ratio along the day-7 direction (genes chosen in one day-7 half, measured in the other half and at day 13).
Cells with exactly one guide only. Pre-registered (Phase Q) before results."""
import os, re
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
X = os.path.join(D, "ext/dixit/")
DEPMAP = os.path.expanduser("data/depmap/CRISPRGeneEffect.csv")
EPS, MIN_CTRL, TOP = 0.1, 0.5, 200


def load(tag):
    base = X + tag
    genes = pd.read_csv(base + "_genenames.csv.gz", header=None, skiprows=1).iloc[:, 1].str.split("_", n=1).str[1].values
    cells = pd.read_csv(base + "_cellnames.csv.gz", header=None, skiprows=1).iloc[:, 1].values
    tri = pd.read_csv(base + ".mtx.txt.gz", sep=" ", comment="%", header=None, skiprows=3, dtype=np.float32, engine="c")
    nr, nc = len(genes), len(cells)
    m = sp.csc_matrix((tri[2].values, (tri[0].values.astype(np.int64) - 1, tri[1].values.astype(np.int64) - 1)), shape=(nr, nc))
    lib = np.asarray(m.sum(0)).ravel()
    m = m @ sp.diags(1e4 / np.maximum(lib, 1))                    # CP10k
    guides = {}
    for line in pd.read_csv(base + "_cbc_gbc_dict.csv.gz", header=None).itertuples(index=False):
        for c in str(line[1]).split(","):
            guides.setdefault(c.strip(), []).append(line[0])
    gpos = {c: i for i, c in enumerate(cells)}
    target = {}
    for c, gs in guides.items():
        if len(gs) == 1 and c in gpos:
            g = gs[0]
            target[gpos[c]] = "CONTROL" if re.search(r"INTERGENIC|NO_SITE|non.?target", g, re.I) else g.split("_")[1].replace("sg", "")
    print(f"{tag}: {nr} genes x {nc} cells, {len(target)} single-guide cells")
    return genes, m.tocsc(), pd.Series(target)


def pseudobulk(m, cols):
    cols = np.asarray(cols)
    return np.asarray(m[:, cols].mean(1)).ravel()


g7, m7, t7 = load("GSM2396858_k562_tfs_7")
g13, m13, t13 = load("GSM2396859_k562_tfs_13")
print("targets day 7:", t7.value_counts().to_dict())
common = sorted(set(g7) & set(g13))
i7 = pd.Series(range(len(g7)), index=g7).groupby(level=0).first()[common].values
i13 = pd.Series(range(len(g13)), index=g13).groupby(level=0).first()[common].values
c7 = pseudobulk(m7, t7.index[t7 == "CONTROL"])[i7]
c13 = pseudobulk(m13, t13.index[t13 == "CONTROL"])[i13]
keep = (c7 >= MIN_CTRL) & (c13 >= MIN_CTRL)
genes = np.array(common)[keep]
c7, c13, i7, i13 = c7[keep], c13[keep], i7[keep], i13[keep]
lfc = lambda x, c: np.log2((x + EPS) / (c + EPS))

ge = pd.read_csv(DEPMAP, index_col=0)
ge.columns = [c.split(" (")[0] for c in ge.columns]
chronos = ge.loc["ACH-000551"]
rows = []
for tg in sorted(set(t7) - {"CONTROL"}):
    cells7 = t7.index[t7 == tg].values
    cells13 = t13.index[t13 == tg].values
    if len(cells7) < 20:
        continue
    h0, h1 = cells7[0::2], cells7[1::2]
    l0 = lfc(pseudobulk(m7, h0)[i7], c7)
    l1 = lfc(pseudobulk(m7, h1)[i7], c7)
    l13 = lfc(pseudobulk(m13, cells13)[i13], c13) if len(cells13) else np.full(len(genes), np.nan)
    mask = genes != tg
    m_7, m_13, reps = [], [], []
    for ls, lm in ((l0, l1), (l1, l0)):
        idx = np.where(mask)[0][np.argsort(-np.abs(ls[mask]))[:TOP]]
        u = ls[idx] / np.linalg.norm(ls[idx])
        m_7.append(lm[idx] @ u); m_13.append(l13[idx] @ u); reps.append(np.corrcoef(ls[idx], lm[idx])[0, 1])
    kd7 = lfc(pseudobulk(m7, cells7)[i7], c7)[genes == tg]
    kd13 = lfc(pseudobulk(m13, cells13)[i13], c13)[genes == tg] if len(cells13) else []
    rows.append(dict(target=tg, n7=len(cells7), n13=len(cells13), m7=np.mean(m_7), m13=np.mean(m_13), rep7=np.mean(reps),
                     kd7=kd7[0] if len(kd7) else np.nan, kd13=kd13[0] if len(kd13) else np.nan, chronos=chronos.get(tg, np.nan)))
d = pd.DataFrame(rows)
nc7 = (t7 == "CONTROL").sum(); nc13 = (t13 == "CONTROL").sum()
d["rel_n"] = (d.n13 / nc13) / (d.n7 / nc7)          # change in share of cells, relative to control cells
d["R"] = d.m13 / d.m7
d.to_csv(os.path.join(D, "outputs/phaseQ_dixit.tsv"), sep="\t", index=False)
print(d.sort_values("chronos").round(2).to_string(index=False))
s = d[(d.rep7 > 0.3) & (d.m7 > 0)]
for v in ["R", "rel_n"]:
    z = s.dropna(subset=[v, "chronos"]) if v == "R" else d.dropna(subset=[v, "chronos"])
    r, p = spearmanr(z[v], z.chronos)
    print(f"Spearman {v} vs Chronos: rho {r:.2f}, p {p:.2g}, n {len(z)}")
