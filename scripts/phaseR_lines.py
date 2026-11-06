"""Phase R: trimming bounds in four essential-gene CRISPRi screens with one timepoint each
(K562 day 6, RPE1 day 7 from Replogle 2022, HepG2 and Jurkat from Nadig 2025).

Per line and perturbation: direction u from the even-position half (pseudobulk part h0, top 200 genes by |lfc|,
target excluded, genes with control mean >= 0.5 CP10k), observed shift m and trimming bound d_top measured with
the odd-position half of the perturbed cells (survb half 1) against all control cells.
S = cells / median cells of perturbations with Chronos > -0.3 in that line, capped at 1.
Chronos: K562 ACH-000551, HepG2 ACH-000739, Jurkat ACH-000995, RPE1 = mean of the RPE1SS lines if screened,
else the median over all lines (reported).
"""
import os, subprocess
import numpy as np, pandas as pd
from scipy.stats import spearmanr

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
PB = os.path.expanduser("data/pseudobulk/")
CACHE = os.path.expanduser("data/dosecurve_cache/")
DEP = os.path.expanduser("data/depmap/")
SURVB = os.path.join(D, "tools/survb/target/release/survb")
OUT = os.path.join(D, "outputs/")
EPS, MIN_CTRL, TOP = 0.1, 0.5, 200

ge = pd.read_csv(DEP + "CRISPRGeneEffect.csv", index_col=0)
ge.columns = [c.split(" (")[0] for c in ge.columns]
model = pd.read_csv(DEP + "Model.csv")
rpe = [i for i in model[model.StrippedCellLineName.str.startswith("RPE1", na=False)].ModelID if i in ge.index]
CHRONOS = {"k562ess": ge.loc["ACH-000551"], "hepg2": ge.loc["ACH-000739"], "jurkat": ge.loc["ACH-000995"],
           "rpe1": ge.loc[rpe].mean() if rpe else ge.median()}
print("RPE1 Chronos source:", f"mean of {len(rpe)} RPE1SS lines" if rpe else "median over all lines")

summary = []
for line in ["k562ess", "rpe1", "hepg2", "jurkat"]:
    genes = open(PB + f"{line}.genes.txt").read().split()
    gr = pd.read_csv(PB + f"{line}.groups.tsv", sep="\t", header=None, names=["g", "n"])
    a = np.memmap(PB + f"{line}.f32", dtype=np.float32, mode="r").reshape(3, len(gr), len(genes))
    assert gr.g[0] == "control"
    ctrl = np.asarray(a[0, 0], dtype=np.float64)
    ok = ctrl >= MIN_CTRL
    gidx = np.where(ok)[0]
    garr = np.array(genes)
    ch = CHRONOS[line]
    gr["chronos"] = gr.g.map(ch)
    ne = gr[(gr.g != "control") & (gr.chronos > -0.3)]
    ref = ne.n.median()
    dirs, perts = [], []
    for r in range(1, len(gr)):
        p = gr.g[r]
        if gr.n[r] < 20:
            continue
        l = np.log2((np.asarray(a[1, r, gidx], dtype=np.float64) + EPS) / (ctrl[gidx] + EPS))
        m = garr[gidx] != p
        idx = np.where(m)[0][np.argsort(-np.abs(l[m]))[:TOP]]
        u = l[idx] / np.linalg.norm(l[idx])
        dirs += [(p, int(gidx[j]), float(w)) for j, w in zip(idx, u)]
        perts.append((p, int(gr.n[r]), min(1.0, gr.n[r] / ref), ch.get(p, np.nan)))
    fd, fp, fb = (OUT + f"phaseR_{line}_{k}.tsv" for k in ("dirs", "pert", "bounds"))
    pd.DataFrame(dirs, columns=["pert", "gene_index", "weight"]).to_csv(fd, sep="\t", index=False)
    pd.DataFrame(perts, columns=["pert", "n", "S", "chronos"]).to_csv(fp, sep="\t", index=False)
    subprocess.run([SURVB, CACHE + line, fd, fp, fb, "1"], check=True, stderr=subprocess.DEVNULL)
    b = pd.read_csv(fb, sep="\t").merge(pd.read_csv(fp, sep="\t")[["pert", "chronos"]], on="pert")
    b["se"] = b.sd_pert / np.sqrt(b.n_pert)
    s = b[(b.m > 2 * b.se) & b.chronos.notna()].copy()
    s["share"] = (s.d_top / s.m).clip(upper=1)
    q = pd.qcut(s.chronos, 4, labels=False)
    r, pv = spearmanr(s.share, s.chronos)
    summary.append(dict(line=line, ref_cells=ref, n_ref=len(ne), perts=len(s), S_lt1=(s.S < 1).mean(),
                        share_median=s.share.median(), half_or_more=(s.share >= 0.5).mean(), full=(s.d_top >= s.m).mean(),
                        share_Q1_most_ess=s.share[q == 0].median(), share_Q4_least=s.share[q == 3].median(), rho=r, p=pv))
    os.remove(fd)                                   # large and rebuildable
t = pd.DataFrame(summary)
t.to_csv(OUT + "phaseR_summary.tsv", sep="\t", index=False)
print(t.round(3).to_string(index=False))
