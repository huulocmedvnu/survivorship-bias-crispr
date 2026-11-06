"""Reference check of tools/survb against plain numpy on a few perturbations (uses the local gwps cache)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from selection_bound import d_max_from_projections
C = os.path.expanduser("data/dosecurve_cache/gwps/")
O = os.path.join(os.path.dirname(__file__), "..", "outputs/")
ip = np.memmap(C + "indptr.bin", dtype=np.uint64, mode="r")
idx = np.memmap(C + "idx.bin", dtype=np.uint16, mode="r")
cnt = np.memmap(C + "cnt.bin", dtype=np.uint16, mode="r")
lib = np.fromfile(C + "lib.bin", dtype=np.float32)
mem = np.fromfile(C + "members.bin", dtype=np.uint32)
gr = pd.read_csv(C + "groups.tsv", sep="\t", header=None, names=["g", "off", "n"]).set_index("g")
dirs = pd.read_csv(O + "phaseP_dirs.tsv", sep="\t")
pert = pd.read_csv(O + "phaseP_pert.tsv", sep="\t").set_index("pert")
rust = pd.read_csv(O + "phaseP_bounds.tsv", sep="\t").set_index("pert")
ng = sum(1 for _ in open(C + "genes.txt"))
def scores(cells, w):
    out = np.empty(len(cells))
    for k, c in enumerate(cells):
        a, b = int(ip[c]), int(ip[c + 1])
        x = np.zeros(ng); x[idx[a:b]] = np.log1p(cnt[a:b] * (1e4 / lib[c]))
        out[k] = x @ w
    return out
o, n = gr.loc["control"]
ctrl = mem[o:o + n]
for p in sys.argv[1:]:
    w = np.zeros(ng); d = dirs[dirs.pert == p]; w[d.gene_index.values] = d.weight.values
    yc = scores(ctrl, w)
    o2, n2 = gr.loc[p]; yp = scores(mem[o2:o2 + n2], w)
    S = pert.S[p]
    ref = dict(m=yp.mean() - yc.mean(), d_top=d_max_from_projections(yc, S), d_bot=-d_max_from_projections(-yc, S))
    print(p, {k: (round(v, 5), round(rust.loc[p, k], 5)) for k, v in ref.items()})
