"""Phase O1: X-Atlas/Orion cells per target gene from metadata columns only (HTTP range reads, IPv4),
checkpointed per file in tahoe/orion_meta/, then compared with DepMap Chronos."""
import os, sys, json, time
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np, pandas as pd, pyarrow.parquet as pq, requests
from scipy.stats import spearmanr
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from httpfile import HTTPFile, S
REPO = "https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion"
CK = os.path.join(H, "orion_meta"); os.makedirs(CK, exist_ok=True)
COLS = ["gene_target", "guide_target", "total_counts", "pass_guide_filter"]

files = [x["path"] for x in S.get(REPO.replace("/datasets/", "/api/datasets/") + "/tree/main/data", timeout=60).json() if x["type"] == "file"]
todo = [f for f in files if not os.path.exists(os.path.join(CK, os.path.basename(f) + ".tsv"))]
print(f"{len(files)} files, {len(todo)} to read")

def one(path):
    f = HTTPFile(f"{REPO}/resolve/main/{path}")
    t = pq.ParquetFile(f).read(columns=COLS).to_pandas()
    t = t[t.pass_guide_filter == 1]
    g = t.groupby(["gene_target", "guide_target"]).agg(n=("total_counts", "size"), umi=("total_counts", "median")).reset_index()
    out = os.path.join(CK, os.path.basename(path) + ".tsv")
    g.to_csv(out + ".tmp", sep="\t", index=False); os.replace(out + ".tmp", out)        # atomic: a killed run leaves no half file
    return path, f.bytes_read

t0 = time.time(); got = 0
with ThreadPoolExecutor(8) as ex:
    for k, fut in enumerate(as_completed([ex.submit(one, p) for p in todo]), 1):
        try:
            p, b = fut.result(); got += b
        except Exception as e:
            print("failed (will retry next run):", e)
        if k % 25 == 0:
            print(f"  {k}/{len(todo)} files, {got/1e6:.0f} MB, {time.time()-t0:.0f}s", flush=True)

parts = []
for fn in sorted(os.listdir(CK)):
    if fn.endswith(".tsv"):
        x = pd.read_csv(os.path.join(CK, fn), sep="\t"); x["line"] = fn.split("_")[0]; parts.append(x)
a = pd.concat(parts)
print(f"files read: {len(parts)} of {len(files)}")
ge = pd.read_csv("data/depmap/CRISPRGeneEffect.csv", index_col=0)
ge.columns = [c.split(" (")[0] for c in ge.columns]
IDS = {"HCT116": "ACH-000971", "HEK293T": sys.argv[1] if len(sys.argv) > 1 else None}
rows = []
for line, d in a.groupby("line"):
    g = d.groupby("gene_target").n.sum()
    ntc = g.get("Non-Targeting", np.nan); g = g.drop("Non-Targeting", errors="ignore")
    if IDS.get(line) in ge.index:
        ch = ge.loc[IDS[line]]
    else:                                                     # HEK293T has no DepMap CRISPR profile
        ch = ge.median(); print(f"{line}: no own DepMap profile, median Chronos over {len(ge)} lines used")
    x = pd.DataFrame({"n": g, "chronos": g.index.map(ch)}).dropna()
    ref = x.n[x.chronos > -0.3].median()
    x["rel"] = x.n / ref
    r, p = spearmanr(np.log(x.n), x.chronos)
    q = pd.qcut(x.chronos, 4, labels=["Q1 most ess", "Q2", "Q3", "Q4 least"])
    print(f"\n{line}: {len(x)} genes with Chronos, cells {int(x.n.sum()):,}, Non-Targeting {ntc:,.0f}, reference {ref:.0f} cells")
    print(f"  Spearman log cells vs Chronos: rho {r:.3f}, p {p:.2g}")
    print(x.groupby(q, observed=True).agg(genes=("n", "size"), chronos=("chronos", "median"), cells=("n", "median"), rel=("rel", "median")).round(2).to_string())
    x.assign(line=line).to_csv(os.path.join(H, "..", "outputs", f"orion_o1_{line}.tsv"), sep="\t")
