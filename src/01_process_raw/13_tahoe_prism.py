"""Phase T1: Tahoe-100M survival per (well, cell line) from cell counts, checked against PRISM viability."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import DATA, SOURCE, SCORER  # noqa: E402
import os, re
import numpy as np, pandas as pd, pyarrow.parquet as pq
from scipy.stats import spearmanr
H = str(DATA / "tahoe")
DEP = f"{DATA}/depmap/"
OUT = f"{SOURCE}/"
norm = lambda s: re.sub(r"[^a-z0-9]", "", re.sub(r"\(.*?\)", "", str(s).lower()))

pf = pq.ParquetFile(os.path.join(H, "obs_metadata.parquet"))
print("obs columns:", pf.schema_arrow.names, "rows", pf.metadata.num_rows)
import pyarrow as pa, pyarrow.compute as pc, time
t0 = time.time()
parts = []
for rg in range(pf.num_row_groups):                                   # Arrow group_by in C++, one row group at a time
    tb = pf.read_row_group(rg, columns=["sample", "plate", "drug", "cell_name", "pass_filter"])
    tb = tb.filter(pc.equal(pc.cast(tb["pass_filter"], pa.string()), "full"))
    tb = tb.set_column(3, "cell_name", pc.cast(tb["cell_name"], pa.string())).drop(["pass_filter"])
    parts.append(tb.group_by(["sample", "plate", "drug", "cell_name"]).aggregate([("sample", "count")]))
counts = pa.concat_tables(parts).group_by(["sample", "plate", "drug", "cell_name"]).aggregate([("sample_count", "sum")]).to_pandas()
counts = counts.rename(columns={"sample_count_sum": "n"})
print(f"counted {counts.n.sum():,} passing cells in {time.time() - t0:.1f}s")
cols = ["sample", "plate", "drug", "cell_name"]
line_col = next(c for c in ["cell_line_id", "cell_line", "cell_name"] if c in cols)
sm = pd.read_parquet(os.path.join(H, "sample_metadata.parquet"))
sm["conc"] = sm.drugname_drugconc.str.extract(r",\s*([0-9.eE+-]+)\s*,").astype(float)
c = counts.merge(sm[["sample", "conc"]], on="sample", how="left")
if "drug" not in c: c = c.merge(sm[["sample", "drug"]], on="sample")
if "plate" not in c: c = c.merge(sm[["sample", "plate"]], on="sample")
c["well_n"] = c.groupby("sample").n.transform("sum")
c["share"] = c.n / c.well_n
dm = c[c.drug.str.contains("DMSO", na=False)]
ref = dm.groupby(["plate", line_col]).agg(dmso_share=("share", "mean"), dmso_n=("n", "mean")).reset_index()
c = c.merge(ref, on=["plate", line_col], how="left")
c["rel"] = c.share / c.dmso_share
c["S"] = c.n / c.dmso_n
cl = pd.read_parquet(os.path.join(H, "cell_line_metadata.parquet")).drop_duplicates("cell_name")
idmap = dict(zip(cl.Cell_ID_Cellosaur, cl.Cell_ID_DepMap)) if line_col == "cell_line_id" else dict(zip(cl.cell_name, cl.Cell_ID_DepMap))
c["depmap_id"] = c[line_col].map(idmap)
c["key"] = c.drug.map(norm)
t = c[~c.drug.str.contains("DMSO", na=False)].copy()
print(f"treated wells x lines: {len(t)}, median cells {t.n.median():.0f}, DMSO wells x lines {len(dm)}")
p = pd.read_csv(DEP + "secondary-screen-dose-response-curve-parameters.csv", low_memory=False)
p = p[(p.r2 >= 0.3) & np.isfinite(p.ec50) & (p.ec50 > 0)].copy()
p["key"] = p.name.map(norm)
m = t.merge(p[["key", "depmap_id", "upper_limit", "lower_limit", "slope", "ec50"]], on=["key", "depmap_id"])
m["v"] = (m.lower_limit + (m.upper_limit - m.lower_limit) / (1 + (m.conc / m.ec50) ** m.slope)).clip(0, 1.2)
m = m.groupby(["sample", line_col, "key", "conc", "n", "rel", "S"], as_index=False).v.mean()
m = m[np.isfinite(m.rel) & np.isfinite(m.S)]
print(f"matched (well x line) with PRISM: {len(m)}, drugs {m.key.nunique()}, lines {m[line_col].nunique()}")
for k in ["rel", "S"]:
    r, pv = spearmanr(m[k], m.v)
    w = m.groupby("key").apply(lambda g: spearmanr(g[k], g.v)[0] if g[line_col].nunique() >= 8 and g.v.std() > 0.05 else np.nan).dropna()
    print(f"{k}: overall rho {r:.3f} (p {pv:.2g}) | within-drug median rho {w.median():.3f} over {len(w)} drugs, positive in {(w > 0).mean():.0%}")
m["vq"] = pd.cut(m.v, [-0.01, 0.25, 0.5, 0.75, 0.95, 1.3])
print(m.groupby("vq", observed=True).agg(n=("rel", "size"), rel=("rel", "median"), S=("S", "median")).round(2).to_string())
c.to_csv(OUT + "tahoe_t1_survival.tsv.gz", sep="\t", index=False)
m.to_csv(OUT + "tahoe_t1_prism.tsv.gz", sep="\t", index=False)
