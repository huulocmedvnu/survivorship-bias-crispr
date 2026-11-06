"""Do dose-response labels (Song 2025 buffered/sensitive, Mixscale, joint EM from dosecurve) track survivorship?
Reads data/dosecurve/regime_v4.tsv (K562 essential, 1,033 phenotype groups, labels recomputed in
the dosecurve project) read-only, joined with Phase P outputs."""
import os
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
O = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outputs/")
lab = pd.read_csv(os.path.expanduser("data/dosecurve/regime_v4.tsv"), sep="\t")
k = pd.read_csv(O + "phaseP_killtest.tsv", sep="\t")[["pert", "R", "rep", "m6", "chronos", "kd6", "kd8"]]
b = pd.read_csv(O + "phaseP_bounds.tsv", sep="\t")[["pert", "S", "m", "d_top", "sd_pert", "n_pert"]]
d = lab.merge(k, left_on="gene", right_on="pert").merge(b, on="pert")
d["share"] = (d.d_top / d.m).where((d.m > 2 * d.sd_pert / np.sqrt(d.n_pert)) & (d.m6 > 0)).clip(upper=1)
d["dkd"] = d.kd8 - d.kd6
print(f"joined: {len(d)} genes")
for meth in ["song", "mix", "em"]:
    print(f"\n== {meth}")
    g = d.groupby(d[meth].fillna("unlabelled"))
    print(g.agg(n=("gene", "size"), chronos=("chronos", "median"), S=("S", "median"), R=("R", "median"),
                share=("share", "median"), dkd=("dkd", "median")).round(3).to_string())
    x = d[d[meth].isin(["buffered", "sensitive"])].copy()
    if x[meth].nunique() == 2 and (x[meth] == "sensitive").sum() >= 10:
        x["buf"] = (x[meth] == "buffered").astype(int)
        for y in ["S", "R", "share"]:
            z = x.dropna(subset=[y, "chronos"])
            f = smf.ols(f"{y} ~ buf + chronos", data=z).fit()
            print(f"  {y} ~ buffered + chronos: buffered beta {f.params['buf']:.3f}, p {f.pvalues['buf']:.2g}, n {len(z)}")
d.to_csv(O + "phaseP_dose_labels.tsv", sep="\t", index=False)
