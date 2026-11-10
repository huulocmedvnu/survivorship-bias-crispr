"""Figure 3: cells per guide and on-target knockdown over the Drepanos timecourse harvested on one day."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plotly.graph_objects as go  # noqa: E402
from plotly.subplots import make_subplots  # noqa: E402
from analysis import QL, BL, compute_all  # noqa: E402
from plot_style import C, BLACK, style, save  # noqa: E402

N, P = compute_all()
fig = make_subplots(rows=1, cols=3, subplot_titles=("a  Cells per guide, relative", "b  On-target knockdown, CRISPR interference", "c  On-target knockdown, CRISPR knockout"),
                    horizontal_spacing=0.09)
for g in P["fig3"]:
    col, dash = (C[1] if g["essential"] else C[0]), ("solid" if g["mod"] == "CRISPRi" else "dash")
    rel, sk = g["rel"], g["kd"]
    fig.add_scatter(x=rel.index, y=rel.values, mode="lines+markers", line=dict(color=col, dash=dash), name=g["label"], row=1, col=1)
    fig.add_scatter(x=sk.index, y=sk.kd.values, mode="lines+markers+text", text=[str(v) for v in sk.guides],
                    textposition="top center" if g["essential"] else "bottom center", textfont=dict(size=11, color=BLACK),
                    line=dict(color=col, dash=dash), showlegend=False, row=1, col=2 if g["mod"] == "CRISPRi" else 3)
fig.update_xaxes(title_text="days after transduction", tickvals=[4, 7, 10, 14], range=[3, 15])
fig.update_yaxes(title_text="cells relative to non-essential guides", range=[-0.05, 1.1], row=1, col=1)
fig.update_yaxes(title_text="target log2 fold change", range=[-2.4, -0.3], row=1, col=2)
fig.update_yaxes(range=[-1.3, -0.3], row=1, col=3)
fig.update_layout(legend=dict(orientation="h", y=-0.28, x=0))
save(style(fig, 1300, 500), "Figure3")
