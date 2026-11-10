"""Figs. S1 and S2: cell numbers in Tahoe-100M against PRISM viability, and cells per perturbation by essentiality."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plotly.graph_objects as go  # noqa: E402
from plotly.subplots import make_subplots  # noqa: E402
from analysis import QL, BL, compute_all  # noqa: E402
from plot_style import C, BLACK, style, save  # noqa: E402

N, P = compute_all()
tp, wt = P["figS1"]["prism"], P["figS1"]["wells"]
fig = make_subplots(rows=1, cols=2, subplot_titles=("a  Tahoe-100M, cell share against PRISM viability", "b  Tahoe-100M, cells per well against dose"), horizontal_spacing=0.12)
for j, (col, lab) in enumerate([("rel", "within-well share / DMSO share"), ("S", "cells / plate DMSO mean")]):
    fig.add_box(x=tp.vb, y=tp[col], name=lab, marker_color=C[j], boxpoints=False, row=1, col=1)
fig.add_box(x=wt.conc.map(lambda c: f"{c:g} µM"), y=wt.S, marker_color=C[2], boxpoints=False, showlegend=False, row=1, col=2)
fig.update_layout(boxmode="group"); fig.update_yaxes(range=[0, 2.5])
fig.update_xaxes(title_text="PRISM viability at the Tahoe dose", categoryorder="array", categoryarray=["> 0.95", "0.75 to 0.95", "0.5 to 0.75", "0.25 to 0.5", "< 0.25"], row=1, col=1)
fig.update_xaxes(title_text="dose", categoryorder="array", categoryarray=[f"{c:g} µM" for c in sorted(wt.conc.dropna().unique())], row=1, col=2)
fig.update_yaxes(title_text="relative cells", row=1, col=1)
save(style(fig, 1100, 450), "FigureS1")

fig = go.Figure()
for i, (lab, m) in enumerate(P["figS2"].items()):
    fig.add_scatter(x=list(m.index), y=m.values, mode="lines+markers", name=lab, marker=dict(symbol="diamond", size=9, color=C[[0, 4, 5][i]]), line=dict(color=C[[0, 4, 5][i]]))
fig.add_hline(y=1, line=dict(color="#7F7F7F", dash="dot", width=1))
fig.update_xaxes(title_text="Chronos gene effect of the target"); fig.update_yaxes(title_text="cells relative to non-essential targets", range=[0, 1.15])
save(style(fig, 900, 450), "FigureS2")
