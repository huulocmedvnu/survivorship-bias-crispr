"""Figure 2: cells, on-target knockdown and downstream response of K562 perturbations at day 6 and day 8."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plotly.graph_objects as go  # noqa: E402
from plotly.subplots import make_subplots  # noqa: E402
from analysis import QL, BL, compute_all  # noqa: E402
from plot_style import C, BLACK, style, save  # noqa: E402

N, P = compute_all()
p = P["fig2"]
fig = make_subplots(rows=1, cols=3, subplot_titles=("a  Cells per perturbation", "b  On-target knockdown (authors)", "c  Downstream response, day 8 / day 6"),
                    horizontal_spacing=0.09)
for j, (day, col) in enumerate([("n6", C[0]), ("n8", C[1])]):
    fig.add_bar(x=QL, y=p["cells"][day], name="day 6" if j == 0 else "day 8", marker_color=col, row=1, col=1)
for j, (day, col) in enumerate([("d6", C[0]), ("d8", C[1])]):
    fig.add_bar(x=QL, y=p["knockdown"][day], marker_color=col, showlegend=False, row=1, col=2)
for key, name, sym, size, col in [("ratio", "day 8 / day 6", "diamond", 11, C[1]), ("placebo", "day 6 half / half", "circle", 9, "#7F7F7F")]:
    med = p[key]
    fig.add_scatter(x=QL, y=[m[0] for m in med], mode="markers", marker=dict(symbol=sym, size=size, color=col), name=name,
                    error_y=dict(type="data", symmetric=False, array=[m[2] - m[0] for m in med], arrayminus=[m[0] - m[1] for m in med], color=col, thickness=1.2),
                    row=1, col=3)
fig.add_hline(y=1, line=dict(color="#7F7F7F", dash="dot", width=1), row=1, col=3)
fig.update_yaxes(title_text="median cells", row=1, col=1); fig.update_yaxes(title_text="median on-target knockdown (%)", range=[0, 100], row=1, col=2)
fig.update_yaxes(title_text="median ratio", range=[0, 1.3], row=1, col=3)
fig.update_xaxes(title_text="Chronos quartile in K562", tickangle=-30)
fig.update_layout(barmode="group")
save(style(fig, 1300, 460), "Figure2")
