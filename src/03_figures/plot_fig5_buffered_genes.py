"""Figure 5: buffered minus other genes in the decrease of knockdown and in the response ratio, day 6 to day 8."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plotly.graph_objects as go  # noqa: E402
from plotly.subplots import make_subplots  # noqa: E402
from analysis import QL, BL, compute_all  # noqa: E402
from plot_style import C, BLACK, style, save  # noqa: E402

N, P = compute_all()
t5 = P["fig5"]
fig = make_subplots(rows=1, cols=2, subplot_titles=("a  Decrease in on-target knockdown, day 6 to day 8", "b  Downstream response ratio, day 8 / day 6"), horizontal_spacing=0.25)
for j, oc in enumerate(t5.outcome.unique()):
    s = t5[t5.outcome == oc]
    fig.add_scatter(x=s.beta, y=s.method, mode="markers", marker=dict(symbol="diamond", size=11, color=C[j]),
                    error_x=dict(type="data", symmetric=False, array=s.hi - s.beta, arrayminus=s.beta - s.lo, color=C[j], thickness=1.2), showlegend=False, row=1, col=j + 1)
    fig.add_vline(x=0, line=dict(color="#7F7F7F", dash="dot", width=1), row=1, col=j + 1)
fig.update_xaxes(title_text="buffered minus other genes (95% CI)")
save(style(fig, 1100, 380), "Figure5")
