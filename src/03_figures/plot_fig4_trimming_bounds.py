"""Figure 4: median fraction of the transcriptional shift explainable by selective survival, by Chronos bin, in seven screens."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plotly.graph_objects as go  # noqa: E402
from plotly.subplots import make_subplots  # noqa: E402
from analysis import QL, BL, compute_all  # noqa: E402
from plot_style import C, BLACK, style, save  # noqa: E402

N, P = compute_all()
fig = go.Figure()
for i, (lab, st) in enumerate(P["fig4"].items()):
    ok = [j for j, v in enumerate(st) if v[3] >= 15]
    fig.add_scatter(x=[BL[j] for j in ok], y=[st[j][0] for j in ok], mode="lines+markers", name=lab, line=dict(color=C[i % 7], width=1.5),
                    marker=dict(symbol="diamond", size=9),
                    error_y=dict(type="data", symmetric=False, array=[st[j][2] - st[j][0] for j in ok], arrayminus=[st[j][0] - st[j][1] for j in ok], thickness=1, width=3))
fig.update_xaxes(title_text="Chronos gene effect of the target in the screened line (bins with at least 15 targets)")
fig.update_yaxes(title_text="median fraction explainable by selective survival", range=[-0.05, 1.05])
save(style(fig, 1000, 520), "Figure4")
