"""Shared figure style: plotly_white, default Plotly colours, Arial, black text, no gridlines."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import FIGURES  # noqa: E402
import plotly.io as pio  # noqa: E402

pio.kaleido.scope.mathjax = None   # otherwise kaleido can print "Loading [MathJax]" into PDF exports

C = ["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3", "#FF6692"]
BLACK = "#111111"


def style(fig, w, h):
    fig.update_layout(template="plotly_white", width=w, height=h, font=dict(family="Arial", size=13, color=BLACK),
                      margin=dict(l=60, r=20, t=40, b=60), legend=dict(font=dict(color=BLACK)))
    for ax in (fig.update_xaxes, fig.update_yaxes):
        ax(showgrid=False, showline=True, linecolor=BLACK, ticks="outside", tickfont=dict(color=BLACK),
           title_font=dict(color=BLACK))
    for a in fig.layout.annotations:
        a.font = dict(family="Arial", size=14, color=BLACK)
    return fig


def save(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.write_image(FIGURES / f"{name}.png", scale=4)
    fig.write_image(FIGURES / f"{name}.pdf")
    print(f"wrote output/figures/{name}.png and .pdf")
