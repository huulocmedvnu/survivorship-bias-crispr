"""Figure 1: study design as an SVG schematic, exported to PNG (600 dpi) and PDF.

Each row: data (left), a pictogram of the comparison (middle) and a small plot of the readout (right).
Marks in the readout plots use the estimates reported in the Results.
Canvas units: 10 per mm, 180 mm wide, so 8 pt is about 28 units.
Run: python3 figures/make_flowchart.py
"""
import json, math, os, random, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
F = os.path.join(HERE, "figures")
REL = json.load(open(os.path.join(HERE, "numbers.json")))["drepano"]["CRISPRi, RPL9 and DCAF13"]["rel_n"]
W = 1800
INK, GREY, LINE = "#111111", "#8C96A8", "#DCE1EA"
BLUE, RED, GREEN, PURPLE, ORANGE = "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A"
FONT = "Arial, Helvetica, sans-serif"
svg = []
add = svg.append


def text(x, y, s, size=27, weight="normal", anchor="start", fill=INK):
    add(f"<text x='{x:.1f}' y='{y:.1f}' font-family='{FONT}' font-size='{size}' font-weight='{weight}' "
        f"text-anchor='{anchor}' fill='{fill}'>{s}</text>")


def line(x1, y1, x2, y2, col=INK, w=3, dash=None, head=False):
    d = f" stroke-dasharray='{dash}'" if dash else ""
    h = " marker-end='url(#ah)'" if head else ""
    add(f"<line x1='{x1:.1f}' y1='{y1:.1f}' x2='{x2:.1f}' y2='{y2:.1f}' stroke='{col}' stroke-width='{w}'{d}{h}/>")


def cell(x, y, col, r=13, dead=False):
    if dead:
        add(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{r}' fill='none' stroke='{col}' stroke-width='2.5' stroke-dasharray='4 4'/>")
    else:
        add(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{r}' fill='{col}' fill-opacity='0.25' stroke='{col}' stroke-width='3'/>"
            f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{r * 0.35:.1f}' fill='{col}'/>")


def group(cx, cy, col, n=6, dead=0, cols=3, step=32, r=13):
    rows = math.ceil(n / cols)
    for k in range(n):
        i, j = k % cols, k // cols
        cell(cx + (i - (cols - 1) / 2) * step, cy + (j - (rows - 1) / 2) * step, col, r, dead=k >= n - dead)


def axes(x0, y0, w, h, xlab=None, ylab=None):
    """Left and bottom axis lines of a mini plot whose bottom-left corner is (x0, y0)."""
    line(x0, y0, x0 + w, y0, INK, 2.5); line(x0, y0, x0, y0 - h, INK, 2.5)
    if xlab:
        text(x0 + w / 2, y0 + 30, xlab, 24, anchor="middle")
    if ylab:
        add(f"<text transform='translate({x0 - 14},{y0 - h / 2}) rotate(-90)' font-family='{FONT}' font-size='24' "
            f"text-anchor='middle' fill='{INK}'>{ylab}</text>")


add("<svg xmlns='http://www.w3.org/2000/svg' width='180mm' height='@HMM@mm' viewBox='0 0 1800 @H@'>")
add(f"<defs><marker id='ah' viewBox='0 0 10 10' refX='8' refY='5' markerWidth='4' markerHeight='4' orient='auto'>"
    f"<path d='M0,0 L10,5 L0,10 Z' fill='{GREY}'/></marker></defs>")
add("<rect width='1800' height='@H@' fill='white'/>")

C1, C2, C3 = 40, 560, 1250
text(C1, 58, "Data", 34, "bold"); text(C2, 58, "Comparison", 34, "bold"); text(C3, 58, "Readout", 34, "bold")
RH, Y0 = 250, 90
FIGX = 1760


def row_frame(k, col, title, lines_, fig):
    y = Y0 + k * RH
    add(f"<rect x='{C1 - 20}' y='{y + 12}' width='8' height='{RH - 34}' fill='{col}'/>")
    text(C1, y + 44, title, 29, "bold", fill=col)
    for i, s in enumerate(lines_):
        text(C1, y + 84 + i * 34, s, 26)
    text(FIGX, y + 44, fig, 27, "bold", "end", col)
    if k:
        line(C1, y, W - 40, y, LINE, 2)
    return y


# ---------------------------------------------------------------- row 1: two timepoints, same guides
y = row_frame(0, BLUE, "K562, two screens", ["Replogle, essential genes, day 6", "Replogle, genome-wide, day 8", "Same guides for 2,062 genes"], "Fig. 2")
cy = y + 110
text(C2 + 105, y + 44, "day 6", 26, "bold", "middle")
group(C2 + 40, cy, BLUE, 4, cols=2); group(C2 + 170, cy, BLUE, 4, cols=2)
line(C2 + 105, cy - 45, C2 + 105, cy + 45, GREY, 2.5, "6 5")
text(C2 + 40, cy + 74, "half A", 24, anchor="middle"); text(C2 + 170, cy + 74, "half B", 24, anchor="middle")
text(C2 + 40, cy + 102, "200 genes", 24, anchor="middle"); text(C2 + 170, cy + 102, "magnitude", 24, anchor="middle")
line(C2 + 245, cy, C2 + 320, cy, GREY, 4, head=True)
text(C2 + 410, y + 44, "day 8", 26, "bold", "middle")
group(C2 + 410, cy, BLUE, 6, dead=2)
text(C2 + 410, cy + 74, "same 200 genes", 24, anchor="middle")
text(C2 + 410, cy + 102, "in survivors", 24, anchor="middle")
# readout: response ratio and knockdown in the most essential quartile
bx, by, bh = C3 + 30, y + 210, 120
axes(bx, by, 330, bh)
for i, (v, op) in enumerate([(1.0, 0.35), (0.68, 1.0)]):
    add(f"<rect x='{bx + 20 + i * 70}' y='{by - v * bh:.1f}' width='54' height='{v * bh:.1f}' fill='{BLUE}' fill-opacity='{op}'/>")
text(bx + 82, by + 30, "response", 24, anchor="middle")
for i, (v, op) in enumerate([(0.76, 0.35), (0.63, 1.0)]):
    add(f"<rect x='{bx + 180 + i * 70}' y='{by - v * bh:.1f}' width='54' height='{v * bh:.1f}' fill='{BLUE}' fill-opacity='{op}'/>")
text(bx + 242, by + 30, "knockdown", 24, anchor="middle")
text(bx + 117, by - 0.68 * bh - 10, "0.68", 24, anchor="middle")
text(bx + 207, by - 0.76 * bh - 10, "76%", 24, anchor="middle"); text(bx + 277, by - 0.63 * bh - 10, "63%", 24, anchor="middle")
add(f"<rect x='{bx + 360}' y='{y + 106}' width='22' height='22' fill='{BLUE}' fill-opacity='0.35'/>"); text(bx + 392, y + 125, "day 6", 24)
add(f"<rect x='{bx + 360}' y='{y + 140}' width='22' height='22' fill='{BLUE}'/>"); text(bx + 392, y + 159, "day 8", 24)

# ---------------------------------------------------------------- row 2: one-day harvest timecourse
y = row_frame(1, RED, "K562 timecourse", ["Drepanos, CRISPRi and knockout", "Harvested on one day"], "Fig. 3")
days = [4, 7, 10, 14]
xs = [C2 + 200 + i * 115 for i in range(4)]
text(C2, y + 108, "essential", 24, fill=RED); text(C2, y + 178, "non-essential", 24, fill=BLUE)
for x, d, a in zip(xs, days, [5, 1, 1, 0]):
    group(x, y + 100, RED, 6, dead=6 - a, cols=3, step=28, r=11)
    group(x, y + 170, BLUE, 6, cols=3, step=28, r=11)
    text(x, y + 232, f"day {d}", 24, anchor="middle")
# readout: median cells per guide relative to non-essential guides (CRISPRi), from numbers.json
bx, by, bh, bw = C3 + 50, y + 200, 120, 280
axes(bx, by, bw, bh, ylab="cells")
px = lambda d: bx + 20 + (bw - 40) * (d - 4) / 10
for col, vals in ((BLUE, [1, 1, 1, 1]), (RED, [REL[str(d)] for d in days])):
    add(f"<polyline points='" + " ".join(f"{px(d):.1f},{by - v * bh * 0.9:.1f}" for d, v in zip(days, vals))
        + f"' fill='none' stroke='{col}' stroke-width='4'/>")
    for d, v in zip(days, vals):
        add(f"<circle cx='{px(d):.1f}' cy='{by - v * bh * 0.9:.1f}' r='7' fill='{col}'/>")
for d in days:
    text(px(d), by + 30, str(d), 24, anchor="middle")
text(bx + bw + 16, by - 0.9 * bh + 8, "non-essential", 24, fill=BLUE); text(bx + bw + 16, by - 4, "essential", 24, fill=RED)

# ---------------------------------------------------------------- row 3: trimming bound
y = row_frame(2, GREEN, "Seven screens, six lines", ["Replogle: K562 (2), RPE1", "Nadig: Jurkat, HepG2", "X-Atlas/Orion: HCT116, HEK293T"], "Fig. 4, S2")
gx0, gx1, gy, gh = C2 + 10, C2 + 560, y + 170, 95
mu, sd = C2 + 150, 55
dens = lambda x: math.exp(-0.5 * ((x - mu) / sd) ** 2)
cut = mu + 0.67 * sd                       # top quarter of the control distribution (S = 0.25)
pts = [(x, gy - gh * dens(x)) for x in range(int(gx0), int(mu + 3.2 * sd), 3)]
tail = [(x, yy) for x, yy in pts if x >= cut]
add(f"<path d='M{cut:.1f},{gy} " + " ".join(f"L{x:.1f},{yy:.1f}" for x, yy in tail) + f" L{tail[-1][0]:.1f},{gy} Z' fill='{GREEN}' fill-opacity='0.35'/>")
add("<polyline points='" + " ".join(f"{x:.1f},{yy:.1f}" for x, yy in pts) + f"' fill='none' stroke='{INK}' stroke-width='3'/>")
line(gx0, gy, gx1, gy, INK, 2.5)
trim = mu + sd * 1.27                      # mean of the top quarter of a normal distribution
pert = mu + sd * 4.6
for x, col in ((mu, INK), (trim, GREEN), (pert, RED)):
    line(x, gy, x, gy - gh - 6, col, 3, "6 4")
text(mu - 6, y + 62, "control cells", 24, anchor="end")
text(trim + 14, y + 100, "top S of controls", 24, fill=GREEN)
text(pert + 10, y + 120, "perturbed", 24, fill=RED)
add(f"<path d='M{mu},{gy + 10} L{mu},{gy + 18} L{trim},{gy + 18} L{trim},{gy + 10}' fill='none' stroke='{GREEN}' stroke-width='3'/>")
text(trim + 10, gy + 26, "bound", 24, fill=GREEN)
add(f"<path d='M{mu},{gy + 40} L{mu},{gy + 48} L{pert},{gy + 48} L{pert},{gy + 40}' fill='none' stroke='{RED}' stroke-width='3'/>")
text(pert + 10, gy + 56, "observed shift", 24, fill=RED)
# readout: explainable fraction against Chronos; medians at Chronos below -2 for HCT116, K562 genome-wide and Jurkat, 0 above -0.3
bx, by, bh, bw = C3 + 50, y + 200, 120, 280
axes(bx, by, bw, bh, ylab="fraction")
text(bx + 24, by + 30, "&lt; -2", 24, anchor="middle"); text(bx + bw - 24, by + 30, "&gt; -0.3", 24, anchor="middle")
text(bx + bw / 2, by + 30, "Chronos", 24, anchor="middle")
for v in (0.38, 0.24, 0.13):
    add(f"<line x1='{bx + 24}' y1='{by - v * bh * 2:.1f}' x2='{bx + bw - 24}' y2='{by - 2:.1f}' stroke='{GREEN}' stroke-width='3'/>")
    add(f"<circle cx='{bx + 24}' cy='{by - v * bh * 2:.1f}' r='6' fill='{GREEN}'/>")
text(bx + bw + 16, by - 80, "explainable", 24); text(bx + bw + 16, by - 50, "by selective", 24); text(bx + bw + 16, by - 20, "survival", 24)

# ---------------------------------------------------------------- row 4: buffered labels
y = row_frame(3, PURPLE, "Dose-response labels", ["K562 essential genes,", "three methods"], "Fig. 5")
kx0, kx1, ky = C2 + 10, C2 + 560, y + 200
line(kx0, ky, kx1, ky, INK, 2.5, head=True)
text(kx1, ky + 32, "on-target knockdown", 24, anchor="end")
text(kx0, ky + 32, "partial", 24)
x_on = kx0 + 0.72 * (kx1 - kx0)
curve = [(x, ky - 12 - 90 * min(1, max(0, (x - x_on) / 50))) for x in range(int(kx0), int(kx1) - 20, 4)]
add("<polyline points='" + " ".join(f"{x:.1f},{yy:.1f}" for x, yy in curve) + f"' fill='none' stroke='{PURPLE}' stroke-width='4'/>")
text(kx0, y + 60, "buffered: downstream response only", 24, fill=PURPLE)
text(kx0, y + 90, "at near-complete knockdown", 24, fill=PURPLE)
for i in range(6):
    cell(kx0 + 40 + i * 58, ky - 40, PURPLE, 12)
for i in range(2):
    cell(x_on + 85 + i * 40, ky - 32, PURPLE, 12, dead=True)
text(x_on + 105, ky - 62, "removed", 22, anchor="middle", fill=PURPLE)
# readout: decrease in on-target knockdown, buffered minus other genes (estimate, 95% CI)
bx, by, bw = C3 + 120, y + 200, 260
xs_ = lambda v: bx + bw * (v + 0.05) / 0.25        # axis from -0.05 to 0.20
line(xs_(0), by - 140, xs_(0), by, GREY, 2, "5 4")
line(bx, by, bx + bw, by, INK, 2.5)
for v in (0, 0.1, 0.2):
    text(xs_(v), by + 30, f"{v:g}", 24, anchor="middle")
for i, (lab, e, lo, hi) in enumerate([("EM", 0.11, 0.04, 0.18), ("Song", 0.085, 0.02, 0.15), ("Mixscale", 0.081, 0.02, 0.15)]):
    yy = by - 115 + i * 38
    line(xs_(lo), yy, xs_(hi), yy, PURPLE, 3)
    add(f"<path d='M{xs_(e):.1f},{yy - 11} L{xs_(e) + 11:.1f},{yy} L{xs_(e):.1f},{yy + 11} L{xs_(e) - 11:.1f},{yy} Z' fill='{PURPLE}'/>")
    text(bx - 14, yy + 8, lab, 24, anchor="end")
text(bx + bw + 16, by - 80, "decrease", 24); text(bx + bw + 16, by - 50, "in on-target", 24); text(bx + bw + 16, by - 20, "knockdown", 24)

# ---------------------------------------------------------------- row 5: drug atlas
y = row_frame(4, ORANGE, "Tahoe-100M", ["Drug atlas, 50 cell lines", "per well, 24 hours"], "Fig. S1")
wx, wy = C2 + 100, y + 140
add(f"<circle cx='{wx}' cy='{wy}' r='72' fill='white' stroke='{INK}' stroke-width='3'/>")
pal = [BLUE, RED, GREEN, PURPLE, ORANGE]
for k in range(12):
    a, rr = k * 2.4, 18 + (k % 3) * 17
    cell(wx + rr * math.cos(a), wy + rr * math.sin(a), pal[k % 5], 9)
add(f"<path d='M{wx},{y + 22} C{wx - 13},{y + 42} {wx - 13},{y + 54} {wx},{y + 56} C{wx + 13},{y + 54} {wx + 13},{y + 42} {wx},{y + 22} Z' fill='{ORANGE}'/>")
text(wx + 24, y + 50, "drug", 24)
line(wx + 95, wy, wx + 170, wy, GREY, 4, head=True)
text(wx + 190, wy - 30, "cells of each line", 24)
text(wx + 190, wy + 4, "against PRISM viability", 24)
text(wx + 190, wy + 38, "of the same line and dose", 24)
# readout: flat cloud, rank correlation 0.005
bx, by, bh, bw = C3 + 50, y + 200, 120, 280
axes(bx, by, bw, bh, xlab="PRISM viability", ylab="cells")
random.seed(3)
for _ in range(70):
    add(f"<circle cx='{bx + 15 + random.random() * (bw - 30):.1f}' cy='{by - 20 - random.random() * (bh - 30):.1f}' r='4' fill='{ORANGE}' fill-opacity='0.7'/>")
text(bx + bw + 16, by - 60, "rho 0.005", 24)

y_end = Y0 + 5 * RH
for x in (C2 - 30, C3 - 30):
    line(x, 30, x, y_end - 10, LINE, 3)
add("</svg>")
Hf = y_end + 10
out = "\n".join(svg).replace("@HMM@", f"{Hf / 10:.1f}").replace("@H@", str(Hf))
os.makedirs(F, exist_ok=True)
p = os.path.join(F, "Figure1.svg")
open(p, "w").write(out + "\n")
subprocess.run(["rsvg-convert", "-w", str(round(180 / 25.4 * 600)), "-b", "white", "-o", os.path.join(F, "Figure1.png"), p], check=True)
subprocess.run(["rsvg-convert", "-f", "pdf", "-o", os.path.join(F, "Figure1.pdf"), p], check=True)
print("wrote", p)
