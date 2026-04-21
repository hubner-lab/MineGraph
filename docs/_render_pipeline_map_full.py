#!/usr/bin/env python3
"""Enhanced MineGraph metro-map.

Everything in the simple overview (`_render_pipeline_map.py`) plus:

    • CSV metadata as a required second input
    • PGGB internals: wfmash → seqwish → smoothxg → odgi
    • Stage-5 sub-stations for every stand-alone command:
        extract  : subgraph · paths · stats
        inject   : anchor subgraph · odgi untangle · gggenes · bin coverage
        convert  : GFA→FASTA · GFA→VG
        tubemap  : SequenceTubeMap server · :3210 browser

Outputs:
    docs/pipeline_metro_map_full.png
    docs/pipeline_metro_map_full.svg
"""
from __future__ import annotations

from pathlib import Path as _Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle
from matplotlib.path import Path
from matplotlib.patches import PathPatch

# ─────────────── palette ───────────────
C_GREEN  = "#00A05F"
C_BLUE   = "#14284B"
C_ORANGE = "#F7941E"
C_PINK   = "#E6007E"

C_PANEL      = "#E5E5E5"
C_PANEL_NUM  = "#D2D2D2"
C_STATION_FC = "#FFFFFF"
C_STATION_EC = "#000000"
C_LABEL      = "#111111"
C_ICON_BG    = "#FFFFFF"
C_ICON_BAR   = "#000000"

LW_LINE    = 7.0
LW_STATION = 2.2
R_STATION  = 9.0
R_STATION_SMALL = 7.5

OFFSET = {
    C_GREEN:  -9,
    C_BLUE:    0,
    C_ORANGE: +9,
    C_PINK:    0,
}

# ─────────────── canvas ───────────────
W_IN, H_IN = 25.0, 16.0
DPI = 200
W_PX, H_PX = 2500, 1600

fig, ax = plt.subplots(figsize=(W_IN, H_IN))
ax.set_xlim(0, W_PX)
ax.set_ylim(0, H_PX)
ax.invert_yaxis()
ax.set_aspect("equal")
ax.axis("off")
fig.patch.set_facecolor("white")


# ─────────────── helpers ───────────────
FONT = "DejaVu Sans"


def panel(x, y, w, h, number):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0,rounding_size=30",
        facecolor=C_PANEL, edgecolor="none", zorder=1,
    ))
    ax.text(
        x + 30, y + 18, str(number),
        fontsize=68, color=C_PANEL_NUM,
        ha="left", va="top", fontweight="bold", zorder=2,
        family=FONT,
    )


def station(x, y, label, label_dy=-30, size=11, small=False, weight="bold"):
    r = R_STATION_SMALL if small else R_STATION
    ax.add_patch(Circle(
        (x, y), r,
        facecolor=C_STATION_FC, edgecolor=C_STATION_EC, linewidth=LW_STATION,
        zorder=6,
    ))
    va = "bottom" if label_dy < 0 else "top"
    ax.text(
        x, y + label_dy, label,
        fontsize=size, color=C_LABEL,
        ha="center", va=va,
        family=FONT, fontweight=weight,
        zorder=7,
    )


def terminus(x, y, vertical=True, length=30, width=14):
    if vertical:
        ax.add_patch(FancyBboxPatch(
            (x - width / 2, y - length / 2), width, length,
            boxstyle="round,pad=0,rounding_size=3",
            facecolor="#000000", edgecolor="none", zorder=6,
        ))
    else:
        ax.add_patch(FancyBboxPatch(
            (x - length / 2, y - width / 2), length, width,
            boxstyle="round,pad=0,rounding_size=3",
            facecolor="#000000", edgecolor="none", zorder=6,
        ))


def dot(x, y, r=11):
    ax.add_patch(Circle((x, y), r, facecolor="#000000", edgecolor="none", zorder=6))


def file_icon(x, y, tag, w=62, h=80):
    fold = 15
    body = Path(
        [
            (x, y),
            (x + w - fold, y),
            (x + w, y + fold),
            (x + w, y + h),
            (x, y + h),
            (x, y),
        ],
        [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY],
    )
    ax.add_patch(PathPatch(body, facecolor=C_ICON_BG, edgecolor=C_ICON_BAR, linewidth=2.2, zorder=8))
    corner = Path(
        [(x + w - fold, y), (x + w - fold, y + fold), (x + w, y + fold)],
        [Path.MOVETO, Path.LINETO, Path.LINETO],
    )
    ax.add_patch(PathPatch(corner, facecolor="none", edgecolor=C_ICON_BAR, linewidth=2.2, zorder=9))
    bar_h = 20
    bar_y = y + h - bar_h - 9
    ax.add_patch(Rectangle(
        (x + 4, bar_y), w - 8, bar_h,
        facecolor=C_ICON_BAR, edgecolor="none", zorder=9,
    ))
    ax.text(
        x + w / 2, bar_y + bar_h / 2, tag,
        fontsize=10, color="#FFFFFF", ha="center", va="center",
        fontweight="bold", zorder=10, family=FONT,
    )


def line(color, points, offset=None, zorder=3, width=None):
    off = OFFSET[color] if offset is None else offset
    shifted = [(x, y + off) for (x, y) in points]
    R = 26
    verts = [shifted[0]]
    codes = [Path.MOVETO]
    for i in range(1, len(shifted) - 1):
        prev, cur, nxt = shifted[i - 1], shifted[i], shifted[i + 1]
        from_prev = (cur[0] - prev[0], cur[1] - prev[1])
        to_next = (nxt[0] - cur[0], nxt[1] - cur[1])
        lp = max((from_prev[0] ** 2 + from_prev[1] ** 2) ** 0.5, 1e-6)
        ln = max((to_next[0] ** 2 + to_next[1] ** 2) ** 0.5, 1e-6)
        up = (from_prev[0] / lp, from_prev[1] / lp)
        un = (to_next[0] / ln, to_next[1] / ln)
        if up[0] * un[0] + up[1] * un[1] > 0.995:
            verts.append(cur)
            codes.append(Path.LINETO)
            continue
        r = min(R, lp * 0.45, ln * 0.45)
        p1 = (cur[0] - up[0] * r, cur[1] - up[1] * r)
        p2 = (cur[0] + un[0] * r, cur[1] + un[1] * r)
        verts.append(p1)
        codes.append(Path.LINETO)
        verts.extend([cur, p2])
        codes.extend([Path.CURVE3, Path.CURVE3])
    verts.append(shifted[-1])
    codes.append(Path.LINETO)
    path = Path(verts, codes)
    ax.add_patch(PathPatch(
        path, facecolor="none", edgecolor=color,
        linewidth=LW_LINE if width is None else width,
        capstyle="round", joinstyle="round", zorder=zorder,
    ))


# ─────────────── panels ───────────────
panel(  60, 340, 1120, 320, 1)
panel(1200, 340,  620, 320, 2)
panel(  60, 760, 1760, 320, 3)
panel(1840, 340,  540, 320, 4)
panel(  60, 1160, 2320, 400, 5)


# ─────────────── key y-rows ───────────────
y1 = 520         # Stage 1 & 2 main spine
y3 = 920         # Stage 3 trunk
y4 = 520         # Stage 4 main (same row as 1/2)
y5_hub = 1300    # Stage 5 hub row (extract/inject/convert/tubemap)
y5_sub = 1440    # Stage 5 sub-station row


# ─────────────── Stage 1 stations ───────────────
S1 = [
    ( 210, y1, "cat fastq",         -34),
    ( 350, y1, "prepare_and\n_mash_input", +32),
    ( 510, y1, "mash distance\nauto-tune", -34),
    ( 690, y1, "RepeatMasker",      +32),
    ( 870, y1, "run_repeat\nmask",   -34),
    (1050, y1, "params.yaml",        +32),
]

# Stage 2 — PGGB sub-stations
pggb_x = 1260
wfm_x = 1400
seq_x = 1500
smo_x = 1600
odg_x = 1700

# Stage 3 stations
S3 = [
    ( 180, y3, "ODGI stats\n(5a)",          -36),
    ( 340, y3, "node\nhistogram",           +32),
    ( 500, y3, "path\nheatmap",             -36),
    ( 660, y3, "similarity\nwindow",        +32),
    ( 820, y3, "Cloud / Shell\n/ Core",     -36),
    ( 980, y3, "consensus\nFASTA",          +32),
    (1160, y3, "ODGI PAV\nmatrix",          -36),
    (1360, y3, "Jaccard\ndistances",        +32),
    (1540, y3, "top-N\ninteractive",        -36),
]

# Stage 4 stations
up_x,  up_y  = 1920, y4
fel_x, fel_y = 2060, y4 - 90
uf_x,  uf_y  = 2060, y4 + 90
rax_x, rax_y = 2220, y4
tree_icon_x  = 2298
tree_term_x  = 2330

# Stage 5 hubs (bigger, main)
H5 = {
    "extract": ( 430, y5_hub),
    "inject":  (1040, y5_hub),
    "convert": (1600, y5_hub),
    "tubemap": (2020, y5_hub),
}

# Stage 5 sub-stations (smaller)
SUB_EXTRACT = [
    (310, y5_sub, "subgraph",  -32, 10),
    (430, y5_sub, "paths",     +30, 10),
    (550, y5_sub, "stats",     -32, 10),
]
SUB_INJECT = [
    ( 820, y5_sub, "anchor\nsubgraph", +30, 10),
    ( 950, y5_sub, "odgi\nuntangle",   -34, 10),
    (1080, y5_sub, "gggenes",          +30, 10),
    (1210, y5_sub, "bin\ncoverage",    -34, 10),
]
SUB_CONVERT = [
    (1520, y5_sub, "GFA → FASTA", -32, 10),
    (1700, y5_sub, "GFA → VG",    +30, 10),
]
SUB_TUBEMAP = [
    (1950, y5_sub, "TubeMap\nserver",  -34, 10),
    (2110, y5_sub, "http://\nlocalhost:3210", +30, 9.5),
]


# ─────────────── metro lines ───────────────
trunk_s1 = [(140, y1)] + [(x, y) for (x, y, *_ ) in S1]
# Stage 2 sub-spine: PGGB → wfmash → seqwish → smoothxg → odgi
trunk_s2 = [(pggb_x, y1), (wfm_x, y1), (seq_x, y1), (smo_x, y1), (odg_x, y1)]
trunk_12 = trunk_s1 + trunk_s2

# drop from odgi → stage 3 (enters from top-right of stage 3)
drop_to_s3 = [
    (odg_x, y1),
    (odg_x, 620),
    (1680, 620),
    (1680, y3 - 60),
    (S3[-1][0] + 80, y3 - 60),
    (S3[-1][0] + 80, y3),
    (S3[-1][0], y3),
]
s3_trunk_reverse = [(x, y) for (x, y, *_) in reversed(S3)]

# From stage 3 leftmost, climb up to UPGMA in stage 4 via a dogleg
s3_to_s4 = [
    (S3[0][0], y3),
    (S3[0][0] - 60, y3),
    (S3[0][0] - 60, 600),
    (odg_x + 80, 600),
    (odg_x + 80, y4),
    (up_x, up_y),
]

# Green — default (graph-UPGMA · Felsenstein)
green_path = trunk_12 + drop_to_s3[1:] + s3_trunk_reverse[1:] + s3_to_s4[1:] + [
    (up_x + 60, up_y),
    (fel_x - 40, up_y),
    (fel_x, fel_y),
    (rax_x - 60, fel_y),
    (rax_x - 60, rax_y),
    (rax_x + 80, rax_y),
    (tree_term_x, rax_y),
]

# Blue — UFBoot-PAV
blue_path = trunk_12 + drop_to_s3[1:] + s3_trunk_reverse[1:] + s3_to_s4[1:] + [
    (up_x + 60, up_y),
    (uf_x - 40, up_y),
    (uf_x, uf_y),
    (rax_x - 60, uf_y),
    (rax_x - 60, rax_y),
    (rax_x + 80, rax_y),
    (tree_term_x, rax_y),
]

# Orange — MSA / RAxML-NG (skips stage 3)
orange_path = trunk_12 + [
    (odg_x, y1),
    (rax_x, rax_y),
    (rax_x + 80, rax_y),
    (tree_term_x, rax_y),
]

line(C_GREEN,  green_path,  zorder=3)
line(C_BLUE,   blue_path,   zorder=3)
line(C_ORANGE, orange_path, zorder=3)


# Pink — stand-alone tools
# main hub spine
pink_hub_path = [
    (140, y5_hub),
    (H5["extract"][0], y5_hub),
    (H5["inject"][0], y5_hub),
    (H5["convert"][0], y5_hub),
    (H5["tubemap"][0], y5_hub),
    (2200, y5_hub),
]
line(C_PINK, pink_hub_path, zorder=4)

# spurs for each tool — drop down and traverse sub-stations
def pink_spur(hub_xy, subs):
    hx, hy = hub_xy
    first_x = subs[0][0]
    last_x  = subs[-1][0]
    path = [(hx, hy), (hx, y5_sub - 50)]
    if hx > first_x:
        path += [(first_x, y5_sub - 50)]
    elif hx < last_x:
        path += [(last_x, y5_sub - 50)]
    path += [(first_x, y5_sub), (last_x, y5_sub)]
    line(C_PINK, path, zorder=4, width=5.5)

pink_spur(H5["extract"], SUB_EXTRACT)
pink_spur(H5["inject"],  SUB_INJECT)
pink_spur(H5["convert"], SUB_CONVERT)
pink_spur(H5["tubemap"], SUB_TUBEMAP)


# ─────────────── origin + termini + file icons ───────────────
# Stage 1 inputs — FASTA + CSV merge
dot(140, y1)
# vertical coupling connecting the two input icons to the trunk
ax.add_patch(FancyBboxPatch(
    (132, y1 - 80), 16, 160,
    boxstyle="round,pad=0,rounding_size=4",
    facecolor="#000000", edgecolor="none", zorder=5,
))
file_icon(40, y1 - 110, "FASTA")
file_icon(40, y1 + 30,  "CSV")
# small label "required inputs"
ax.text(72, y1 - 128, "required inputs",
        fontsize=10, color="#222", ha="center", va="bottom",
        style="italic", family=FONT, fontweight="bold")

# PGGB outputs: GFA / OG / VCF / HTML (above panel 2, between header and panel)
terminus(odg_x, 300, vertical=False, length=400, width=12)
for i, tag in enumerate(["GFA", "OG", "VCF", "HTML"]):
    file_icon(odg_x - 200 + i * 100, 210, tag)
# short connector from PGGB spine (odgi) up to the terminus bar
ax.add_patch(Rectangle((odg_x - 6, 300), 12, y1 - 300,
                       facecolor="#000000", edgecolor="none", zorder=5))

# Stage 3 outputs — XLSX / PNG×3 / FASTA / TSV / HTML (interactive)
s3_icon_xs = [
    ( 180, "XLSX"),
    ( 340, "PNG"),
    ( 500, "PNG"),
    ( 660, "PNG"),
    ( 980, "FASTA"),
    (1160, "TSV"),
    (1540, "HTML"),
]
terminus((s3_icon_xs[0][0] + s3_icon_xs[-1][0]) / 2, 750,
         vertical=False, length=s3_icon_xs[-1][0] - s3_icon_xs[0][0] + 80, width=12)
for (xc, tag) in s3_icon_xs:
    file_icon(xc - 31, 660, tag)

# Stage 4 tree outputs (shifted inward so they stay inside the 2400 px canvas)
file_icon(tree_icon_x, rax_y - 40, "NWK")
file_icon(tree_icon_x, rax_y + 50, "PNG")
file_icon(tree_icon_x, rax_y - 130, "TSV")
terminus(tree_term_x, rax_y, vertical=True, length=230, width=12)

# Stage 5 inputs — pre-built graph origin + BED + representatives
dot(140, y5_hub)
file_icon(50, y5_hub - 40, "GFA")

# BED + representatives for inject
file_icon(H5["inject"][0] - 80, y5_hub - 130, "BED")
file_icon(H5["inject"][0] + 20, y5_hub - 130, "TXT")
ax.text(H5["inject"][0] - 30, y5_hub - 150,
        "inject inputs",
        fontsize=10, color="#222", ha="center", va="bottom",
        style="italic", family=FONT, fontweight="bold")

# Stage 5 outputs — below sub-station row
# extract outputs
file_icon(SUB_EXTRACT[0][0] - 31, y5_sub + 50, "GFA")
file_icon(SUB_EXTRACT[2][0] - 31, y5_sub + 50, "TSV")

# inject outputs
file_icon(SUB_INJECT[2][0] - 31, y5_sub + 50, "PNG")
file_icon(SUB_INJECT[3][0] - 31, y5_sub + 50, "PNG")

# convert outputs
file_icon(SUB_CONVERT[0][0] - 31, y5_sub + 50, "FASTA")
file_icon(SUB_CONVERT[1][0] - 31, y5_sub + 50, "VG")

# tubemap output
file_icon(SUB_TUBEMAP[1][0] - 31, y5_sub + 50, "HTML")


# ─────────────── stations ───────────────
for (x, y, label, dy) in S1:
    station(x, y, label, label_dy=dy)

# Stage 2 PGGB + internals (PGGB is a "major" station, internals are smaller)
station(pggb_x, y1, "PGGB", label_dy=-32, size=14)
station(wfm_x, y1, "wfmash",  label_dy=+30, size=10, small=True)
station(seq_x, y1, "seqwish", label_dy=-28, size=10, small=True)
station(smo_x, y1, "smoothxg",label_dy=+30, size=10, small=True)
station(odg_x, y1, "odgi",    label_dy=-28, size=10, small=True)

for (x, y, label, dy) in S3:
    station(x, y, label, label_dy=dy, size=10.5)

station(up_x,  up_y,  "UPGMA",               label_dy=-32, size=13)
station(fel_x, fel_y, "Felsenstein\n(default)", label_dy=-32, size=10.5)
station(uf_x,  uf_y,  "UFBoot-PAV",          label_dy=+32, size=10.5)
station(rax_x, rax_y, "RAxML-NG\n(MSA)",     label_dy=-32, size=10.5)

# Stage 5 hubs
for name, (x, y) in H5.items():
    station(x, y, name, label_dy=-32, size=14)

# Stage 5 sub-stations
for group in (SUB_EXTRACT, SUB_INJECT, SUB_CONVERT, SUB_TUBEMAP):
    for (x, y, label, dy, sz) in group:
        station(x, y, label, label_dy=dy, size=sz, small=True)


# ─────────────── title + legend (top header) ───────────────
ax.text(
    60, 150, "MineGraph",
    fontsize=64, color="#111111",
    ha="left", va="bottom", fontweight="bold",
    family=FONT,
)
ax.text(
    60, 198, "plastid & mitochondrial graph-pangenome pipeline  —  full workflow",
    fontsize=17, color="#444444", style="italic",
    ha="left", va="bottom",
    family=FONT, fontweight="bold",
)

# STAGE legend (top-right, left column)
lx = 1280
ax.text(lx, 70, "STAGE", fontsize=14, fontweight="bold", color="#111111",
        family=FONT, ha="left", va="bottom")
for i, s in enumerate([
    "1.  Pre-processing (merge · mash · RepeatMask · auto-tune)",
    "2.  Graph construction (PGGB: wfmash → seqwish → smoothxg → odgi)",
    "3.  Statistics & PAV (ODGI · consensus · Jaccard)",
    "4.  Phylogenetics (UPGMA · RAxML-NG)",
    "5.  Stand-alone tools (extract · inject · convert · tubemap)",
]):
    ax.text(lx, 92 + i * 20, s, fontsize=11, color="#222",
            family=FONT, ha="left", va="bottom", fontweight="bold")

# METHOD legend (top-right, right column)
lx2 = 1880
ax.text(lx2, 70, "METHOD", fontsize=14, fontweight="bold", color="#111111",
        family=FONT, ha="left", va="bottom")
for i, (col, s) in enumerate([
    (C_GREEN,  "Default — graph-UPGMA · Felsenstein bootstrap"),
    (C_BLUE,   "Graph-UPGMA · UFBoot-PAV adaptive bootstrap"),
    (C_ORANGE, "MSA — RAxML-NG  (--tree_type msa)"),
    (C_PINK,   "Stand-alone tools (on a pre-built graph)"),
]):
    yy = 92 + i * 20
    ax.add_patch(Rectangle((lx2, yy - 10), 46, 7, facecolor=col, edgecolor="none"))
    ax.text(lx2 + 58, yy, s, fontsize=11, color="#222",
            family=FONT, ha="left", va="bottom", fontweight="bold")


# ─────────────── save ───────────────
out_dir = _Path(__file__).resolve().parent
png_path = out_dir / "pipeline_metro_map_full.png"
svg_path = out_dir / "pipeline_metro_map_full.svg"
fig.tight_layout(pad=0)
fig.savefig(png_path, dpi=DPI, bbox_inches="tight", facecolor="white", pad_inches=0.15)
fig.savefig(svg_path, bbox_inches="tight", facecolor="white", pad_inches=0.15)
print(f"wrote {png_path}")
print(f"wrote {svg_path}")
