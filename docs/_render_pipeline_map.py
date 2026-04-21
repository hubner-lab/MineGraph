#!/usr/bin/env python3
"""Render the MineGraph pipeline as an nf-core style metro-map.

Produces:
    docs/pipeline_metro_map.png
    docs/pipeline_metro_map.svg

Run:
    python docs/_render_pipeline_map.py
"""
from __future__ import annotations

from pathlib import Path as _Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle
from matplotlib.path import Path
from matplotlib.patches import PathPatch

# ─────────────── palette ───────────────
C_GREEN  = "#00A05F"   # default — graph-UPGMA · Felsenstein bootstrap
C_BLUE   = "#14284B"   # graph-UPGMA · UFBoot-PAV bootstrap
C_ORANGE = "#F7941E"   # MSA — RAxML-NG
C_PINK   = "#E6007E"   # stand-alone commands (extract / inject / convert / tubemap)

C_PANEL      = "#E5E5E5"
C_PANEL_NUM  = "#FFFFFF"
C_STATION_FC = "#FFFFFF"
C_STATION_EC = "#000000"
C_LABEL      = "#111111"
C_ICON_BG    = "#FFFFFF"
C_ICON_BAR   = "#000000"

LW_LINE      = 7.0
LW_STATION   = 2.2
R_STATION    = 9.0

OFFSET = {
    C_GREEN:  -9,
    C_BLUE:    0,
    C_ORANGE: +9,
    C_PINK:    0,
}

# ─────────────── canvas ───────────────
W_IN, H_IN = 22.0, 12.5
DPI = 200
W_PX, H_PX = 2200, 1250

fig, ax = plt.subplots(figsize=(W_IN, H_IN))
ax.set_xlim(0, W_PX)
ax.set_ylim(0, H_PX)
ax.invert_yaxis()
ax.set_aspect("equal")
ax.axis("off")
fig.patch.set_facecolor("white")


# ─────────────── helpers ───────────────
def panel(x, y, w, h, number):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0,rounding_size=30",
        facecolor=C_PANEL, edgecolor="none", zorder=1,
    ))
    ax.text(
        x + 28, y + 14, str(number),
        fontsize=52, color=C_PANEL_NUM,
        ha="left", va="top", fontweight="bold", zorder=2,
        family="DejaVu Sans",
    )


def station(x, y, label, label_dy=-30, size=11):
    ax.add_patch(Circle(
        (x, y), R_STATION,
        facecolor=C_STATION_FC, edgecolor=C_STATION_EC, linewidth=LW_STATION,
        zorder=6,
    ))
    va = "bottom" if label_dy < 0 else "top"
    ax.text(
        x, y + label_dy, label,
        fontsize=size, color=C_LABEL,
        ha="center", va=va,
        family="DejaVu Sans",
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
    """Document icon (folded corner) with a black tag bar near the bottom."""
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
        fontweight="bold", zorder=10, family="DejaVu Sans",
    )


def line(color, points, offset=None, zorder=3):
    """Smoothed poly-line routed through (x,y) waypoints with an orthogonal offset."""
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
        path, facecolor="none", edgecolor=color, linewidth=LW_LINE,
        capstyle="round", joinstyle="round", zorder=zorder,
    ))


# ─────────────── stages ───────────────
# Stage 1 — Preprocessing
panel(  60, 220, 1060, 300, 1)
# Stage 2 — Graph construction
panel(1140, 220,  400, 300, 2)
# Stage 3 — Statistics & PAV (wraps under 1+2)
panel(  60, 620, 1480, 300, 3)
# Stage 4 — Phylogenetics
panel(1560, 220,  580, 300, 4)
# Stage 5 — Stand-alone tools
panel(1560, 620,  580, 300, 5)

# ─────────────── stations / coordinates ───────────────
y1 = 380
# Stage 1 stations
S1 = [
    ( 180, y1, "cat fastq",        -34),
    ( 320, y1, "prepare_and\n_mash_input", +32),
    ( 480, y1, "mash distance\nauto-tune",  -34),
    ( 660, y1, "RepeatMasker",              +32),
    ( 840, y1, "run_repeat\nmask",          -34),
    (1020, y1, "params.yaml",               +32),
]

# Stage 2 — single station: PGGB
pggb_x, pggb_y = 1340, y1

# Stage 3 stations (wraps under)
y3 = 780
S3 = [
    ( 180, y3, "ODGI stats\n(5a)",          -36),
    ( 340, y3, "node\nhistogram",           +32),
    ( 500, y3, "path\nheatmap",             -36),
    ( 660, y3, "similarity\nwindow",        +32),
    ( 820, y3, "Cloud / Shell\n/ Core",     -36),
    ( 980, y3, "consensus\nFASTA",          +32),
    (1160, y3, "ODGI PAV\nmatrix",          -36),
    (1340, y3, "Jaccard\ndistances",        +32),
]

# Stage 4 stations — phylogenetics
y4 = 380
up_x,  up_y  = 1680, y4
fel_x, fel_y = 1820, y4 - 86
uf_x,  uf_y  = 1820, y4 + 86
rax_x, rax_y = 1980, y4

# Stage 5 stations — stand-alone tools
y5 = 780
ext_x, inj_x, conv_x, tube_x = 1680, 1800, 1920, 2060


# ─────────────── metro lines ───────────────
# Shared trunk: FASTA → ... → PGGB
trunk = [(100, y1)] + [(x, y) for (x, y, *_ ) in S1] + [(pggb_x, pggb_y)]

# Stage 2 → Stage 3 bridge (bends down and left to enter Stage 3 trunk at its left end)
drop_to_s3 = [
    (pggb_x, pggb_y),
    (pggb_x - 40, pggb_y),
    (pggb_x - 40, 560),
    (160, 560),
    (160, y3),
    (180, y3),
]

# Stage 3 trunk: left → right
s3_trunk = [(x, y) for (x, y, *_) in S3]

# Stage 3 → Stage 4 bridge
s3_to_s4 = [
    (1340, y3),
    (1440, y3),
    (1440, 560),
    (1620, 560),
    (1620, y4),
    (up_x, up_y),
]

# Green path (Felsenstein default)
green_path = trunk + drop_to_s3[1:] + s3_trunk[1:] + s3_to_s4[1:] + [
    (up_x + 70, up_y),
    (fel_x - 40, up_y),
    (fel_x, fel_y),
    (rax_x - 60, fel_y),
    (rax_x - 60, rax_y),
    (rax_x + 70, rax_y),
    (2140, rax_y),
]

# Blue path (UFBoot-PAV)
blue_path = trunk + drop_to_s3[1:] + s3_trunk[1:] + s3_to_s4[1:] + [
    (up_x + 70, up_y),
    (uf_x - 40, up_y),
    (uf_x, uf_y),
    (rax_x - 60, uf_y),
    (rax_x - 60, rax_y),
    (rax_x + 70, rax_y),
    (2140, rax_y),
]

# Orange path (MSA / RAxML-NG) — skips Stage 3, goes straight through the inter-panel gap to RAxML
orange_path = trunk + [
    (pggb_x, pggb_y),
    (rax_x, rax_y),
    (rax_x + 70, rax_y),
    (2140, rax_y),
]

# Pink path (stand-alone tools)
pink_path = [
    (1580, 680),          # enters Stage 5 from the left edge, above the tools row
    (1580, y5),
    (ext_x, y5),
    (inj_x, y5),
    (conv_x, y5),
    (tube_x, y5),
    (2140, y5),
]

line(C_GREEN,  green_path,  zorder=3)
line(C_BLUE,   blue_path,   zorder=3)
line(C_ORANGE, orange_path, zorder=3)
line(C_PINK,   pink_path,   zorder=4)


# ─────────────── origin + termini ───────────────
# FASTA input → spine origin
dot(100, y1)
file_icon(28, y1 - 40, "FASTA")

# PGGB outputs: GFA / OG / VCF / HTML — icons above the Stage-2 panel
stage2_icons_y = 120
terminus(pggb_x, 210, vertical=False, length=340, width=12)
for i, tag in enumerate(["GFA", "OG", "VCF", "HTML"]):
    file_icon(pggb_x - 170 + i * 90, stage2_icons_y, tag)

# Stage 3 terminus — XLSX / PNG / FASTA / TSV icons placed above their source stations
stage3_icons_y = 520
s3_icon_xs = [
    (180, "XLSX"),    # ODGI stats 5a
    (340, "PNG"),     # node histogram 5b
    (500, "PNG"),     # path heatmap 5b
    (660, "PNG"),     # similarity window 5b
    (980, "FASTA"),   # consensus FASTA 5c
    (1160, "TSV"),    # ODGI PAV matrix 5c
]
# terminus spans the covered stations
terminus((s3_icon_xs[0][0] + s3_icon_xs[-1][0]) / 2, 612,
         vertical=False, length=s3_icon_xs[-1][0] - s3_icon_xs[0][0] + 80, width=12)
for (xc, tag) in s3_icon_xs:
    file_icon(xc - 31, stage3_icons_y, tag)

# Stage 4 tree outputs — to the right of the rax station
file_icon(2140, rax_y - 40, "NWK")
file_icon(2140, rax_y + 50, "PNG")
terminus(2140, rax_y, vertical=True, length=110, width=12)

# Stage 5 — pre-built graph input, BED input, tool outputs
file_icon(1480, 660, "GFA")    # pre-built graph input to stage 5
dot(1580, 680)
file_icon(inj_x - 32, y5 - 112, "BED")    # BED input above 'inject'
file_icon(conv_x - 60, y5 + 60, "FASTA")  # convert FASTA output
file_icon(conv_x + 10, y5 + 60, "VG")     # convert VG output
file_icon(2140, y5 - 40, "HTML")          # tubemap browser endpoint


# ─────────────── stations ───────────────
for (x, y, label, dy) in S1:
    station(x, y, label, label_dy=dy)
station(pggb_x, pggb_y, "PGGB", label_dy=-30, size=14)
for (x, y, label, dy) in S3:
    station(x, y, label, label_dy=dy, size=10.5)
station(up_x,  up_y,  "UPGMA",               label_dy=-30, size=13)
station(fel_x, fel_y, "Felsenstein\n(default)", label_dy=-30, size=10.5)
station(uf_x,  uf_y,  "UFBoot-PAV",          label_dy=+32, size=10.5)
station(rax_x, rax_y, "RAxML-NG\n(MSA)",     label_dy=-30, size=10.5)
station(ext_x,  y5, "extract",  label_dy=-30, size=12)
station(inj_x,  y5, "inject",   label_dy=+32, size=12)
station(conv_x, y5, "convert",  label_dy=-30, size=12)
station(tube_x, y5, "tubemap",  label_dy=+32, size=12)


# ─────────────── title + logo + legend ───────────────
ax.text(
    80, 1010, "MineGraph",
    fontsize=58, color="#111111",
    ha="left", va="top", fontweight="bold",
    family="DejaVu Sans",
)
ax.text(
    80, 1082, "plastid · mitochondrial graph-pangenome pipeline",
    fontsize=17, color="#444444", style="italic",
    ha="left", va="top", family="DejaVu Sans",
)

# STAGE list
ax.text(80,  1140, "STAGE", fontsize=13, fontweight="bold", color="#111111", family="DejaVu Sans")
for i, s in enumerate([
    "1. Pre-processing (merge, mash, RepeatMask, auto-tune)",
    "2. Graph construction (PGGB)",
    "3. Statistics & PAV (ODGI, consensus, Jaccard)",
    "4. Phylogenetics (UPGMA · RAxML-NG)",
    "5. Stand-alone tools (extract · inject · convert · tubemap)",
]):
    ax.text(80, 1164 + i * 17, s, fontsize=10.5, color="#222222", family="DejaVu Sans")

# METHOD legend
ax.text(860, 1140, "METHOD", fontsize=13, fontweight="bold", color="#111111", family="DejaVu Sans")
for i, (col, s) in enumerate([
    (C_GREEN,  "Default — graph-UPGMA · Felsenstein bootstrap"),
    (C_BLUE,   "Graph-UPGMA · UFBoot-PAV adaptive bootstrap"),
    (C_ORANGE, "MSA — RAxML-NG (--tree_type msa)"),
    (C_PINK,   "Stand-alone tools (operate on a pre-built graph)"),
]):
    yy = 1164 + i * 17
    ax.add_patch(Rectangle((860, yy - 4), 46, 6, facecolor=col, edgecolor="none"))
    ax.text(918, yy, s, fontsize=10.5, color="#222222", va="center", family="DejaVu Sans")


# ─────────────── save ───────────────
out_dir = _Path(__file__).resolve().parent
png_path = out_dir / "pipeline_metro_map.png"
svg_path = out_dir / "pipeline_metro_map.svg"

fig.tight_layout(pad=0)
fig.savefig(png_path, dpi=DPI, bbox_inches="tight", facecolor="white", pad_inches=0.15)
fig.savefig(svg_path, bbox_inches="tight", facecolor="white", pad_inches=0.15)
print(f"wrote {png_path}")
print(f"wrote {svg_path}")
