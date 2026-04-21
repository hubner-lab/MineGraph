#!/usr/bin/env python3
"""
Figure 2 — 5-panel composite v3 (priority layout)
==================================================
Panel priority: circular map + compartment provenance most prominent (top row).
Gene PAV heatmap placed last (bottom-right, lowest visual weight).

Z-path reading order (a→b→c→d→e):
  a  Row 1 left    : Consensus circular genome map  (square, 1:1)
  b  Row 1 right   : Compartment provenance composite (wide, 2:1)
  c  Row 2 full    : Node frequency histogram (full-width, 3:1)
  d  Row 3 left    : IR consensus position map (wide, 1.5:1)
  e  Row 3 right   : Variant gene PAV heatmap (portrait, 0.557:1)  ← LAST

Content → position mapping vs v2:
  v2-a (histogram)   → v3-c
  v2-b (circular)    → v3-a  (elevated: top-left, most prominent)
  v2-c (compartment) → v3-b  (elevated: top-right, large)
  v2-d (gene PAV)    → v3-e  (demoted: last)
  v2-e (IR map)      → v3-d

Output:
  manuscript/figures/Figure2/Figure2_full_v3.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpecFromSubplotSpec
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE = "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
PANELS = {
    "a": f"{BASE}/manuscript/figures/Figure2/Figure2b_extracted.png",  # circular map
    "b": f"{BASE}/manuscript/figures/Figure2/Figure2c_consensus_compartment.png",  # compartment
    "c": f"{BASE}/manuscript/figures/Figure2/Figure2a.png",            # histogram
    "d": f"{BASE}/manuscript/figures/FigSX_ir_consensus_map.png",      # IR map
    "e": f"{BASE}/manuscript/figures/Figure2/Figure2d_gene_pav_compact.png",  # gene PAV
}
OUT = f"{BASE}/manuscript/figures/Figure2/Figure2_full_v3.png"

# ── Auto-crop helper ───────────────────────────────────────────────────────────
def autocrop(img: Image.Image, pad: int = 22, threshold: int = 248) -> np.ndarray:
    rgba = img.convert("RGBA")
    bg   = Image.new("RGBA", rgba.size, (255, 255, 255, 255))
    bg.paste(rgba, mask=rgba.split()[3])
    arr  = np.array(bg.convert("L"))
    rows = np.any(arr < threshold, axis=1)
    cols = np.any(arr < threshold, axis=0)
    if not rows.any():
        return np.array(bg.convert("RGB"))
    r0, r1 = np.where(rows)[0][[0, -1]]
    c0, c1 = np.where(cols)[0][[0, -1]]
    r0 = max(0, r0 - pad);  r1 = min(arr.shape[0] - 1, r1 + pad)
    c0 = max(0, c0 - pad);  c1 = min(arr.shape[1] - 1, c1 + pad)
    return np.array(bg.crop((c0, r0, c1 + 1, r1 + 1)).convert("RGB"))

# ── Load and crop all panels ───────────────────────────────────────────────────
print("Loading panels …")
imgs   = {}
ratios = {}
for key, path in PANELS.items():
    raw  = Image.open(path)
    arr  = autocrop(raw, pad=20)
    imgs[key]   = arr
    h, w = arr.shape[:2]
    ratios[key] = w / h
    print(f"  [{key}] {raw.size[0]}×{raw.size[1]} → cropped {w}×{h}  ratio {w/h:.3f}")

# ── Layout mathematics ─────────────────────────────────────────────────────────
# W = 14 inches (356 mm — wide double-column)
W = 14.0

ra, rb, rc, rd, re = (ratios[k] for k in "abcde")

# Row 1 (a + b):  H1 = W / (ra + rb)
# Row 2 (c full): H2 = W / rc
# Row 3 (d + e):  H3 = W / (rd + re)
H1 = W / (ra + rb)
H2 = W / rc
H3 = W / (rd + re)
H_total = H1 + H2 + H3

gap   = 0.12          # 3 mm per gap, 2 gaps
H_fig = H_total + 2 * gap

print(f"\nLayout (W={W:.2f}\"):  H1={H1:.2f}\"  H2={H2:.2f}\"  H3={H3:.2f}\"  → total {H_fig:.2f}\"")
print(f"  Row1 widths:  a={H1*ra:.2f}\"  b={H1*rb:.2f}\"")
print(f"  Row3 widths:  d={H3*rd:.2f}\"  e={H3*re:.2f}\"")

# ── Figure ─────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(W, H_fig), facecolor="white", dpi=300)

outer = gridspec.GridSpec(
    3, 1,
    figure=fig,
    height_ratios=[H1, H2, H3],
    hspace=gap / H_fig * 3,
)

# Row 1: a (circular, left) + b (compartment, right)
inner1 = GridSpecFromSubplotSpec(
    1, 2,
    subplot_spec=outer[0],
    width_ratios=[ra, rb],
    wspace=0.025,
)
ax_a = fig.add_subplot(inner1[0])
ax_b = fig.add_subplot(inner1[1])

# Row 2: c (histogram, full-width)
ax_c = fig.add_subplot(outer[1])

# Row 3: d (IR map, left) + e (gene PAV, right)
inner3 = GridSpecFromSubplotSpec(
    1, 2,
    subplot_spec=outer[2],
    width_ratios=[rd, re],
    wspace=0.025,
)
ax_d = fig.add_subplot(inner3[0])
ax_e = fig.add_subplot(inner3[1])

# ── Display helper ─────────────────────────────────────────────────────────────
def show(ax, arr, label):
    ax.imshow(arr, interpolation="lanczos", aspect="auto")
    ax.axis("off")
    ax.set_title(
        label,
        loc="left",
        fontsize=11,
        fontweight="bold",
        pad=3,
        fontfamily="sans-serif",
        color="black",
    )

show(ax_a, imgs["a"], "a")
show(ax_b, imgs["b"], "b")
show(ax_c, imgs["c"], "c")
show(ax_d, imgs["d"], "d")
show(ax_e, imgs["e"], "e")

# ── Save ───────────────────────────────────────────────────────────────────────
plt.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white", pad_inches=0.05)
print(f"\nSaved → {OUT}")
plt.close()
print("Done.")
