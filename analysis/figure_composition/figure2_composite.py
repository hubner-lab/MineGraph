#!/usr/bin/env python3
"""
Figure 2 — 4-panel composite v4 (priority layout, large PAV legend)
====================================================================
Panel priority: circular map and compartment provenance most prominent.
Gene PAV heatmap shows top-2 most variable genes per category with large legend.
IR consensus map removed (moved to supplement).

Z-path reading order (a→b→c→d):
  a  Row 1 full   : Node frequency histogram (full-width, 3:1)
  b  Row 2 left   : Consensus circular genome map (square, 1:1)  ← LARGE 6.3"
  c  Row 2 right  : Gene PAV top-2/category with category legend (1.22:1)
  d  Row 3 full   : Compartment provenance composite (full-width, 2:1)

Output:
  manuscript/figures/Figure2/Figure2_full_v4.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
})
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpecFromSubplotSpec
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE = "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
PANELS = {
    "a": f"{BASE}/manuscript/figures/Figure2/Figure2a.png",                      # histogram
    "b": f"{BASE}/manuscript/figures/Figure2/c2.png",                            # circular map (enlarged legend)
    "c": f"{BASE}/manuscript/figures/Figure2/Figure2d_gene_pav_top2cat.png",     # gene PAV top-2
    "d": f"{BASE}/manuscript/figures/Figure2/Figure2c_consensus_compartment.png",# compartment
}
OUT = f"{BASE}/manuscript/figures/Figure2/Figure2_full_v4.png"

# ── Auto-crop helper ───────────────────────────────────────────────────────────
def autocrop(img: Image.Image, pad: int = 20, threshold: int = 248) -> np.ndarray:
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
    r0 = max(0, r0 - pad);  r1 = min(arr.shape[0]-1, r1 + pad)
    c0 = max(0, c0 - pad);  c1 = min(arr.shape[1]-1, c1 + pad)
    return np.array(bg.crop((c0, r0, c1+1, r1+1)).convert("RGB"))

# ── Load and crop ──────────────────────────────────────────────────────────────
print("Loading panels …")
imgs   = {}
ratios = {}
for key, path in PANELS.items():
    raw = Image.open(path)
    arr = autocrop(raw, pad=20)
    imgs[key]   = arr
    h, w = arr.shape[:2]
    ratios[key] = w / h
    print(f"  [{key}] {raw.size[0]}×{raw.size[1]} → cropped {w}×{h}  ratio {w/h:.3f}")

# ── Layout mathematics ─────────────────────────────────────────────────────────
W = 14.0   # figure width in inches

ra, rb, rc, rd = (ratios[k] for k in "abcd")

# Row 1 (a, full-width histogram):  H1 = W / ra
# Row 2 (b circular + c PAV):       H2 = W / (rb + rc)
# Row 3 (d, full-width compartment):H3 = W / rd
H1 = W / ra
H2 = W / (rb + rc)
H3 = W / rd

gap   = 0.35
H_fig = H1 + H2 + H3 + 2*gap

print(f"\nLayout (W={W:.2f}\"): H1={H1:.2f}\"  H2={H2:.2f}\"  H3={H3:.2f}\"  total={H_fig:.2f}\"")
print(f"  Row2: b(circular)={H2*rb:.2f}\"×{H2:.2f}\"  c(PAV)={H2*rc:.2f}\"×{H2:.2f}\"")

# ── Figure ─────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(W, H_fig), facecolor="white", dpi=300)

outer = gridspec.GridSpec(
    3, 1,
    figure=fig,
    height_ratios=[H1, H2, H3],
    hspace=gap / H_fig * 3,
)

# Row 1: a (histogram, full-width)
ax_a = fig.add_subplot(outer[0])

# Row 2: b (circular, left) + c (PAV, right)
inner2 = GridSpecFromSubplotSpec(
    1, 2,
    subplot_spec=outer[1],
    width_ratios=[rb, rc],
    wspace=0.025,
)
ax_b = fig.add_subplot(inner2[0])
ax_c = fig.add_subplot(inner2[1])

# Row 3: d (compartment, full-width)
ax_d = fig.add_subplot(outer[2])

# ── Display helper ─────────────────────────────────────────────────────────────
def show(ax, arr, label):
    ax.imshow(arr, interpolation="lanczos", aspect="auto")
    ax.axis("off")
    ax.set_title(
        label,
        loc="left",
        fontsize=16,
        fontweight="bold",
        pad=4,
        fontfamily="sans-serif",
        color="black",
    )

show(ax_a, imgs["a"], "a")
show(ax_b, imgs["b"], "b")
show(ax_c, imgs["c"], "c")
show(ax_d, imgs["d"], "d")

# ── Save ───────────────────────────────────────────────────────────────────────
plt.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white", pad_inches=0.05)
print(f"\nSaved → {OUT}")
plt.close()
print("Done.")
