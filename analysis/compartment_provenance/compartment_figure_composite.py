#!/usr/bin/env python3
"""
Figure 2c — Consensus vs Pangenome BLAST composite panel
=========================================================
Three-component panel:
  LEFT   : Auto-cropped Bandage visualization (192_vs_consensus.png)
            + BLAST identity rainbow color bar annotation
  TOP-R  : Donut chart — overall compartment composition of BLAST hits
  BOT-R  : Horizontal stacked bars per family (FAMILY_ORDER for consistency
            with Figure 2d); Restionaceae anomaly flagged

Output:
  manuscript/figures/Figure2/Figure2c_consensus_compartment.png
  rebuttal/exports/R2_2/Figure2c_consensus_compartment.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colorbar import ColorbarBase
from matplotlib.cm import ScalarMappable
from matplotlib.patches import FancyArrowPatch
from PIL import Image

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE    = "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
BANDAGE = f"{BASE}/manuscript/figures/Figure2/192_vs_consensus.png"
FREQ    = f"{BASE}/rebuttal/exports/R2_2/per_family_compartment.tsv"
SUMM    = f"{BASE}/rebuttal/exports/R2_2/compartment_coverage_summary.tsv"
OUTDIR  = f"{BASE}/rebuttal/exports/R2_2"
FIG2C   = f"{BASE}/manuscript/figures/Figure2/Figure2c_consensus_compartment.png"

# ── Compartment color palette (consistent with existing R2.2 figures) ──────────
COMP_COLORS = {
    "LSC": "#4a86c8",
    "IRb": "#e05c5c",
    "SSC": "#3a9e5f",
    "IRa": "#d4b33a",
}
COMP_ORDER = ["LSC", "IRb", "SSC", "IRa"]

# Family order — matches Figure 2d for cross-panel consistency
FAMILY_ORDER = [
    "Poaceae", "Bromeliaceae", "Cyperaceae", "Juncaceae", "Typhaceae",
    "Restionaceae", "Eriocaulaceae", "Xyridaceae", "Joinvilleaceae",
    "Flagellariaceae", "Rapateaceae",
]

# ── Load data ──────────────────────────────────────────────────────────────────
summ = pd.read_csv(SUMM, sep="\t").set_index("Compartment")
fam  = pd.read_csv(FREQ, sep="\t").set_index("Family")
fam  = fam.loc[[f for f in FAMILY_ORDER if f in fam.index]]

# Overall hit percentages
overall_pct = {c: summ.loc[c, "Pct_of_all_hits"] for c in COMP_ORDER}

# ── Auto-crop Bandage image (remove white margins) ─────────────────────────────
bim = Image.open(BANDAGE).convert("RGBA")
# Create a white-background image for bbox detection
bg = Image.new("RGBA", bim.size, (255, 255, 255, 255))
bg.paste(bim, mask=bim.split()[3])
gray = bg.convert("L")
# Find bounding box of non-white pixels (threshold 250)
arr = np.array(gray)
non_white = np.where(arr < 250)
if non_white[0].size > 0:
    r0, r1 = non_white[0].min(), non_white[0].max()
    c0, c1 = non_white[1].min(), non_white[1].max()
    pad = 8  # px padding — tight crop so Bandage fills the axes
    r0 = max(0, r0 - pad); r1 = min(arr.shape[0], r1 + pad)
    c0 = max(0, c0 - pad); c1 = min(arr.shape[1], c1 + pad)
    bim_crop = bg.crop((c0, r0, c1, r1))
else:
    bim_crop = bg
bim_arr = np.array(bim_crop.convert("RGB"))

# ── Figure layout ──────────────────────────────────────────────────────────────
# 14 × 6.8 inches; left 62% = Bandage (zoomed, no donut), right 38% = family bars
fig = plt.figure(figsize=(14, 6.8), facecolor="white")
gs  = gridspec.GridSpec(
    1, 2,
    figure=fig,
    width_ratios=[1.65, 1.0],
    wspace=0.10,
)

ax_band = fig.add_subplot(gs[0, 0])   # Bandage — full height left
ax_bar  = fig.add_subplot(gs[0, 1])   # H-bars  — full height right

# ══════════════════════════════════════════════════════════════════════════════
# PANEL A — Bandage visualization
# ══════════════════════════════════════════════════════════════════════════════
ax_band.imshow(bim_arr, interpolation="lanczos")
ax_band.axis("off")

# BLAST identity rainbow color bar (horizontal, below image)
cmap_blast = LinearSegmentedColormap.from_list(
    "blast_rainbow",
    ["#3333ff", "#00ccff", "#00ff66", "#ffff00", "#ff6600", "#ff0000"],
    N=256,
)
# Colorbar — larger, more prominent, inset at bottom of Bandage panel
cax = ax_band.inset_axes([0.08, 0.03, 0.84, 0.045])
cb  = ColorbarBase(
    cax, cmap=cmap_blast,
    orientation="horizontal",
    norm=plt.Normalize(vmin=55, vmax=100),
)
cb.set_label("BLAST nucleotide identity (%)", fontsize=10, labelpad=4)
cb.ax.tick_params(labelsize=9)
cb.set_ticks([55, 70, 85, 100])

# Panel label
ax_band.text(-0.01, 1.01, "i", transform=ax_band.transAxes,
             fontsize=13, fontweight="bold", va="bottom", ha="right")

# Annotation: grey = no consensus hit
grey_patch = mpatches.Patch(color="#c8c8c8", label="No BLAST hit")
ax_band.legend(handles=[grey_patch], fontsize=9, loc="upper right",
               framealpha=0.88, borderpad=0.5)

# ══════════════════════════════════════════════════════════════════════════════
# PANEL ii — Horizontal stacked bars per family (compartment breakdown)
# ══════════════════════════════════════════════════════════════════════════════
# Reverse family order so Poaceae is at top (matching Figure 2d top-to-bottom)
fam_plot = fam.loc[::-1]
y_pos    = np.arange(len(fam_plot))

lefts = np.zeros(len(fam_plot))
for comp in COMP_ORDER:
    vals = fam_plot[comp].values
    bars = ax_bar.barh(
        y_pos, vals,
        left=lefts,
        color=COMP_COLORS[comp],
        height=0.65,
        label=comp,
    )
    lefts += vals

# Family labels
ax_bar.set_yticks(y_pos)
yticklabels = []
for fam_name in fam_plot.index:
    if fam_name == "Restionaceae":
        yticklabels.append(f"{fam_name} *")
    else:
        yticklabels.append(fam_name)
ax_bar.set_yticklabels(yticklabels, fontsize=9.5)

# Mark Restionaceae bar with asterisk annotation
rest_idx = list(fam_plot.index).index("Restionaceae")
ax_bar.annotate(
    "* elevated IRa (33%)",
    xy=(100, rest_idx),
    xytext=(102, rest_idx),
    fontsize=8,
    color="#c84040",
    va="center",
    ha="left",
    annotation_clip=False,
)

ax_bar.set_xlim(0, 100)
ax_bar.set_xlabel("% of BLAST hits", fontsize=10)
ax_bar.set_title("Per-family compartment composition", fontsize=10, pad=5, color="#333333")
ax_bar.tick_params(axis="x", labelsize=9)
ax_bar.tick_params(axis="y", length=0)
ax_bar.spines["top"].set_visible(False)
ax_bar.spines["right"].set_visible(False)
ax_bar.axvline(67.8, color="#888888", linewidth=0.8, linestyle="--", alpha=0.6)
ax_bar.text(67.8, -0.85, "mean\n67.8%", fontsize=7,
            color="#666666", ha="center", va="top")

# Compartment legend (replaces donut)
legend_patches = [
    mpatches.Patch(color=COMP_COLORS[c], label=f"{c}  {overall_pct[c]:.1f}%")
    for c in COMP_ORDER
]
ax_bar.legend(handles=legend_patches, fontsize=8.5, loc="lower right",
              framealpha=0.9, handlelength=1.2, labelspacing=0.35,
              title="Compartment", title_fontsize=9)

# Panel label
ax_bar.text(0.0, 1.04, "ii", transform=ax_bar.transAxes,
            fontsize=13, fontweight="bold", va="bottom", ha="left")

# ── Save ───────────────────────────────────────────────────────────────────────
for out_path in [f"{OUTDIR}/Figure2c_consensus_compartment.png", FIG2C]:
    plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved → {out_path}")
plt.close()
print("Done.")
