#!/usr/bin/env python3
"""
Figure 2c — Gene PAV: top-2 most variable genes per functional category
========================================================================
From the 60 variant genes, selects the 2 genes with highest std(% presence
across families) per functional category. Genes are grouped by category
(not re-clustered) to keep the category visual front-and-centre.

Fonts are set ~2× larger than typical to remain legible when composited
at ~0.49× scale in Figure2_full_v4.

Output:
  manuscript/figures/Figure2/Figure2d_gene_pav_top2cat.png
  rebuttal/exports/R2_3/gene_pav_top2cat.png
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
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

BASE   = "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
FREQ   = f"{BASE}/rebuttal/exports/R2_3/gene_pav_family_freq.tsv"
OUTDIR = f"{BASE}/rebuttal/exports/R2_3"
OUT    = f"{BASE}/manuscript/figures/Figure2/Figure2d_gene_pav_top2cat.png"

# ── Gene categories ───────────────────────────────────────────────────────────
GENE_CATEGORIES = {
    "psaa":"PSI","psab":"PSI","psac":"PSI","psai":"PSI","psaj":"PSI",
    "psba":"PSII","psbb":"PSII","psbc":"PSII","psbd":"PSII","psbe":"PSII",
    "psbf":"PSII","psbh":"PSII","psbi":"PSII","psbj":"PSII","psbk":"PSII",
    "psbl":"PSII","psbm":"PSII","psbn":"PSII","psbt":"PSII","psbz":"PSII","psbg":"PSII",
    "peta":"Cyt_b6f","petb":"Cyt_b6f","petd":"Cyt_b6f","pete":"Cyt_b6f",
    "petg":"Cyt_b6f","petl":"Cyt_b6f","petn":"Cyt_b6f","ccsa":"Cyt_b6f",
    "rbcl":"RuBisCO",
    "atpa":"ATPase","atpb":"ATPase","atpe":"ATPase","atpf":"ATPase",
    "atph":"ATPase","atpi":"ATPase",
    "ndha":"NADH","ndhb":"NADH","ndhc":"NADH","ndhd":"NADH","ndhe":"NADH",
    "ndhf":"NADH","ndhg":"NADH","ndhh":"NADH","ndhi":"NADH","ndhj":"NADH","ndhk":"NADH",
    "rpl2":"RPL","rpl14":"RPL","rpl16":"RPL","rpl20":"RPL","rpl22":"RPL",
    "rpl23":"RPL","rpl32":"RPL","rpl33":"RPL","rpl36":"RPL",
    "rps2":"RPS","rps3":"RPS","rps4":"RPS","rps7":"RPS","rps8":"RPS",
    "rps11":"RPS","rps12":"RPS","rps14":"RPS","rps15":"RPS","rps16":"RPS",
    "rps18":"RPS","rps19":"RPS",
    "rpoa":"RpoABCD","rpob":"RpoABCD","rpoc1":"RpoABCD","rpoc2":"RpoABCD",
    "clpp":"Other","clpp1":"Other","infa":"Other","matk":"Other",
    "accd":"Other","cema":"Other",
    "ycf1":"Ycf","ycf2":"Ycf","ycf3":"Ycf","ycf4":"Ycf",
    "rrn16":"rRNA","rrn23":"rRNA","rrn4.5":"rRNA","rrn5":"rRNA",
}

CAT_ORDER = ["PSI","PSII","Cyt_b6f","NADH","RPL","RPS","RpoABCD","Ycf","Other",
             "ATPase","RuBisCO","rRNA"]

CATEGORY_COLORS = {
    "PSI":     "#2ecc71",
    "PSII":    "#27ae60",
    "RuBisCO": "#16a085",
    "Cyt_b6f": "#1abc9c",
    "ATPase":  "#f39c12",
    "NADH":    "#e67e22",
    "RPL":     "#3498db",
    "RPS":     "#2980b9",
    "RpoABCD": "#9b59b6",
    "rRNA":    "#c0392b",
    "Ycf":     "#7f8c8d",
    "Other":   "#95a5a6",
}

FAMILY_ORDER = [
    "Poaceae","Bromeliaceae","Cyperaceae","Juncaceae","Typhaceae",
    "Restionaceae","Eriocaulaceae","Xyridaceae","Joinvilleaceae",
    "Flagellariaceae","Rapateaceae",
]

# ── Load & filter variant genes ───────────────────────────────────────────────
freq_raw = pd.read_csv(FREQ, sep="\t", index_col="family")
hm = freq_raw.T
ordered_fams = [f for f in FAMILY_ORDER if f in hm.columns]
hm = hm[ordered_fams]

min_pct = hm.min(axis=1)
max_pct = hm.max(axis=1)
hm_var  = hm.loc[~((min_pct == 100) | (max_pct == 0))].copy()
print(f"Variant genes: {len(hm_var)}")

# ── Select top-2 most variable per category (highest std across families) ─────
selected_genes = []
cat_of_gene    = {}   # gene → category for colour strip
n_top = 2

for cat in CAT_ORDER:
    cat_genes = [g for g in hm_var.index if GENE_CATEGORIES.get(g, "Other") == cat]
    if not cat_genes:
        continue
    var_scores = hm_var.loc[cat_genes].std(axis=1).sort_values(ascending=False)
    picked = var_scores.head(n_top).index.tolist()
    selected_genes.extend(picked)
    for g in picked:
        cat_of_gene[g] = cat

hm_sel = hm_var.loc[selected_genes]
n_genes = len(hm_sel)
n_fams  = len(ordered_fams)
print(f"Selected (top-{n_top}/cat): {n_genes} genes across {len(set(cat_of_gene.values()))} categories")

# ── Figure layout ─────────────────────────────────────────────────────────────
# Fonts set ~2× print size — figure embedded at ~0.49× in composite,
# so source 14pt → ~7pt at final output.
FS_GENE   = 15   # gene name y-axis labels
FS_FAM    = 16   # family x-axis labels
FS_CB     = 14   # colorbar label
FS_LEG    = 15   # category legend patches
FS_LEG_TI = 14   # legend section title

ROW_H  = 0.42   # inches per gene row
COL_W  = 0.90   # inches per family column
MARG_H = 2.8    # vertical margin (title + x-labels + top pad)
MARG_W = 3.2    # left margin (cat strip + gene labels) + right (colorbar)

fig_h = n_genes * ROW_H + MARG_H
fig_w = n_fams  * COL_W + MARG_W

fig, (ax_cat, ax_hm) = plt.subplots(
    1, 2,
    figsize=(fig_w, fig_h),
    gridspec_kw={"width_ratios": [0.09, 1.0], "wspace": 0.03},
    facecolor="white",
)

# ── Category colour strip ─────────────────────────────────────────────────────
cat_colors = [CATEGORY_COLORS.get(cat_of_gene.get(g, "Other"), "#bdc3c7")
              for g in selected_genes]
for i, col in enumerate(cat_colors):
    ax_cat.barh(i, 1, color=col, height=1.0, align="center")

# Category group labels centred on each block
seen, block_start, prev_cat = {}, 0, None
blocks = []
for i, g in enumerate(selected_genes):
    cat = cat_of_gene.get(g, "Other")
    if cat != prev_cat:
        if prev_cat is not None:
            blocks.append((prev_cat, block_start, i - 1))
        block_start = i
        prev_cat = cat
blocks.append((prev_cat, block_start, len(selected_genes) - 1))

ax_cat.set_xlim(0, 1)
ax_cat.set_ylim(-0.5, n_genes - 0.5)
ax_cat.set_yticks(range(n_genes))
ax_cat.set_yticklabels(selected_genes, fontsize=FS_GENE, fontfamily="monospace")
ax_cat.set_xticks([])
ax_cat.invert_yaxis()
ax_cat.tick_params(left=False, pad=2)
for sp in ax_cat.spines.values():
    sp.set_visible(False)

# ── Heatmap ───────────────────────────────────────────────────────────────────
cmap = LinearSegmentedColormap.from_list(
    "pav", ["#FFFFFF", "#d4e9f7", "#2c7fb8", "#0a2558"], N=256
)
mat = hm_sel.values.astype(float)
im  = ax_hm.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=100,
                   interpolation="none")

# Family x-axis labels
ax_hm.set_xticks(range(n_fams))
ax_hm.set_xticklabels(ordered_fams, fontsize=FS_FAM, rotation=40, ha="right")
ax_hm.set_yticks(range(n_genes))
ax_hm.set_yticklabels([""] * n_genes)
ax_hm.invert_yaxis()

# Grid lines
for j in range(1, n_fams):
    ax_hm.axvline(j - 0.5, color="white", linewidth=0.8)
for i in range(1, n_genes):
    ax_hm.axhline(i - 0.5, color="white", linewidth=0.35, alpha=0.6)

# Category block separators (thicker white lines between groups)
block_boundaries = [b[2] + 0.5 for b in blocks[:-1]]
for y in block_boundaries:
    ax_hm.axhline(y, color="white", linewidth=2.5)
    ax_cat.axhline(y, color="white", linewidth=2.5)

# Colorbar — inset at top-right of heatmap axes
cax = ax_hm.inset_axes([1.02, 0.0, 0.045, 1.0])
cb  = plt.colorbar(im, cax=cax)
cb.set_label("% of family genomes\nwith gene", fontsize=FS_CB, labelpad=6)
cb.ax.tick_params(labelsize=FS_CB - 2)
cb.set_ticks([0, 25, 50, 75, 100])

# ── Compact inline category legend (inside heatmap axes, bottom-left) ────────
selected_cats = list(dict.fromkeys(cat_of_gene.values()))  # preserve insertion order
legend_patches = [
    mpatches.Patch(color=CATEGORY_COLORS.get(c, "#bdc3c7"), label=c)
    for c in selected_cats
]
ax_hm.legend(
    handles=legend_patches,
    fontsize=12,
    ncol=2,
    loc="lower right",
    title="Functional category",
    title_fontsize=13,
    framealpha=0.93,
    handlelength=1.1,
    borderpad=0.6,
    labelspacing=0.35,
    columnspacing=0.9,
)

# ── Save ─────────────────────────────────────────────────────────────────────
plt.tight_layout(pad=0.6)
for path in [f"{OUTDIR}/gene_pav_top2cat.png", OUT]:
    plt.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved → {path}")
plt.close()

print(f"\nSelected genes ({n_genes}):")
for b in blocks:
    cat, s, e = b
    print(f"  {cat:12s}: {selected_genes[s:e+1]}")
