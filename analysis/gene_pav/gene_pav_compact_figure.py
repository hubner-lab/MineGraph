#!/usr/bin/env python3
"""
Figure 2d — Compact Gene PAV Heatmap (variant genes only)
==========================================================
Loads pre-computed gene_pav_family_freq.tsv from R2.3 analysis.
Filters to variant genes: excludes universally present (core) and
universally absent genes. Keeps singletons (present in ≥1 family).
Clusters genes by shared loss/gain pattern (hierarchical clustering).

Output:
  manuscript/figures/Figure2/Figure2d_gene_pav_compact.png  ← Figure 2d
  rebuttal/exports/R2_3/gene_pav_compact.png                ← export copy
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

BASE   = "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
FREQ   = f"{BASE}/rebuttal/exports/R2_3/gene_pav_family_freq.tsv"
OUTDIR = f"{BASE}/rebuttal/exports/R2_3"
FIG2D  = f"{BASE}/manuscript/figures/Figure2/Figure2d_gene_pav_compact.png"

# ── Gene category map ─────────────────────────────────────────────────────────
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
    "ndhf":"NADH","ndhg":"NADH","ndhh":"NADH","ndhi":"NADH","ndhj":"NADH",
    "ndhk":"NADH",
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

CATEGORY_COLORS = {
    "PSI":"#2ecc71","PSII":"#27ae60","RuBisCO":"#16a085","Cyt_b6f":"#1abc9c",
    "ATPase":"#f39c12","NADH":"#e67e22","RPL":"#3498db","RPS":"#2980b9",
    "RpoABCD":"#9b59b6","rRNA":"#c0392b","Ycf":"#95a5a6","Other":"#7f8c8d",
    "ORF":"#bdc3c7",
}

FAMILY_ORDER = [
    "Poaceae","Bromeliaceae","Cyperaceae","Juncaceae","Typhaceae",
    "Restionaceae","Eriocaulaceae","Xyridaceae","Joinvilleaceae",
    "Flagellariaceae","Rapateaceae",
]

# ── Load frequency table (genes × families, values = % presence) ──────────────
freq_raw = pd.read_csv(FREQ, sep="\t", index_col="family")  # families × genes
hm = freq_raw.T  # genes × families

# Keep only families in FAMILY_ORDER that are present
ordered_fams = [f for f in FAMILY_ORDER if f in hm.columns]
hm = hm[ordered_fams]

# ── Filter: variant genes only ────────────────────────────────────────────────
# Exclude universally present (min == 100 across all families)
# Exclude universally absent  (max == 0  across all families)
min_pct = hm.min(axis=1)
max_pct = hm.max(axis=1)
variant_mask = ~((min_pct == 100) | (max_pct == 0))
hm_var = hm.loc[variant_mask].copy()
print(f"Total genes: {len(hm)} | Core (excl): {(min_pct==100).sum()} | "
      f"Absent all (excl): {(max_pct==0).sum()} | Variant (kept): {len(hm_var)}")

# ── Cluster genes by loss/gain pattern ───────────────────────────────────────
# Binary-ise at 50% threshold for clustering (presence = >50% of family has gene)
hm_bin = (hm_var > 50).astype(int)
# Use Hamming distance on binary matrix; ward linkage on gene axis
dist = pdist(hm_bin.values, metric="hamming")
Z = linkage(dist, method="ward")
gene_order = hm_var.index[leaves_list(Z)].tolist()
hm_var = hm_var.loc[gene_order]

n_genes = len(hm_var)
n_fams  = len(ordered_fams)

# ── Figure layout ─────────────────────────────────────────────────────────────
# Compact: target ~0.28 per gene row; narrow category strip
fig_h = max(6, n_genes * 0.28 + 1.5)
fig_w = max(5, n_fams * 0.7 + 2.5)

fig, (ax_cat, ax_hm) = plt.subplots(
    1, 2,
    figsize=(fig_w, fig_h),
    gridspec_kw={"width_ratios": [0.12, 1]},
    facecolor="white",
)

# ── Category colour strip ──────────────────────────────────────────────────────
cat_colors = [CATEGORY_COLORS.get(GENE_CATEGORIES.get(g, "Other"), "#bdc3c7")
              for g in gene_order]
for i, col in enumerate(cat_colors):
    ax_cat.barh(i, 1, color=col, height=1.0, align="center")
ax_cat.set_xlim(0, 1)
ax_cat.set_ylim(-0.5, n_genes - 0.5)
ax_cat.set_yticks(range(n_genes))
ax_cat.set_yticklabels(gene_order, fontsize=6)
ax_cat.set_xticks([])
ax_cat.invert_yaxis()
ax_cat.tick_params(left=False)
for sp in ax_cat.spines.values():
    sp.set_visible(False)
ax_cat.set_title("Category", fontsize=7, pad=4)

# ── Heatmap ───────────────────────────────────────────────────────────────────
cmap = LinearSegmentedColormap.from_list(
    "pav", ["#FFFFFF", "#d4e9f7", "#2c7fb8", "#0a2558"], N=256
)
mat = hm_var.values.astype(float)
im = ax_hm.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=100,
                  interpolation="none")

ax_hm.set_xticks(range(n_fams))
ax_hm.set_xticklabels(ordered_fams, fontsize=7, rotation=45, ha="right")
ax_hm.set_yticks(range(n_genes))
ax_hm.set_yticklabels([""] * n_genes)
ax_hm.invert_yaxis()

# Vertical lines between families
for j in range(1, n_fams):
    ax_hm.axvline(j - 0.5, color="white", linewidth=0.6, alpha=0.8)

# Thin horizontal lines between genes for readability
for i in range(1, n_genes):
    ax_hm.axhline(i - 0.5, color="white", linewidth=0.25, alpha=0.5)

cb = plt.colorbar(im, ax=ax_hm, orientation="vertical",
                  fraction=0.03, pad=0.02)
cb.set_label("% genomes in family with gene", fontsize=7)
cb.ax.tick_params(labelsize=6)

# Category legend (compact, top-right)
seen_cats = sorted({GENE_CATEGORIES.get(g, "Other") for g in gene_order},
                   key=lambda c: list(CATEGORY_COLORS.keys()).index(c)
                   if c in CATEGORY_COLORS else 99)
legend_patches = [mpatches.Patch(color=CATEGORY_COLORS.get(c, "#bdc3c7"), label=c)
                  for c in seen_cats]
ax_hm.legend(handles=legend_patches, fontsize=5.5, loc="lower right",
             title="Category", title_fontsize=6, framealpha=0.9, ncol=2)

plt.tight_layout(pad=0.8)

for out_path in [f"{OUTDIR}/gene_pav_compact.png", FIG2D]:
    plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved → {out_path}")
plt.close()

# ── Console summary ────────────────────────────────────────────────────────────
print(f"\nVariant genes in compact heatmap ({n_genes} total):")
for g in gene_order:
    row = hm_var.loc[g]
    absent = [f for f in ordered_fams if row[f] == 0]
    partial = [f for f in ordered_fams if 0 < row[f] < 100]
    print(f"  {g:12s}  absent={absent or '-'}  partial={partial or '-'}")
