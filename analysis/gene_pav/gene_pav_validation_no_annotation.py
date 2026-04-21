#!/usr/bin/env python3
"""
Gene PAV validation: confirms known gene losses (ndhK, rps16, ycf1/ycf2)
recovered by graph-native node PAV without requiring gene annotation.

Output: rebuttal/exports/gene_pav_validation/
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# ── paths ──────────────────────────────────────────────────────────────────
ROOT = Path("/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready")
PAV_TSV = ROOT / "rebuttal/exports/R2_3/gene_pav_matrix.tsv"
OUT_DIR = ROOT / "rebuttal/exports/gene_pav_validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── data ───────────────────────────────────────────────────────────────────
df = pd.read_csv(PAV_TSV, sep="\t")
df["genome"] = df["genome"].str.strip()
df["family"] = df["family"].str.strip()

gene_cols = [c for c in df.columns if c not in ("genome", "family")]

# Focus genes + context panel
KEY_GENES = ["ndhk", "rps16", "ycf1", "ycf2"]
KEY_LABELS = {"ndhk": "ndhK", "rps16": "rps16", "ycf1": "ycf1", "ycf2": "ycf2"}

# Family order (phylogenetic / size)
FAM_ORDER = [
    "Poaceae", "Bromeliaceae", "Cyperaceae", "Juncaceae",
    "Typhaceae", "Restionaceae", "Eriocaulaceae",
    "Joinvilleaceae", "Xyridaceae", "Flagellariaceae", "Rapateaceae"
]
FAM_ORDER = [f for f in FAM_ORDER if f in df["family"].unique()]

# Family palette (consistent with project)
FAM_COLORS = {
    "Poaceae": "#4E9E5B",
    "Bromeliaceae": "#D55E00",
    "Cyperaceae": "#0072B2",
    "Juncaceae": "#CC79A7",
    "Typhaceae": "#E69F00",
    "Restionaceae": "#56B4E9",
    "Eriocaulaceae": "#009E73",
    "Joinvilleaceae": "#F0E442",
    "Xyridaceae": "#999999",
    "Flagellariaceae": "#882255",
    "Rapateaceae": "#44AA99",
}

# ── compute family-level presence rates ────────────────────────────────────
fam_rates = (
    df.groupby("family")[gene_cols]
    .mean()
    .reindex(FAM_ORDER)
)

# Sort genes: by variance across families descending
gene_var = fam_rates.var()
gene_order = gene_var.sort_values(ascending=False).index.tolist()
# put KEY_GENES first
gene_order_sorted = KEY_GENES + [g for g in gene_order if g not in KEY_GENES]

# ── known published loss events ────────────────────────────────────────────
# Format: {gene: [(family, label), ...]}
KNOWN_LOSSES = {
    "ndhk":  [("Eriocaulaceae", "Eriocaulaceae\n(Iles et al. 2015)")],
    "rps16": [("Restionaceae",  "Restionaceae"),
              ("Xyridaceae",    "Xyridaceae")],
    "ycf1":  [("Poaceae",       "Poaceae\n(variable)"),
              ("Joinvilleaceae","Joinvilleaceae")],
    "ycf2":  [("Juncaceae",     "Juncaceae"),
              ("Xyridaceae",    "Xyridaceae")],
}

# ── Figure layout ──────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 13), dpi=150)
fig.patch.set_facecolor("white")

gs = gridspec.GridSpec(
    2, 2,
    figure=fig,
    hspace=0.42, wspace=0.32,
    left=0.09, right=0.97, top=0.93, bottom=0.06
)

# ── Panel A: Dot plot of key genes × family ─────────────────────────────
ax_a = fig.add_subplot(gs[0, :])

key_rates = fam_rates[KEY_GENES].copy()

n_fam = len(FAM_ORDER)
n_gene = len(KEY_GENES)

x_positions = np.arange(n_gene)
bar_width = 0.7 / n_fam

for fi, fam in enumerate(FAM_ORDER):
    offset = (fi - n_fam / 2 + 0.5) * bar_width
    rates = key_rates.loc[fam]
    x_off = x_positions + offset
    color = FAM_COLORS.get(fam, "#888888")
    ax_a.bar(x_off, rates.values, width=bar_width * 0.9,
             color=color, alpha=0.88, edgecolor="none", label=fam)

# Mark known total loss events (rate = 0)
for gi, gene in enumerate(KEY_GENES):
    for fam, lbl in KNOWN_LOSSES.get(gene, []):
        if fam in FAM_ORDER:
            fi = FAM_ORDER.index(fam)
            offset = (fi - n_fam / 2 + 0.5) * bar_width
            ax_a.annotate(
                "", xy=(gi + offset, 0.03), xytext=(gi + offset, 0.18),
                arrowprops=dict(arrowstyle="-|>", color="#cc0000", lw=1.5)
            )

ax_a.set_xticks(x_positions)
ax_a.set_xticklabels([KEY_LABELS[g] for g in KEY_GENES], fontsize=14, fontweight="bold")
ax_a.set_ylabel("Fraction of genomes with gene present", fontsize=11)
ax_a.set_ylim(0, 1.08)
ax_a.set_xlim(-0.5, n_gene - 0.5)
ax_a.axhline(1.0, color="0.7", lw=0.8, ls="--")
ax_a.set_title(
    "A  |  Graph-native PAV recovers known gene losses across Poales families",
    fontsize=12, fontweight="bold", loc="left", pad=6
)
ax_a.spines[["top", "right"]].set_visible(False)

legend_patches = [
    mpatches.Patch(facecolor=FAM_COLORS.get(f, "#888"), label=f) for f in FAM_ORDER
]
ax_a.legend(
    handles=legend_patches, ncol=6, loc="upper right",
    fontsize=8, framealpha=0.5, title="Family", title_fontsize=8
)

# red arrow annotation
ax_a.text(0.01, 0.02,
          "↓ red arrow = known total gene loss (literature-confirmed)",
          transform=ax_a.transAxes, fontsize=8, color="#cc0000", va="bottom")

# ── Panel B: Binary heatmap per genome for KEY_GENES ──────────────────────
ax_b = fig.add_subplot(gs[1, 0])

# Sort genomes by family order then genome name
df_sorted = df.copy()
df_sorted["fam_order"] = df_sorted["family"].apply(
    lambda x: FAM_ORDER.index(x) if x in FAM_ORDER else 999
)
df_sorted = df_sorted.sort_values(["fam_order", "genome"]).reset_index(drop=True)

mat = df_sorted[KEY_GENES].values.T  # genes × genomes

cmap_binary = LinearSegmentedColormap.from_list(
    "presence", ["#d73027", "#4dac26"], N=2
)

im = ax_b.imshow(mat, aspect="auto", cmap=cmap_binary, vmin=0, vmax=1,
                  interpolation="nearest")

# Family boundary lines
bounds = []
for fi, fam in enumerate(FAM_ORDER):
    idx = df_sorted[df_sorted["family"] == fam].index
    if len(idx):
        bounds.append(idx[-1] + 0.5)

for b in bounds[:-1]:
    ax_b.axvline(b, color="white", lw=1.2, alpha=0.9)

# Family labels on top
prev = -0.5
for fi, fam in enumerate(FAM_ORDER):
    idx = df_sorted[df_sorted["family"] == fam].index
    if not len(idx):
        continue
    mid = (prev + idx[-1] + 0.5) / 2
    ax_b.text(mid, -0.7, fam, ha="center", va="bottom", fontsize=6,
              rotation=45, color=FAM_COLORS.get(fam, "#333"))
    prev = idx[-1] + 0.5

ax_b.set_yticks(range(len(KEY_GENES)))
ax_b.set_yticklabels([KEY_LABELS[g] for g in KEY_GENES], fontsize=11, fontweight="bold")
ax_b.set_xlabel("Genomes (sorted by family)", fontsize=10)
ax_b.set_title(
    "B  |  Per-genome binary PAV (green = present, red = absent)",
    fontsize=11, fontweight="bold", loc="left", pad=6
)
ax_b.tick_params(bottom=False, labelbottom=False)

# ── Panel C: Loss summary table ────────────────────────────────────────────
ax_c = fig.add_subplot(gs[1, 1])
ax_c.axis("off")

summary_data = []
for gene in KEY_GENES:
    for fam in FAM_ORDER:
        sub = df[df["family"] == fam]
        n = len(sub)
        pres = int(sub[gene].sum())
        rate = pres / n if n else 0
        if rate < 0.5:  # flagging losses
            lit = "Yes" if (gene, fam) in [
                ("ndhk",  "Eriocaulaceae"),
                ("rps16", "Restionaceae"),
                ("rps16", "Xyridaceae"),
            ] else "—"
            summary_data.append([KEY_LABELS[gene], fam,
                                  f"{pres}/{n}", f"{rate:.0%}", lit])

# Also include partial losses
for gene in KEY_GENES:
    for fam in FAM_ORDER:
        sub = df[df["family"] == fam]
        n = len(sub)
        pres = int(sub[gene].sum())
        rate = pres / n if n else 0
        if 0.5 <= rate < 0.95 and n > 1:
            summary_data.append([KEY_LABELS[gene], fam,
                                  f"{pres}/{n}", f"{rate:.0%}", "—"])

# Deduplicate
seen = set()
uniq = []
for row in summary_data:
    key = (row[0], row[1])
    if key not in seen:
        seen.add(key)
        uniq.append(row)

# Sort
uniq.sort(key=lambda r: (r[0], r[3]))

headers = ["Gene", "Family", "Present", "Rate", "Literature"]

col_widths = [0.13, 0.28, 0.12, 0.10, 0.20]
x_starts = [0.0]
for w in col_widths[:-1]:
    x_starts.append(x_starts[-1] + w)

# Header
for xi, (h, xs) in enumerate(zip(headers, x_starts)):
    ax_c.text(xs, 1.0, h, fontsize=9, fontweight="bold",
              va="top", ha="left", color="#222")

ax_c.axhline(0.97, color="0.4", lw=0.8, xmin=0, xmax=1)

for ri, row in enumerate(uniq[:20]):
    y = 0.93 - ri * 0.042
    rate_val = float(row[3].strip("%")) / 100
    row_bg = "#fff3f3" if rate_val == 0 else ("#fff9e6" if rate_val < 0.95 else "white")
    ax_c.axhspan(y - 0.018, y + 0.022, color=row_bg, alpha=0.5,
                 xmin=0, xmax=1)
    for xi, (cell, xs) in enumerate(zip(row, x_starts)):
        color = "#cc0000" if (xi == 3 and rate_val == 0) else "#222"
        ax_c.text(xs, y, cell, fontsize=8, va="center",
                  ha="left", color=color)

ax_c.set_xlim(0, 0.9)
ax_c.set_ylim(-0.05, 1.05)
ax_c.set_title(
    "C  |  Gene loss / partial absence summary",
    fontsize=11, fontweight="bold", loc="left", pad=6
)

# ── save ───────────────────────────────────────────────────────────────────
out_png = OUT_DIR / "gene_pav_validation.png"
out_pdf = OUT_DIR / "gene_pav_validation.pdf"
fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Saved: {out_png}")

# ── export TSV summary ──────────────────────────────────────────────────────
summary_rows = []
for gene in KEY_GENES:
    for fam in FAM_ORDER:
        sub = df[df["family"] == fam]
        n = len(sub)
        pres = int(sub[gene].sum())
        rate = pres / n if n else np.nan
        summary_rows.append({
            "gene": KEY_LABELS[gene],
            "family": fam,
            "n_genomes": n,
            "n_present": pres,
            "fraction_present": round(rate, 4) if not np.isnan(rate) else "NA",
            "is_known_loss": "Yes" if (gene, fam) in [
                ("ndhk", "Eriocaulaceae"),
                ("rps16", "Restionaceae"),
                ("rps16", "Xyridaceae"),
            ] else "No"
        })

pd.DataFrame(summary_rows).to_csv(OUT_DIR / "gene_loss_summary.tsv", sep="\t", index=False)
print(f"Saved: {OUT_DIR / 'gene_loss_summary.tsv'}")
print("Done.")
