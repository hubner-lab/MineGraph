#!/usr/bin/env python3
"""
R2.3 extension — Add 'Consensus' column to gene PAV top-2 heatmap
==================================================================
For each of the 18 selected genes, determine whether the gene sequence
has a BLASTN match in the 77,746 bp pangenome consensus using multiple
reference genomes as query sources.

Presence scale:
  100 = PRESENT   (pident >= 75%, coverage >= 50%)
   50 = PARTIAL   (hit exists: pident >= 70%, coverage >= 20%)
    0 = ABSENT    (no qualifying hit)

Outputs:
  rebuttal/exports/R2_3/gene_pav_top2cat_consensus.png
  rebuttal/exports/R2_3/consensus_gene_presence.tsv
"""

import os, subprocess, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("BioPython required: pip install biopython")

BASE      = "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
FREQ      = f"{BASE}/rebuttal/exports/R2_3/gene_pav_family_freq.tsv"
CONSENSUS = f"{BASE}/consensus/Consensus_sequence.fasta"
GB_DIR    = f"{BASE}/data/wanted_gb"
OUTDIR    = f"{BASE}/rebuttal/exports/R2_3"

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
}

CAT_ORDER = ["PSI","PSII","Cyt_b6f","NADH","RPL","RPS","RpoABCD","Ycf","Other",
             "ATPase","RuBisCO","rRNA"]

CATEGORY_COLORS = {
    "PSI":     "#2ecc71", "PSII":    "#27ae60", "RuBisCO": "#16a085",
    "Cyt_b6f": "#1abc9c", "ATPase":  "#f39c12", "NADH":    "#e67e22",
    "RPL":     "#3498db", "RPS":     "#2980b9", "RpoABCD": "#9b59b6",
    "rRNA":    "#c0392b", "Ycf":     "#7f8c8d", "Other":   "#95a5a6",
}

FAMILY_ORDER = [
    "Poaceae","Bromeliaceae","Cyperaceae","Juncaceae","Typhaceae",
    "Restionaceae","Eriocaulaceae","Xyridaceae","Joinvilleaceae",
    "Flagellariaceae","Rapateaceae",
]

# ── Step 1: Reproduce top-2 selection ──────────────────────────────────────────
freq_raw    = pd.read_csv(FREQ, sep="\t", index_col="family")
hm          = freq_raw.T
ordered_fams = [f for f in FAMILY_ORDER if f in hm.columns]
hm          = hm[ordered_fams]

min_pct = hm.min(axis=1)
max_pct = hm.max(axis=1)
hm_var  = hm.loc[~((min_pct == 100) | (max_pct == 0))].copy()

selected_genes, cat_of_gene = [], {}
for cat in CAT_ORDER:
    cat_genes = [g for g in hm_var.index if GENE_CATEGORIES.get(g, "Other") == cat]
    if not cat_genes:
        continue
    picked = hm_var.loc[cat_genes].std(axis=1).sort_values(ascending=False).head(2).index.tolist()
    selected_genes.extend(picked)
    for g in picked:
        cat_of_gene[g] = cat

hm_sel  = hm_var.loc[selected_genes]
n_genes = len(selected_genes)
n_fams  = len(ordered_fams)
print(f"Selected genes ({n_genes}): {selected_genes}")

# ── Step 2: Extract gene sequences from reference GenBanks ─────────────────────
target_refs = [
    "Zea_mays_subsp._parviglumis.gb", "Carex_pendula.gb",
    "Aechmea_fasciata.gb", "Juncus_effusus.gb", "Typha_latifolia.gb",
    "Joinvillea_ascendens.gb", "Eriocaulon_buergerianum.gb",
]
gene_seqs = {}
for ref_file in target_refs:
    ref_path = os.path.join(GB_DIR, ref_file)
    if not os.path.exists(ref_path):
        continue
    try:
        rec = next(SeqIO.parse(ref_path, "genbank"))
    except Exception:
        continue
    for feat in rec.features:
        if feat.type not in ("CDS", "rRNA", "tRNA"):
            continue
        raw = feat.qualifiers.get("gene", feat.qualifiers.get("product", []))
        if not raw:
            continue
        gname = raw[0].lower().replace("-","").replace("_","")
        if gname in gene_seqs:
            continue
        try:
            gseq = str(feat.extract(rec.seq))
            if len(gseq) >= 40:
                gene_seqs[gname] = gseq
        except Exception:
            pass

print(f"Gene sequences extracted: {len(gene_seqs)}")
print(f"Missing: {[g for g in selected_genes if g not in gene_seqs]}")

# ── Step 3: BLASTN gene sequences against consensus ────────────────────────────
query_fa = os.path.join(OUTDIR, "_query_genes.fasta")
with open(query_fa, "w") as fh:
    for g in selected_genes:
        if g in gene_seqs:
            fh.write(f">{g}\n{gene_seqs[g]}\n")

db_path = os.path.join(OUTDIR, "_blast_db", "cons")
os.makedirs(os.path.dirname(db_path), exist_ok=True)
subprocess.run(["makeblastdb", "-in", CONSENSUS, "-dbtype", "nucl", "-out", db_path],
               capture_output=True)

blast_out = subprocess.run(
    ["blastn", "-query", query_fa, "-db", db_path,
     "-outfmt", "6 qseqid sseqid pident length qlen",
     "-perc_identity", "60", "-word_size", "11", "-evalue", "1e-5"],
    capture_output=True, text=True
)

best_hits = {}
for line in blast_out.stdout.strip().split("\n"):
    if not line:
        continue
    parts = line.split("\t")
    gene, pident, alnlen, qlen = parts[0], float(parts[2]), int(parts[3]), int(parts[4])
    cov = 100.0 * alnlen / qlen
    if gene not in best_hits or best_hits[gene][0] < pident:
        best_hits[gene] = (pident, cov)

# ── Step 4: Assign presence values ────────────────────────────────────────────
consensus_gene_presence = {}
for gene in selected_genes:
    if gene not in gene_seqs:
        consensus_gene_presence[gene] = float("nan")
    elif gene not in best_hits:
        consensus_gene_presence[gene] = 0.0
    else:
        pident, cov = best_hits[gene]
        if pident >= 75 and cov >= 50:
            consensus_gene_presence[gene] = 100.0
        elif pident >= 70 and cov >= 20:
            consensus_gene_presence[gene] = 50.0
        else:
            consensus_gene_presence[gene] = 0.0

print("\nConsensus gene presence (BLASTN):")
for g, v in consensus_gene_presence.items():
    label = {100:"PRESENT", 50:"PARTIAL", 0:"ABSENT"}.get(v, "no_seq")
    bh = best_hits.get(g)
    detail = f"id={bh[0]:.1f}%, cov={bh[1]:.1f}%" if bh else "no hit"
    print(f"  {g:12s}: {label:8s}  ({detail})")

# Save TSV
presence_df = pd.DataFrame.from_dict(consensus_gene_presence, orient="index", columns=["consensus_pct"])
presence_df.index.name = "gene"
presence_df.to_csv(f"{OUTDIR}/consensus_gene_presence.tsv", sep="\t")
print(f"\nSaved TSV → {OUTDIR}/consensus_gene_presence.tsv")

# ── Step 5: Build heatmap with consensus column ────────────────────────────────
mat_fam  = hm_sel.values.astype(float)
gap_col  = np.full((n_genes, 1), np.nan)
cons_col = np.array([consensus_gene_presence[g] for g in selected_genes]).reshape(-1, 1)
mat_full = np.hstack([mat_fam, gap_col, cons_col])
n_cols_full = mat_full.shape[1]

xlabels         = ordered_fams + ["", "Consensus\n(77 kb)"]
xtick_positions = list(range(n_fams)) + [n_fams, n_fams + 1]

FS_GENE, FS_FAM, FS_CB = 15, 15, 14
ROW_H, COL_W, MARG_H, MARG_W = 0.42, 0.90, 3.0, 3.2

fig_h = n_genes * ROW_H + MARG_H
fig_w = n_cols_full * COL_W + MARG_W

fig, (ax_cat, ax_hm) = plt.subplots(
    1, 2, figsize=(fig_w, fig_h),
    gridspec_kw={"width_ratios": [0.09, 1.0], "wspace": 0.03},
    facecolor="white",
)

# Category colour strip
cat_colors = [CATEGORY_COLORS.get(cat_of_gene.get(g, "Other"), "#bdc3c7") for g in selected_genes]
for i, col in enumerate(cat_colors):
    ax_cat.barh(i, 1, color=col, height=1.0, align="center")

# Block separators
blocks = []
block_start, prev_cat = 0, None
for i, g in enumerate(selected_genes):
    cat = cat_of_gene.get(g, "Other")
    if cat != prev_cat:
        if prev_cat is not None:
            blocks.append((prev_cat, block_start, i - 1))
        block_start, prev_cat = i, cat
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

# Heatmap
cmap = LinearSegmentedColormap.from_list(
    "pav", ["#FFFFFF", "#d4e9f7", "#2c7fb8", "#0a2558"], N=256
)
masked = np.ma.masked_invalid(mat_full)
im = ax_hm.imshow(masked, aspect="auto", cmap=cmap, vmin=0, vmax=100, interpolation="none")

# Grey for gap column
gap_mask = np.where(np.isnan(mat_full), 0.5, np.nan).reshape(n_genes, n_cols_full)
ax_hm.imshow(gap_mask, aspect="auto", cmap=plt.cm.Greys, vmin=0, vmax=1,
             interpolation="none", alpha=0.4)

# Bold separator line before consensus column
ax_hm.axvline(n_fams + 0.5, color="#222222", linewidth=2.2)

# X labels
ax_hm.set_xticks(xtick_positions)
ax_hm.set_xticklabels(xlabels, fontsize=FS_FAM, rotation=40, ha="right")
for lbl in ax_hm.get_xticklabels():
    if "Consensus" in lbl.get_text():
        lbl.set_fontweight("bold")
        lbl.set_color("#1a1a2e")
        lbl.set_fontsize(FS_FAM - 1)

ax_hm.set_yticks(range(n_genes))
ax_hm.set_yticklabels([""] * n_genes)
ax_hm.invert_yaxis()

# Grid lines
for j in range(1, n_fams):
    ax_hm.axvline(j - 0.5, color="white", linewidth=0.8)
for i in range(1, n_genes):
    ax_hm.axhline(i - 0.5, color="white", linewidth=0.35, alpha=0.6)

block_boundaries = [b[2] + 0.5 for b in blocks[:-1]]
for y in block_boundaries:
    ax_hm.axhline(y, color="white", linewidth=2.5)
    ax_cat.axhline(y, color="white", linewidth=2.5)

# Colorbar
cax = ax_hm.inset_axes([1.02, 0.0, 0.045, 1.0])
cb  = plt.colorbar(im, cax=cax)
cb.set_label("% of family genomes\nwith gene  |  Consensus: BLAST coverage",
             fontsize=FS_CB - 1, labelpad=6)
cb.ax.tick_params(labelsize=FS_CB - 2)
cb.set_ticks([0, 50, 100])
cb.set_ticklabels(["0\n(absent)", "50\n(partial)", "100\n(present)"])

# Legend patches
selected_cats = list(dict.fromkeys(cat_of_gene.values()))
legend_patches = [
    mpatches.Patch(color=CATEGORY_COLORS.get(c, "#bdc3c7"), label=c)
    for c in selected_cats
]
ax_hm.legend(
    handles=legend_patches, fontsize=12, ncol=2,
    loc="lower right", title="Functional category", title_fontsize=13,
    framealpha=0.93, handlelength=1.1, borderpad=0.6,
    labelspacing=0.35, columnspacing=0.9,
)

ax_hm.set_title(
    "Gene PAV across Poales families  +  representation in 77 kb pangenome consensus",
    fontsize=13, pad=10, fontweight="bold"
)

plt.savefig(f"{OUTDIR}/gene_pav_top2cat_consensus.png", dpi=300, bbox_inches="tight", facecolor="white")
print(f"Saved → {OUTDIR}/gene_pav_top2cat_consensus.png")
plt.close()
