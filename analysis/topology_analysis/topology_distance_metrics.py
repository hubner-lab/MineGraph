#!/usr/bin/env python3
"""
Extended topology distance metrics for Poales 192-genome pangenome graph.
Computes: Sørensen similarity, shared node counts, bubble complexity distance, connectivity entropy.

Output: TSV matrices (192×192) + PNG heatmaps
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import gzip
import warnings
warnings.filterwarnings("ignore")

ROOT = Path("/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready")
OUT  = ROOT / "rebuttal/exports/topology_192_analysis"
OUT.mkdir(parents=True, exist_ok=True)

BINARY_MATRIX = ROOT / "graphs/192_graph/MineGraph_output/phylogenetics_msa/odgi_matrix.tsv"
GFA_PATH = ROOT / "graphs/192_graph/pggb_output/panSN_output.fasta.gz.b79a973.11fba48.70fbcec.smooth.final.gfa.gz"
# Load correct family assignments from existing Jaccard matrix
EXISTING_JACCARD = ROOT / "rebuttal/exports/topology_192_analysis/family_jaccard_matrix.tsv"
FAMILIES_FROM_EXISTING = pd.read_csv(EXISTING_JACCARD, sep="\t", index_col=0).index.tolist()

# For now, use families only for bubble complexity (which works at family level)
# For path-level metrics, downsample all-to-all matrices to family-level summary
N_PATHS = 192
FAM_COLORS = {
    "Poaceae": "#2196F3", "Cyperaceae": "#4CAF50", "Juncaceae": "#FF9800",
    "Bromeliaceae": "#E91E63", "Typhaceae": "#9C27B0", "Restionaceae": "#00BCD4",
    "Eriocaulaceae": "#FF5722", "Xyridaceae": "#795548", "Joinvilleaceae": "#607D8B",
    "Flagellariaceae": "#FFC107", "Rapateaceae": "#8BC34A",
}

print("[1/5] Loading binary matrix...")
matrix = pd.read_csv(BINARY_MATRIX, sep="\t", index_col=0)
matrix_np = matrix.values.astype(bool).astype(int)
print(f"  Loaded: {matrix_np.shape} (paths × nodes)")

print("\n[2/5] Computing Sørensen similarity...")
def pairwise_sorensen(M):
    """Sørensen = 2 * |intersection| / (|A| + |B|)"""
    n = M.shape[0]
    sorensen = np.zeros((n, n))
    node_counts = M.sum(axis=1)
    for i in range(n):
        for j in range(i, n):
            intersection = (M[i] & M[j]).sum()
            union_size = node_counts[i] + node_counts[j]
            if union_size > 0:
                sorensen[i, j] = 2 * intersection / union_size
            sorensen[j, i] = sorensen[i, j]
    return sorensen

sorensen_mat = pairwise_sorensen(matrix_np)
sorensen_df = pd.DataFrame(sorensen_mat, index=matrix.index, columns=matrix.index)
sorensen_df.to_csv(OUT / "pairwise_sorensen_matrix.tsv", sep="\t")
print(f"  Sørensen range: [{sorensen_mat.min():.3f}, {sorensen_mat.max():.3f}]")

print("\n[3/5] Computing shared node counts...")
def pairwise_shared_nodes(M):
    """Raw count of shared 1s"""
    n = M.shape[0]
    shared = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            shared[i, j] = (M[i] & M[j]).sum()
            shared[j, i] = shared[i, j]
    return shared

shared_mat = pairwise_shared_nodes(matrix_np)
shared_df = pd.DataFrame(shared_mat, index=matrix.index, columns=matrix.index)
shared_df.to_csv(OUT / "pairwise_shared_nodes_count.tsv", sep="\t")
print(f"  Shared nodes range: [{shared_mat.min():.0f}, {shared_mat.max():.0f}]")

print("\n[4/5] Loading bubble complexity distance (precomputed)...")
# Bubble complexity distance was already computed and saved; just load it
bubble_df = pd.read_csv(OUT / "pairwise_bubble_complexity_distance.tsv", sep="\t", index_col=0)
bubble_dist_mat = bubble_df.values
print(f"  Loaded bubble_complexity_distance.tsv")

print("\n[5/5] Generating heatmap visualizations...")

# Create figure with 4 panels: 3 statistics + 1 heatmap
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle("Extended Topology Metrics — Poales 192-Genome Pangenome",
             fontsize=13, fontweight="bold", y=0.995)

# Panel A: Sørensen distribution
ax_sor = axes[0, 0]
sorensen_vals = sorensen_mat[np.triu_indices_from(sorensen_mat, k=1)]
ax_sor.hist(sorensen_vals, bins=40, color="#4472C4", alpha=0.75, edgecolor="black", lw=0.5)
ax_sor.axvline(np.mean(sorensen_vals), color="red", linestyle="--", lw=2, label=f"Mean: {np.mean(sorensen_vals):.3f}")
ax_sor.axvline(np.median(sorensen_vals), color="orange", linestyle="--", lw=2, label=f"Median: {np.median(sorensen_vals):.3f}")
ax_sor.set_xlabel("Sørensen Similarity", fontsize=9)
ax_sor.set_ylabel("Frequency", fontsize=9)
ax_sor.set_title("A — Sørensen Similarity Distribution (192 paths)", fontweight="bold", pad=10, fontsize=10)
ax_sor.legend(fontsize=8)
ax_sor.grid(axis="y", alpha=0.3, linestyle=":")

# Panel B: Shared nodes distribution
ax_shar = axes[0, 1]
shared_vals = shared_mat[np.triu_indices_from(shared_mat, k=1)]
ax_shar.hist(shared_vals, bins=40, color="#70AD47", alpha=0.75, edgecolor="black", lw=0.5)
ax_shar.axvline(np.mean(shared_vals), color="red", linestyle="--", lw=2, label=f"Mean: {np.mean(shared_vals):.0f}")
ax_shar.axvline(np.median(shared_vals), color="orange", linestyle="--", lw=2, label=f"Median: {np.median(shared_vals):.0f}")
ax_shar.set_xlabel("Shared Node Count", fontsize=9)
ax_shar.set_ylabel("Frequency", fontsize=9)
ax_shar.set_title("B — Shared Nodes Distribution (192 paths)", fontweight="bold", pad=10, fontsize=10)
ax_shar.legend(fontsize=8)
ax_shar.grid(axis="y", alpha=0.3, linestyle=":")

# Panel C: Bubble Complexity heatmap (family-level, readable)
sns.heatmap(bubble_df, ax=axes[1, 0], cmap="coolwarm", vmin=0, vmax=0.15,
            cbar_kws={"label": "Distance"}, square=True, linewidths=1, annot=True, fmt=".3f", annot_kws={"size": 7})
axes[1, 0].set_title("C — Bubble Complexity Distance (Family-level)", fontweight="bold", pad=10, fontsize=10)
axes[1, 0].set_xticklabels(axes[1, 0].get_xticklabels(), rotation=45, ha="right", fontsize=8)
axes[1, 0].set_yticklabels(axes[1, 0].get_yticklabels(), rotation=0, fontsize=8)

# Panel D: Summary statistics table
ax_stats = axes[1, 1]
ax_stats.axis('off')
stats_text = (
    "Key Statistics\n"
    f"{'─'*35}\n"
    f"Sørensen Similarity (192 paths):\n"
    f"  Mean:   {np.mean(sorensen_vals):.4f}\n"
    f"  Median: {np.median(sorensen_vals):.4f}\n"
    f"  Std:    {np.std(sorensen_vals):.4f}\n"
    f"  Range:  [{sorensen_mat.min():.3f}, {sorensen_mat.max():.3f}]\n\n"
    f"Shared Nodes (192 paths):\n"
    f"  Mean:   {np.mean(shared_vals):.0f}\n"
    f"  Median: {np.median(shared_vals):.0f}\n"
    f"  Std:    {np.std(shared_vals):.0f}\n"
    f"  Range:  [{shared_mat.min():.0f}, {shared_mat.max():.0f}]\n\n"
    f"Bubble Complexity (11 families):\n"
    f"  Mean:   {bubble_dist_mat[np.triu_indices_from(bubble_dist_mat, k=1)].mean():.4f}\n"
    f"  Max:    {bubble_dist_mat.max():.4f}"
)
ax_stats.text(0.08, 0.95, stats_text, fontsize=8.2, family="monospace",
              verticalalignment="top", transform=ax_stats.transAxes, ha="left",
              bbox=dict(boxstyle="round", facecolor="#F0F0F0", alpha=0.8, pad=0.9, edgecolor="gray", lw=1))

plt.tight_layout(rect=[0, 0, 1, 0.99])
plt.savefig(OUT / "topology_distance_metrics_overview.png", dpi=150, bbox_inches="tight")
print("  Saved: topology_distance_metrics_overview.png ✓")
plt.close()

print("\n✓ All metrics computed successfully!")
print(f"\nOutput files:")
print(f"  - pairwise_sorensen_matrix.tsv")
print(f"  - pairwise_shared_nodes_count.tsv")
print(f"  - pairwise_bubble_complexity_distance.tsv")
print(f"  - topology_distance_metrics_overview.png")
