#!/usr/bin/env python3

# Imports
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
import logging
import re
from scipy.cluster.hierarchy import dendrogram, linkage

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s    [%(levelname)s] - %(filename)s: %(message)s',
    datefmt='%d/%m/%Y %H:%M:%S',
    handlers=[logging.StreamHandler(stream=sys.stderr)]
)

# Global publication-style defaults
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Liberation Sans', 'Arial'],
    'axes.titlesize': 22,
    'axes.labelsize': 20,
    'legend.fontsize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'axes.linewidth': 1.5,
    'font.weight': 'bold'
})

def _clean_name(name: str) -> str:
    """
    Remove PanSN suffixes like '#1#1' (or any '#<num>#<num>' anywhere),
    then strip stray '#' at the ends.
    """
    s = re.sub(r'#\d+#\d+', '', str(name))
    s = s.strip('#')
    # Collapse any accidental double separators (e.g., 'AAA##BBB' -> 'AAA#BBB')
    s = re.sub(r'#{2,}', '#', s)
    return s

def read_data(filename: str) -> pd.DataFrame:
    """
    Parse the 'gretl window' output in a DataFrame, set first column as index,
    and clean haplotype names.
    """
    if filename == "-":
        df = pd.read_csv(sys.stdin, sep="\t", header=None, skiprows=1)
    else:
        df = pd.read_csv(filename, sep="\t", header=None, skiprows=1)

    # Column 0 holds haplotype/path names in your current script
    df[0] = df[0].apply(_clean_name)
    df.set_index(0, inplace=True)

    return df

def _adaptive_scaling(n_rows: int, n_cols: int):
    """
    Return (font_scale, fig_width, fig_height) based on matrix size.
    Larger labels for small matrices; reasonable scaling for large ones.
    """
    # Font scale tuned for publication readability
    if n_rows <= 30:
        font_scale = 2.2   # ~22 pt ticks
    elif n_rows <= 60:
        font_scale = 1.8   # ~18 pt ticks
    elif n_rows <= 120:
        font_scale = 1.4   # ~14 pt ticks
    else:
        font_scale = 1.0   # ~10 pt ticks

    # Figure size adapts to both dimensions but capped to avoid monster PDFs
    fig_width = max(8.0, min(24.0, 0.18 * max(40, n_cols)))
    fig_height = max(8.0, min(24.0, 0.20 * max(40, n_rows)))
    return font_scale, fig_width, fig_height

def plotter_window(df: pd.DataFrame, output: str) -> None:
    """
    Plot heatmap with a left dendrogram (row clustering only),
    cleaned labels, and adaptive font/size.
    """
    # Ensure numeric values
    df = df.apply(pd.to_numeric, errors='coerce')
    df.fillna(0, inplace=True)

    # Row clustering (haplotypes)
    row_linkage = linkage(df.values, method='ward')

    n_rows, n_cols = df.shape
    font_scale, fig_w, fig_h = _adaptive_scaling(n_rows, n_cols)
    sns.set(font_scale=font_scale)

    # Figure layout: left dendrogram, right heatmap
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.20, 1.00])

    # Row dendrogram (no labels on the dendrogram axis)
    ax_row_dendro = fig.add_subplot(gs[0, 0])
    dendrogram(row_linkage, orientation='left', ax=ax_row_dendro, no_labels=True)
    ax_row_dendro.invert_yaxis()
    ax_row_dendro.axis("off")

    # Heatmap
    ax_heatmap = fig.add_subplot(gs[0, 1])
    hm = sns.heatmap(
        df,
        ax=ax_heatmap,
        cmap="viridis",
        cbar_kws={'label': 'Count'}
    )

    # Improve tick label readability
    ax_heatmap.tick_params(length=4, width=1.2)
    for label in ax_heatmap.get_yticklabels():
        label.set_rotation(0)      # keep haplotype labels horizontal
        label.set_ha('right')
    for label in ax_heatmap.get_xticklabels():
        label.set_rotation(90)     # columns (windows) vertical to save space
        label.set_ha('center')

    plt.tight_layout()
    plt.savefig(output + ".similarity_window_with_tree.png", dpi=600, bbox_inches='tight', transparent=False)
    plt.savefig(output + ".similarity_window_with_tree.pdf", dpi=600, bbox_inches='tight', transparent=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Window: Process TSV input and plot a heatmap with a left dendrogram (publication-ready)."
    )
    parser.add_argument('-i', '--input', type=str, help='Path to the input TSV file', required=True)
    parser.add_argument('-o', '--output', type=str, help='Path to the output file prefix', required=True)
    args = parser.parse_args()

    logging.info("Reading data from %s", args.input)
    df = read_data(args.input)

    logging.info("Plotting to %s", args.output)
    plotter_window(df, args.output)

