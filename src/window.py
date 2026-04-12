#!/usr/bin/env python3

# Imports
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
import logging
import re

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
    """Remove PanSN suffixes like '#1#1'."""
    s = re.sub(r'#\d+#\d+', '', str(name))
    s = s.strip('#')
    s = re.sub(r'#{2,}', '#', s)
    return s

def read_data(filename: str) -> pd.DataFrame:
    """Load data and clean haplotype names."""
    if filename == "-":
        df = pd.read_csv(sys.stdin, sep="\t", header=None, skiprows=1)
    else:
        df = pd.read_csv(filename, sep="\t", header=None, skiprows=1)

    df[0] = df[0].apply(_clean_name)
    df.set_index(0, inplace=True)
    return df

def _adaptive_scaling(n_rows: int, n_cols: int):
    """Adaptive font and size depending on matrix size."""
    if n_rows <= 30:
        font_scale = 2.2
    elif n_rows <= 60:
        font_scale = 1.8
    elif n_rows <= 120:
        font_scale = 1.4
    else:
        font_scale = 1.0

    fig_width = max(8.0, min(24.0, 0.18 * max(40, n_cols)))
    fig_height = max(8.0, min(24.0, 0.20 * max(40, n_rows)))
    return font_scale, fig_width, fig_height

def plotter_window(df: pd.DataFrame, output: str) -> None:
    """Plot ONLY the heatmap (NO dendrogram)."""
    df = df.apply(pd.to_numeric, errors='coerce')
    df.fillna(0, inplace=True)

    n_rows, n_cols = df.shape
    font_scale, fig_w, fig_h = _adaptive_scaling(n_rows, n_cols)
    sns.set(font_scale=font_scale)

    # Make the heatmap
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    hm = sns.heatmap(
        df,
        ax=ax,
        cmap="viridis",
        cbar_kws={'label': 'Count'}
    )

    # Improve labels
    ax.tick_params(length=4, width=1.3)

    for label in ax.get_yticklabels():
        label.set_rotation(0)
        label.set_ha('right')

    for label in ax.get_xticklabels():
        label.set_rotation(90)
        label.set_ha('center')

    plt.tight_layout()
    plt.savefig(output + ".similarity_window.png", dpi=600, bbox_inches='tight', transparent=False)
    plt.savefig(output + ".similarity_window.pdf", dpi=600, bbox_inches='tight', transparent=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Window: Process TSV input and plot a heatmap (NO dendrogram)."
    )
    parser.add_argument('-i', '--input', type=str, required=True, help='Input TSV file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output prefix')
    args = parser.parse_args()

    logging.info("Reading data from %s", args.input)
    df = read_data(args.input)

    logging.info("Plotting heatmap to %s", args.output)
    plotter_window(df, args.output)
