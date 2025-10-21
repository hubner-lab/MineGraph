#!/usr/bin/env python3

# Imports
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
import logging
from scipy.cluster.hierarchy import dendrogram, linkage

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s    [%(levelname)s] - %(filename)s: %(message)s',
    datefmt='%d/%m/%Y %H:%M:%S',  # 24-hour format
    handlers=[logging.StreamHandler(stream=sys.stderr)]
)

# Global font and text settings for publication
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Liberation Sans', 'Arial'],
    'font.size': 18,
    'axes.titlesize': 22,
    'axes.labelsize': 20,
    'legend.fontsize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
})

def read_data(filename: str) -> pd.DataFrame:
    """
    Parse the "gretl window" output in DataFrame.

    :param filename: Input file name (or stdin of "-")
    :return: Data in pandas DataFrame
    """
    if filename == "-":  # Read from stdin
        df = pd.read_csv(sys.stdin, sep="\t", header=None, skiprows=1)
        df.set_index(0, inplace=True)
        return df
    else:
        df = pd.read_csv(filename, sep="\t", header=None, skiprows=1)
        df.set_index(0, inplace=True)
        return df

def plotter_window(df: pd.DataFrame, output: str) -> None:
    """
    Plot the heatmap of the data with only a left dendrogram.

    :param df: Input DataFrame
    :param output: Output file prefix
    :return: None
    """
    # Ensure the DataFrame has numeric values
    df = df.apply(pd.to_numeric, errors='coerce')
    df.fillna(0, inplace=True)

    # Compute the hierarchical clustering linkage for rows
    row_linkage = linkage(df.values, method='ward')

    # Create a figure with gridspec for the row dendrogram and heatmap
    fig = plt.figure(figsize=(12, 12))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.2, 1])

    # Plot row dendrogram
    ax_row_dendro = fig.add_subplot(gs[0, 0])
    dendrogram(row_linkage, orientation='left', ax=ax_row_dendro, no_labels=True)
    ax_row_dendro.invert_yaxis()
    ax_row_dendro.axis("off")

    # Plot heatmap
    ax_heatmap = fig.add_subplot(gs[0, 1])
    sns.heatmap(df, ax=ax_heatmap, cmap="coolwarm", cbar_kws={'label': 'Count'})

    # Save the plot
    plt.tight_layout()
    plt.savefig(output + ".similarity_window_with_tree.png", dpi=600, bbox_inches='tight', transparent=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Window: Process TSV input and plot a heatmap with a left dendrogram.")
    parser.add_argument('-i', '--input', type=str, help='Path to the input TSV file', required=True)
    parser.add_argument('-o', '--output', type=str, help='Path to the output file prefix', required=True)
    args = parser.parse_args()

    logging.info("Reading data from %s", args.input)
    df = read_data(args.input)

    logging.info("Plotting to %s", args.output)
    plotter_window(df, args.output)
