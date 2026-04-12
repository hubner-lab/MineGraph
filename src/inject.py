#!/usr/bin/env python3
"""
inject.py — Gene annotation injection and gggenes visualization.

Runs INSIDE the opggb container. Called by MineGraph.py host wrapper.

Pipeline:
  1. Build OG from GFA (if input is .gfa / .gfa.gz)
  2. odgi inject  — embed BED gene annotations as graph traversal paths
  3. odgi extract — subgraph containing representative genomes + gene paths
  4. odgi flip    — normalise path orientations before untangle
  5. odgi untangle -g — map gene paths onto genome coordinate space (gggenes TSV)
  6. odgi bin     — per-path bin coverage table
  7. R gggenes    — gene-arrow visualization (PDF + PNG)
  8. Python       — bin-coverage heatmap (PDF + PNG)
  9. inject_stats.txt — graph statistics summary
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_capture(cmd: List[str], label: str = "") -> str:
    """Run command, return stdout. Raises on non-zero exit."""
    if label:
        print(f"[inject] {label}", file=sys.stderr)
    print("[CMD]", " ".join(cmd), file=sys.stderr)
    res = subprocess.run(cmd, check=True, text=True, capture_output=True)
    return res.stdout


def run_live(cmd: List[str], label: str = "") -> int:
    """Run command with live stderr/stdout. Returns exit code."""
    if label:
        print(f"[inject] {label}", file=sys.stderr)
    print("[CMD]", " ".join(cmd), file=sys.stderr)
    return subprocess.run(cmd).returncode


def die(msg: str) -> None:
    print(f"[inject][ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def read_lines(path: str) -> List[str]:
    with open(path) as f:
        return [ln.strip() for ln in f if ln.strip()]


def write_lines(path: str, lines: List[str]) -> None:
    with open(path, "w") as f:
        for ln in lines:
            f.write(ln + "\n")


def bed_gene_names(bed_path: str) -> List[str]:
    """Return ordered unique gene names from BED column 4."""
    names: List[str] = []
    seen: set = set()
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            name = parts[3]
            if name not in seen:
                names.append(name)
                seen.add(name)
    return names


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

def step_build_og(graph_path: str, out_dir: str, threads: str) -> str:
    """Convert GFA → OG if needed. Returns path to OG."""
    if graph_path.lower().endswith(".gfa") or graph_path.lower().endswith(".gfa.gz"):
        og_path = os.path.join(out_dir, "input_graph.og")
        rc = run_live(["odgi", "build", "-g", graph_path, "-o", og_path, "-t", threads],
                      "Building OG from GFA")
        if rc != 0:
            die("odgi build failed")
        return og_path
    return graph_path


def step_inject(graph_og: str, bed: str, out_og: str, threads: str) -> None:
    rc = run_live(
        ["odgi", "inject", "-i", graph_og, "-b", bed, "-o", out_og, "-t", threads],
        "Injecting gene annotations (odgi inject)"
    )
    if rc != 0:
        die("odgi inject failed")


def get_graph_paths(og: str) -> List[str]:
    out = run_capture(["odgi", "paths", "-i", og, "-L"])
    return [p for p in out.splitlines() if p]


def step_extract_subgraph(
    injected_og: str,
    wanted_file: str,
    out_og: str,
    threads: str,
    anchor: str,
) -> None:
    rc = run_live(
        [
            "odgi", "extract",
            "-i", injected_og,
            "-r", anchor,
            "-p", wanted_file,
            "-R", wanted_file,
            "-o", out_og,
            "-t", threads,
        ],
        f"Extracting subgraph (anchor: {anchor})"
    )
    if rc != 0:
        die("odgi extract failed")


def step_flip(og_in: str, og_out: str, threads: str) -> None:
    rc = run_live(
        ["odgi", "flip", "-i", og_in, "-o", og_out, "-t", threads],
        "Normalising path orientations (odgi flip)"
    )
    if rc != 0:
        die("odgi flip failed")


def step_untangle(
    og: str,
    gene_paths_file: str,
    genome_paths_file: str,
    out_tsv: str,
    threads: str,
) -> None:
    """Run odgi untangle in gggenes mode; strip log lines; write clean TSV."""
    cmd = [
        "odgi", "untangle",
        "-i", og,
        "-R", gene_paths_file,
        "-Q", genome_paths_file,
        "-j", "0.0",
        "-t", threads,
        "-g",
    ]
    print("[inject] Running odgi untangle (gggenes mode)", file=sys.stderr)
    print("[CMD]", " ".join(cmd), file=sys.stderr)
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        die("odgi untangle failed")
    # Strip odgi log lines (start with '[') but keep header and data
    clean = [ln for ln in result.stdout.splitlines() if not ln.startswith("[")]
    with open(out_tsv, "w") as f:
        f.write("\n".join(clean) + "\n")
    data_lines = len([l for l in clean if l and not l.startswith("molecule")])
    print(f"[inject] Untangle TSV: {len(clean)} lines ({data_lines} data rows) → {out_tsv}",
          file=sys.stderr)


def step_bin(og: str, out_tsv: str, bins: str, threads: str) -> None:
    cmd = ["odgi", "bin", "-i", og, "-n", bins, "-t", threads]
    print(f"[inject] Running odgi bin (n={bins})", file=sys.stderr)
    print("[CMD]", " ".join(cmd), file=sys.stderr)
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        die("odgi bin failed")
    with open(out_tsv, "w") as f:
        f.write(result.stdout)
    print(f"[inject] Bin coverage TSV → {out_tsv}", file=sys.stderr)


def step_stats(og: str, out_txt: str) -> None:
    stats = run_capture(["odgi", "stats", "-i", og, "-S"])
    paths_out = run_capture(["odgi", "paths", "-i", og, "-L"])
    with open(out_txt, "w") as f:
        f.write("=== Graph Stats ===\n")
        f.write(stats)
        f.write("\n=== Paths ===\n")
        f.write(paths_out)
    print(f"[inject] Stats → {out_txt}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Embedded R script (gggenes gene-arrow plot)
# ---------------------------------------------------------------------------

_R_GGGENES = r"""
suppressPackageStartupMessages({
  library(ggplot2)
  library(gggenes)
  library(dplyr)
  library(stringr)
})

args        <- commandArgs(trailingOnly=TRUE)
tsv_file    <- args[1]                      # cleaned untangle TSV
out_prefix  <- args[2]                      # output file prefix
gap_thresh  <- as.numeric(args[3])          # gap for merge_close (bp)
label_map_f <- if (length(args) >= 4 && args[4] != "none") args[4] else NULL
color_map_f <- if (length(args) >= 5 && args[5] != "none") args[5] else NULL

# ── Load data ─────────────────────────────────────────────────────────────
x <- read.delim(tsv_file, stringsAsFactors=FALSE)
x$start  <- as.numeric(x$start)
x$end    <- as.numeric(x$end)
x$strand <- as.numeric(x$strand)
x <- x[!is.na(x$start) & !is.na(x$end), ]
if (nrow(x) == 0) { cat("No data rows — nothing to plot.\n"); quit(status=0) }

# ── Gene name ─────────────────────────────────────────────────────────────
# odgi inject names paths as "{BED_name}:{start}-{end}".
# Strip the coordinate suffix to recover the BED name, then extract the
# short gene name: if the name contains a dot (e.g. "Poa.rrn16") use the
# part after the last dot; otherwise use the full BED name.
x$bed_name  <- gsub(":[0-9]+-[0-9]+$", "", x$gene)
x$gene_name <- ifelse(grepl("\\.", x$bed_name),
                      sub(".*\\.", "", x$bed_name),
                      x$bed_name)

# ── Molecule label ────────────────────────────────────────────────────────
# Auto-clean: strip ":0-NNN" subrange suffix and "#sample#hap" panSN suffix;
# replace underscores with spaces.
x$label <- x$molecule
x$label <- gsub(":[0-9]+-[0-9]+$",   "", x$label)
x$label <- gsub("#[^#]+#[^#:]+$",    "", x$label)
x$label <- gsub("_",                 " ", x$label)

# Apply optional label map (TSV: path_prefix<TAB>display_label)
if (!is.null(label_map_f)) {
  lm <- read.delim(label_map_f, header=FALSE, col.names=c("prefix","label"),
                   stringsAsFactors=FALSE)
  for (i in seq_len(nrow(lm))) {
    x$label[grepl(lm$prefix[i], x$molecule, fixed=TRUE)] <- lm$label[i]
  }
}

# ── Merge close fragments (preserve IR double copies) ─────────────────────
# Fragments within gap_thresh bp of each other belong to the same gene copy;
# the two IR copies are >gap_thresh bp apart so they remain separate arrows.
merge_close <- function(df, gap_thr) {
  if (nrow(df) == 0)
    return(data.frame(start=integer(0), end=integer(0), strand=integer(0)))
  df <- df[order(df$start), ]
  result     <- list()
  cur_start  <- df$start[1]
  cur_end    <- df$end[1]
  cur_strand <- df$strand[1]
  for (i in seq_len(nrow(df))) {
    s <- df$start[i]; e <- df$end[i]; st <- df$strand[i]
    if (s - cur_end <= gap_thr) {
      cur_end <- max(cur_end, e)
    } else {
      result[[length(result)+1]] <-
        data.frame(start=cur_start, end=cur_end, strand=cur_strand)
      cur_start <- s; cur_end <- e; cur_strand <- st
    }
  }
  result[[length(result)+1]] <-
    data.frame(start=cur_start, end=cur_end, strand=cur_strand)
  do.call(rbind, result)
}

x_merged <- x %>%
  group_by(label, gene_name) %>%
  group_modify(~ merge_close(.x, gap_thr=gap_thresh)) %>%
  ungroup()

# ── Colour assignment ─────────────────────────────────────────────────────
# 20-colour palette; cycles for >20 distinct genes.
palette_20 <- c(
  "#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
  "#a65628","#f781bf","#999999","#17becf","#bcbd22",
  "#ffbb78","#aec7e8","#98df8a","#c5b0d5","#c49c94",
  "#f7b6d2","#dbdb8d","#9edae5","#393b79","#637939"
)
gene_names_uniq <- sort(unique(x_merged$gene_name))
n_g <- length(gene_names_uniq)

if (!is.null(color_map_f)) {
  cm <- read.delim(color_map_f, header=FALSE, col.names=c("gene","color"),
                   stringsAsFactors=FALSE)
  colors <- setNames(cm$color, cm$gene)
  missing_g <- setdiff(gene_names_uniq, names(colors))
  if (length(missing_g) > 0) {
    extra <- setNames(
      palette_20[((seq_along(missing_g) - 1) %% length(palette_20)) + 1],
      missing_g
    )
    colors <- c(colors, extra)
  }
} else {
  colors <- setNames(
    palette_20[((seq_along(gene_names_uniq) - 1) %% length(palette_20)) + 1],
    gene_names_uniq
  )
}

# ── Plot dimensions ───────────────────────────────────────────────────────
n_genomes <- length(unique(x_merged$label))
plot_h    <- max(3.5, n_genomes * 0.85 + 1.5)

p <- ggplot(x_merged,
       aes(xmin   = start,
           xmax   = end,
           y      = label,
           fill   = gene_name,
           forward = as.logical(strand))) +
  geom_gene_arrow(
    arrowhead_height     = unit(3,   "mm"),
    arrowhead_width      = unit(2,   "mm"),
    arrow_body_height    = unit(2.5, "mm")
  ) +
  scale_fill_manual(values = colors, name = "Gene") +
  scale_x_continuous(
    labels = function(x) paste0(round(x / 1000, 1), " kb"),
    name   = "Position (kb)"
  ) +
  labs(
    title    = "Gene positions across pangenome graph paths",
    subtitle = paste0(
      "Gene paths injected from BED file; positions resolved by odgi untangle.  ",
      "Duplicate arrows = double-copy regions (e.g. inverted repeats).  ",
      "Reversed arrows = opposite strand orientation."
    ),
    y = NULL
  ) +
  theme_genes() +
  theme(
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 8.5, colour = "grey35"),
    legend.position  = "right",
    axis.text.y      = element_text(size = 10),
    axis.text.x      = element_text(size = 9)
  )

ggsave(paste0(out_prefix, "_gene_arrows.pdf"), p, width = 14, height = plot_h)
ggsave(paste0(out_prefix, "_gene_arrows.png"), p, width = 14, height = plot_h, dpi = 250)
cat(sprintf("[R] Saved %s_gene_arrows.pdf and .png\n", out_prefix))

# ── Gene-presence summary table ───────────────────────────────────────────
presence <- x_merged %>%
  group_by(label, gene_name) %>%
  summarise(copies = n(), total_bp = sum(end - start), .groups = "drop")

write.table(presence,
            paste0(out_prefix, "_gene_presence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[R] Saved %s_gene_presence.tsv\n", out_prefix))

# ── Per-gene copy-number check (flag unexpectedly single-copy genes) ──────
copy_summary <- presence %>%
  group_by(gene_name) %>%
  summarise(
    genomes_with_gene       = n(),
    mean_copies_per_genome  = round(mean(copies), 2),
    max_copies              = max(copies),
    .groups = "drop"
  )
cat("\n[R] Copy-number summary per gene:\n")
print(as.data.frame(copy_summary))
"""

# ---------------------------------------------------------------------------
# Embedded Python bin-coverage heatmap
# ---------------------------------------------------------------------------

_PY_BIN_PLOT = """
import sys, re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

tsv_file   = sys.argv[1]
out_prefix = sys.argv[2]
n_bins_arg = int(sys.argv[3]) if len(sys.argv) > 3 else 500

def clean_name(s):
    s = re.sub(r':[0-9]+-[0-9]+$', '', str(s))
    s = re.sub(r'#[^#]+#[^#:]+$', '', s)
    return s.replace('_', ' ')

df = pd.read_csv(tsv_file, sep='\\t')
# Only genome paths (skip gene annotation paths which have very sparse coverage)
df['label'] = df['path.name'].apply(clean_name)
# Rough filter: gene paths are typically sparse (<5% bin occupancy per path)
path_occ = df.groupby('path.name')['bin'].count()
total_bins = df['bin'].max()
genome_paths = path_occ[path_occ > total_bins * 0.05].index
df = df[df['path.name'].isin(genome_paths)]

pivot = df.pivot_table(
    index='label', columns='bin',
    values='mean.cov', aggfunc='mean', fill_value=0
)
if pivot.empty:
    print("[bin_plot] No genome data to plot — skipping.", file=sys.stderr)
    sys.exit(0)

fig, ax = plt.subplots(figsize=(16, max(4, len(pivot) * 0.4 + 2)))
vmax = np.percentile(pivot.values[pivot.values > 0], 95) if pivot.values.max() > 0 else 1
im = ax.imshow(
    pivot.values, aspect='auto', cmap='YlOrRd',
    norm=mcolors.PowerNorm(gamma=0.4, vmin=0, vmax=vmax)
)
ax.set_yticks(range(len(pivot)))
ax.set_yticklabels(pivot.index, fontsize=8)
n = len(pivot.columns)
step = max(1, n // 20)
ax.set_xticks(range(0, n, step))
ax.set_xticklabels(pivot.columns[::step], rotation=45, ha='right', fontsize=7)
ax.set_xlabel('Bin', fontsize=10)
ax.set_title('Per-path bin coverage (genome paths only)', fontsize=12, fontweight='bold')
plt.colorbar(im, ax=ax, label='Mean coverage')
plt.tight_layout()
plt.savefig(f'{out_prefix}_bin_coverage.png', dpi=200)
plt.savefig(f'{out_prefix}_bin_coverage.pdf')
print(f'[bin_plot] Saved {out_prefix}_bin_coverage.png and .pdf')
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="inject.py",
        description="""
MineGraph inject pipeline (runs inside the opggb container).

Injects gene annotations from a BED file into a pangenome graph as traversal
paths, then uses odgi untangle to resolve gene positions in each representative
genome, and produces publication-quality gene-arrow and bin-coverage plots.

BED format expected (tab-separated, 6 columns):
  path_name  start  end  gene_name  score  strand
  - path_name : must match an existing path in the graph
  - gene_name : becomes the annotation path name (col 4)
  - strand    : + or - (used for arrow direction)

Representative genomes file: one graph path name per line (panSN format OK).
""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Examples:

  Basic:
    inject.py --graph graph.og --bed genes.bed \\
              --genomes reps.txt --output ./out --threads 8

  Custom anchor, larger bins, relaxed merge gap:
    inject.py --graph graph.og --bed genes.bed \\
              --genomes reps.txt --output ./out \\
              --anchor "Species_A#1#1" --bins 1000 --gap 10000

  With custom display labels and gene colours:
    inject.py --graph graph.og --bed genes.bed \\
              --genomes reps.txt --output ./out \\
              --label-map labels.tsv --color-map colors.tsv

  Skip subgraph extraction (use full injected graph — may be slow):
    inject.py --graph graph.og --bed genes.bed \\
              --genomes reps.txt --output ./out --no-subgraph
""",
    )

    parser.add_argument("--graph",    "-g", required=True,
                        help="Input pangenome graph (.gfa, .gfa.gz, or .og)")
    parser.add_argument("--bed",      "-b", required=True,
                        help="Gene annotations in BED6 format.\n"
                             "Column 1 must be a path name present in the graph.")
    parser.add_argument("--genomes",  "-r", required=True,
                        help="Text file: one representative genome path name per line.\n"
                             "The first line is used as subgraph anchor by default.")
    parser.add_argument("--output",   "-o", required=True,
                        help="Output directory (created if it does not exist).")
    parser.add_argument("--threads",  "-t", default="4",
                        help="Number of threads (default: 4)")
    parser.add_argument("--anchor",   default=None,
                        help="Path to anchor the subgraph extraction on.\n"
                             "Default: first genome in --genomes.")
    parser.add_argument("--bins",     default="500",
                        help="Number of bins for the coverage plot (default: 500).")
    parser.add_argument("--gap",      default="5000",
                        help="Gap threshold (bp) for merging fragmented gene segments\n"
                             "within the same copy (default: 5000).\n"
                             "Increase for genes with many alignment gaps.")
    parser.add_argument("--label-map", dest="label_map", default=None,
                        help="Optional TSV (path_prefix<TAB>display_label).\n"
                             "Overrides the auto-generated molecule labels in the plot.")
    parser.add_argument("--color-map", dest="color_map", default=None,
                        help="Optional TSV (gene_name<TAB>hex_color).\n"
                             "Overrides the auto-assigned gene colours.")
    parser.add_argument("--no-subgraph", action="store_true", dest="no_subgraph",
                        help="Skip subgraph extraction and run untangle on the full\n"
                             "injected graph. Useful for small graphs; may crash on\n"
                             "large graphs (>100k nodes).")
    parser.add_argument("--skip-viz", action="store_true", dest="skip_viz",
                        help="Skip all visualisation steps (R plot and bin heatmap).\n"
                             "Outputs the TSV files only.")

    return parser


def main() -> None:
    parser = build_parser()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # ── Validate inputs ──────────────────────────────────────────────────
    for flag, path in [("--graph", args.graph), ("--bed", args.bed),
                       ("--genomes", args.genomes)]:
        if not os.path.isfile(path):
            die(f"{flag}: file not found: {path}")

    out = args.output
    Path(out).mkdir(parents=True, exist_ok=True)

    genome_names = read_lines(args.genomes)
    if not genome_names:
        die("--genomes file is empty")

    gene_names = bed_gene_names(args.bed)
    if not gene_names:
        die("--bed file has no entries in column 4 (name field)")

    print(f"[inject] Graph    : {args.graph}", file=sys.stderr)
    print(f"[inject] BED      : {args.bed}  ({len(gene_names)} unique gene annotations)",
          file=sys.stderr)
    print(f"[inject] Genomes  : {len(genome_names)} representative paths", file=sys.stderr)
    print(f"[inject] Output   : {out}", file=sys.stderr)
    print(f"[inject] Threads  : {args.threads}", file=sys.stderr)
    print(f"[inject] Gap      : {args.gap} bp  (fragment merge threshold)", file=sys.stderr)
    print(f"[inject] Bins     : {args.bins}", file=sys.stderr)

    # ── Step 1: Build OG if GFA given ───────────────────────────────────
    og = step_build_og(args.graph, out, args.threads)

    # ── Step 2: odgi inject ──────────────────────────────────────────────
    injected_og = os.path.join(out, "injected.og")
    step_inject(og, args.bed, injected_og, args.threads)

    # ── Validate that at least some genome paths are in the graph ────────
    all_paths = get_graph_paths(injected_og)
    path_set  = set(all_paths)
    present_genome_paths = [p for p in genome_names if p in path_set]
    missing = [p for p in genome_names if p not in path_set]

    if missing:
        print(f"[inject][WARN] {len(missing)} genome path(s) not found in graph:",
              file=sys.stderr)
        for p in missing[:5]:
            print(f"  {p}", file=sys.stderr)
        if missing[5:]:
            print(f"  ... and {len(missing)-5} more", file=sys.stderr)
    if not present_genome_paths:
        die("None of the requested genome paths exist in the graph. "
            "Check that --genomes path names match paths in the graph. "
            "Run: odgi paths -i graph.og -L  to list available paths.")

    # Identify injected gene annotation paths (named "{BED_name}:{start}-{end}")
    injected_gene_paths = [
        p for p in all_paths
        if any(p.startswith(gn) for gn in gene_names)
    ]
    print(f"[inject] Genome paths found : {len(present_genome_paths)}", file=sys.stderr)
    print(f"[inject] Gene paths found   : {len(injected_gene_paths)}", file=sys.stderr)
    if not injected_gene_paths:
        die("No injected gene annotation paths found in graph. "
            "Verify that BED column 1 (path name) matches a path in the graph.")

    # ── Step 3: Subgraph extraction ──────────────────────────────────────
    if args.no_subgraph:
        work_og = injected_og
        print("[inject] Skipping subgraph extraction (--no-subgraph)", file=sys.stderr)
    else:
        anchor = args.anchor if args.anchor else present_genome_paths[0]
        if anchor not in path_set:
            die(f"Anchor path not found in graph: {anchor}\n"
                f"Use --anchor to specify a valid path.")

        wanted = (
            [anchor]
            + [p for p in present_genome_paths if p != anchor]
            + injected_gene_paths
        )
        wanted_file = os.path.join(out, "wanted_paths.txt")
        write_lines(wanted_file, wanted)

        subgraph_og = os.path.join(out, "subgraph.og")
        step_extract_subgraph(injected_og, wanted_file, subgraph_og, args.threads, anchor)
        work_og = subgraph_og

    # ── Step 4: odgi flip ────────────────────────────────────────────────
    flipped_og = os.path.join(out, "subgraph_flipped.og")
    step_flip(work_og, flipped_og, args.threads)

    # ── Identify paths in flipped subgraph (may have :0-N suffix appended) ─
    subgraph_paths = get_graph_paths(flipped_og)

    genome_name_set = set(genome_names)

    def is_genome_path(p: str) -> bool:
        # Match by stripping :0-N suffix to get base name
        base = p.split(":")[0]
        return base in genome_name_set

    def is_gene_path(p: str) -> bool:
        # Gene paths start with one of the BED gene names
        base = p.split(":")[0]
        return any(base.startswith(gn) or base == gn for gn in gene_names)

    sub_genome_paths = [p for p in subgraph_paths if is_genome_path(p)]
    sub_gene_paths   = [p for p in subgraph_paths if is_gene_path(p)]

    if not sub_genome_paths:
        die("No genome paths found in subgraph after extraction. "
            "This may happen if the anchor path is very short and does not overlap "
            "with the representative genomes. Try --no-subgraph or a different --anchor.")
    if not sub_gene_paths:
        die("No gene annotation paths found in subgraph. "
            "Verify that the gene annotation BED coordinates overlap with the "
            "representative genome paths.")

    print(f"[inject] Subgraph genome paths : {len(sub_genome_paths)}", file=sys.stderr)
    print(f"[inject] Subgraph gene paths   : {len(sub_gene_paths)}", file=sys.stderr)

    genome_paths_file = os.path.join(out, "genome_paths_subgraph.txt")
    gene_paths_file   = os.path.join(out, "gene_paths_subgraph.txt")
    write_lines(genome_paths_file, sub_genome_paths)
    write_lines(gene_paths_file,   sub_gene_paths)

    # ── Step 5: odgi untangle ────────────────────────────────────────────
    untangle_tsv = os.path.join(out, "untangle.tsv")
    step_untangle(flipped_og, gene_paths_file, genome_paths_file,
                  untangle_tsv, args.threads)

    # ── Step 6: odgi bin ─────────────────────────────────────────────────
    bin_tsv = os.path.join(out, "bin_coverage.tsv")
    step_bin(flipped_og, bin_tsv, args.bins, args.threads)

    # ── Step 7+8: Visualisations ─────────────────────────────────────────
    out_prefix = os.path.join(out, "inject")

    if not args.skip_viz:
        # Write R script
        r_script = os.path.join(out, "_gggenes_plot.R")
        with open(r_script, "w") as f:
            f.write(_R_GGGENES)

        r_args = [
            untangle_tsv,
            out_prefix,
            args.gap,
            args.label_map if args.label_map else "none",
            args.color_map if args.color_map else "none",
        ]
        rc = subprocess.run(["Rscript", r_script] + r_args).returncode
        if rc != 0:
            print("[inject][WARN] R visualisation failed — TSV outputs are still available",
                  file=sys.stderr)

        # Write and run Python bin-coverage plot
        bin_script = os.path.join(out, "_bin_plot.py")
        with open(bin_script, "w") as f:
            f.write(_PY_BIN_PLOT)
        rc = subprocess.run(
            ["python3", bin_script, bin_tsv, out_prefix, args.bins]
        ).returncode
        if rc != 0:
            print("[inject][WARN] Bin coverage plot failed — TSV output is still available",
                  file=sys.stderr)

    # ── Step 9: Graph stats ───────────────────────────────────────────────
    stats_file = os.path.join(out, "inject_stats.txt")
    step_stats(flipped_og, stats_file)

    # ── Final summary ─────────────────────────────────────────────────────
    print("\n[inject] ══ Done ══", file=sys.stderr)
    print(f"[inject] Outputs in: {out}", file=sys.stderr)
    if not args.skip_viz:
        print(f"[inject]   Gene-arrow plot   : inject_gene_arrows.png / .pdf",
              file=sys.stderr)
        print(f"[inject]   Bin-coverage plot : inject_bin_coverage.png / .pdf",
              file=sys.stderr)
        print(f"[inject]   Gene-presence TSV : inject_gene_presence.tsv",
              file=sys.stderr)
    print(f"[inject]   Untangle TSV      : untangle.tsv", file=sys.stderr)
    print(f"[inject]   Bin-coverage TSV  : bin_coverage.tsv", file=sys.stderr)
    print(f"[inject]   Stats             : inject_stats.txt", file=sys.stderr)


if __name__ == "__main__":
    main()
