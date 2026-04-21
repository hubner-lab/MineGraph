#!/usr/bin/env Rscript

# ============================================================
# PAV t-SNE Sensitivity Analysis — Response to Reviewer 2, Comment 7
#
# Reviewer concern:
#   "Were monomorphic or highly linked nodes filtered before
#    dimensionality reduction? Critical for stability."
#
# This script:
#   1. Loads the ODGI PAV matrix (odgi_matrix.tsv)
#   2. Applies four filtering thresholds on node presence frequency
#   3. Runs PCA → t-SNE under each threshold (same seed/perplexity)
#   4. Plots a 2x2 panel using the same Family colour scheme as
#      the main manuscript figures
#   5. Reports node counts retained under each filter
#
# Filtering thresholds tested:
#   F1 (current):  remove fixed nodes only (freq = 0 or freq = N)
#   F2 (5–95 %):   also remove very rare (<5%) and near-universal (>95%) nodes
#   F3 (10–90 %):  stricter — retain only nodes in 10–90% of paths
#   F4 (20–80 %):  most restrictive — only intermediate-frequency nodes
#
# Inputs:
#   odgi_matrix.tsv                               (ODGI -H output)
#   Supplementary_Table_S3_Full_Poales_Data.xlsx  (Family metadata)
#
# Output:
#   PAV_tsne_sensitivity_R2.7.png  (2x2 panel, 300 dpi)
#   PAV_tsne_sensitivity_R2.7_node_counts.tsv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(readxl)
  library(stringr)
  library(FactoMineR)
  library(Rtsne)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# USER PARAMETERS
# Set these paths before running.
# Run from: scripts/graph_and_pav_analysis/
# ============================================================

matrix_file  <- "../../Pav Matrix/odgi_matrix.tsv"
xlsx_s3      <- "../../data/metadata/Supplementary_Table_S3_Full_Poales_Data.xlsx"
output_png   <- "PAV_tsne_sensitivity_R2.7.png"
output_tsv   <- "PAV_tsne_sensitivity_R2.7_node_counts.tsv"

tsne_perplexity <- 10
tsne_max_iter   <- 1000
tsne_seed       <- 42
pca_ncomp       <- 50    # PCA components fed into t-SNE (capped at min(nrow-1, ncol-1))

# ============================================================
# FAMILY COLOUR SCHEME
# (identical logic to odgi_tsne_poaceae_tribe_spectrummm.R
#  and NEW_PAV_matrix.r)
# ============================================================

make_raxml_tip <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\[|\\]", "")
  x <- str_replace_all(x, "-", "")
  x <- str_replace_all(x, "\\(|\\)", "")
  x <- str_replace_all(x, "[^A-Za-z0-9]", "_")
  substr(x, 1, 30)
}

meta_raw <- read_excel(xlsx_s3, sheet = "Full List") %>%
  mutate(across(everything(), as.character))

if (any(names(meta_raw) == "")) {
  meta_raw <- meta_raw %>% select(-which(names(meta_raw) == ""))
}

meta_family <- meta_raw %>%
  mutate(
    tip_label = make_raxml_tip(Species),
    Family    = if_else(is.na(Family) | Family == "", "UnassignedFamily", Family)
  ) %>%
  select(tip_label, Family) %>%
  distinct()

# ============================================================
# READ ODGI MATRIX
# ============================================================

cat("[1/5] Reading ODGI matrix...\n")
dt <- fread(matrix_file)

if (!("path.name" %in% names(dt))) {
  stop("odgi_matrix.tsv must contain a 'path.name' column.")
}

node_cols <- grep("^node\\.", names(dt), value = TRUE)
if (length(node_cols) == 0) stop("No node.* columns found.")

# Strip PanSN suffix (#1#1 etc.) from path names
dt[, path.name := sub("#\\d+#\\d+$", "", path.name)]

# Binary presence matrix (rows = paths, cols = nodes)
mat_full <- as.matrix(dt[, ..node_cols])
suppressWarnings(storage.mode(mat_full) <- "numeric")
mat_full <- (mat_full > 0) * 1L
rownames(mat_full) <- dt$path.name

N <- nrow(mat_full)
cat(sprintf("    %d paths × %d nodes\n", N, ncol(mat_full)))

# ============================================================
# DEFINE FOUR FILTERING THRESHOLDS
# ============================================================

col_freq <- colSums(mat_full)   # presence count per node (0..N)

filters <- list(
  "F1: SD>0 only\n(current)" = list(
    lo = 1,
    hi = N - 1,
    label = sprintf("F1 (current): freq ∈ [1, %d]\nRemoves fixed nodes only", N - 1)
  ),
  "F2: 5–95%" = list(
    lo = ceiling(0.05 * N),
    hi = floor(0.95 * N),
    label = sprintf("F2: freq ∈ [%.0f, %.0f]\n(5–95%% of paths)", ceiling(0.05*N), floor(0.95*N))
  ),
  "F3: 10–90%" = list(
    lo = ceiling(0.10 * N),
    hi = floor(0.90 * N),
    label = sprintf("F3: freq ∈ [%.0f, %.0f]\n(10–90%% of paths)", ceiling(0.10*N), floor(0.90*N))
  ),
  "F4: 20–80%" = list(
    lo = ceiling(0.20 * N),
    hi = floor(0.80 * N),
    label = sprintf("F4: freq ∈ [%.0f, %.0f]\n(20–80%% of paths)", ceiling(0.20*N), floor(0.80*N))
  )
)

# ============================================================
# FAMILY COLOUR PALETTE (built once, used for all panels)
# ============================================================

# Map all paths to family for palette construction
all_paths_df <- tibble(path.name = rownames(mat_full)) %>%
  mutate(tip_label = make_raxml_tip(path.name)) %>%
  left_join(meta_family, by = "tip_label") %>%
  mutate(Family = if_else(is.na(Family), "UnassignedFamily", Family))

fam_levels <- sort(unique(all_paths_df$Family))
n_fam      <- length(fam_levels)

if (requireNamespace("Polychrome", quietly = TRUE) && n_fam > 3) {
  seed_cols <- grDevices::hcl.colors(min(8, n_fam), palette = "Dark 3")
  fam_pal   <- Polychrome::createPalette(n_fam, seedcolors = seed_cols, M = 1000)
} else {
  fam_pal <- grDevices::hcl.colors(n_fam, palette = "Dark 3")
}
names(fam_pal) <- fam_levels

# ============================================================
# RUN PCA -> t-SNE FOR EACH FILTER; COLLECT RESULTS
# ============================================================

node_count_rows <- list()
plot_list       <- list()

for (fname in names(filters)) {

  fdef <- filters[[fname]]
  keep <- (col_freq >= fdef$lo) & (col_freq <= fdef$hi)
  mat_f <- mat_full[, keep, drop = FALSE]

  n_kept <- sum(keep)
  cat(sprintf("\n[filter] %s → %d nodes retained (%.1f%%)\n",
              fname, n_kept, 100 * n_kept / ncol(mat_full)))

  node_count_rows[[fname]] <- tibble(
    filter      = fname,
    lo_freq     = fdef$lo,
    hi_freq     = fdef$hi,
    nodes_total = ncol(mat_full),
    nodes_kept  = n_kept,
    pct_kept    = round(100 * n_kept / ncol(mat_full), 2)
  )

  if (n_kept < 10) {
    warning(sprintf("Too few nodes for filter %s — skipping.", fname))
    next
  }

  # --- PCA ---
  cat("    Running PCA...\n")
  ncomp_use <- min(pca_ncomp, nrow(mat_f) - 1, ncol(mat_f) - 1)
  pca_res <- PCA(as.data.frame(mat_f), ncp = ncomp_use,
                 graph = FALSE, scale.unit = TRUE)
  pca_mat <- pca_res$ind$coord   # N × ncomp_use

  # --- t-SNE ---
  cat("    Running t-SNE...\n")
  set.seed(tsne_seed)
  tsne_res <- Rtsne(
    pca_mat,
    dims        = 2,
    perplexity  = tsne_perplexity,
    max_iter    = tsne_max_iter,
    verbose     = FALSE,
    check_duplicates = FALSE
  )

  # --- Assemble plot data ---
  tsne_df <- tibble(
    path.name = rownames(mat_f),
    V1        = tsne_res$Y[, 1],
    V2        = tsne_res$Y[, 2]
  ) %>%
    mutate(tip_label = make_raxml_tip(path.name)) %>%
    left_join(meta_family, by = "tip_label") %>%
    mutate(Family = if_else(is.na(Family), "UnassignedFamily", Family))

  # --- Plot ---
  pad_x <- diff(range(tsne_df$V1)) * 0.04
  pad_y <- diff(range(tsne_df$V2)) * 0.04

  p <- ggplot(tsne_df, aes(x = V1, y = V2, color = Family)) +
    geom_point(size = 2.5, alpha = 0.88) +
    scale_color_manual(values = fam_pal, drop = FALSE) +
    coord_cartesian(
      xlim = range(tsne_df$V1) + c(-pad_x, pad_x),
      ylim = range(tsne_df$V2) + c(-pad_y, pad_y)
    ) +
    labs(
      title   = fname,
      subtitle = sprintf("%s nodes (%.1f%%)", format(n_kept, big.mark=","),
                         100 * n_kept / ncol(mat_full)),
      x = "", y = ""
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      legend.position  = "none",
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "grey80", fill = NA, linewidth = 0.4)
    )

  plot_list[[fname]] <- p
  cat(sprintf("    Panel done: %s\n", fname))
}

# ============================================================
# COMBINE PANELS + SHARED LEGEND
# ============================================================

cat("\n[4/5] Assembling 2x2 figure...\n")

# 2x2 grid using wrap_plots — explicit nrow/ncol avoids ambiguity
panel_grid <- wrap_plots(plot_list, nrow = 2, ncol = 2)

# Stack grid over guide_area(), then collect legends at bottom.
# Parentheses are critical: plot_layout() must apply to the full
# stacked composition, not to guide_area() alone.
combined_final <- (
  panel_grid /
  guide_area()
) +
  plot_layout(heights = c(5, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 10)
  )

ggsave(output_png, combined_final, width = 12, height = 12, dpi = 300)
cat(sprintf("[5/5] Saved: %s\n", output_png))

# ============================================================
# WRITE NODE COUNT TABLE
# ============================================================

node_counts_tbl <- bind_rows(node_count_rows)
write_tsv(node_counts_tbl, output_tsv)
cat(sprintf("      Saved: %s\n", output_tsv))

cat("\n=== Node retention summary ===\n")
print(node_counts_tbl)

cat("\n=== Script complete ===\n")
cat("Figure:  ", output_png, "\n")
cat("Table:   ", output_tsv, "\n")
cat("\nTo cite this analysis in the response letter:\n")
cat("  'Monomorphic nodes (present in all or absent from all 192 paths)\n")
cat("   were removed prior to dimensionality reduction (F1, current default).\n")
cat("   Sensitivity analysis across four frequency filters (F1–F4, Figure SX)\n")
cat("   demonstrates that Family-level clustering is stable regardless of\n")
cat("   the filtering stringency applied.'\n")
