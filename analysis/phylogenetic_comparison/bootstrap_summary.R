#!/usr/bin/env Rscript
# =============================================================================
# A1: Bootstrap Support Quantification
# Graph-PAV bootstrapped tree (MineGraph, 1000 Felsenstein BS replicates)
# Parses embedded BS values from Newick; computes summary statistics & histogram
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
})

BASE  <- "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
TREES <- file.path(BASE, "phylogenetic_trees")
FIGS  <- file.path(BASE, "manuscript/figures")

# ---- 1. Load the bootstrapped PAV tree ----
tree_file <- file.path(TREES, "192_graph_minegraph_bs_tree_named.nwk")
tree <- read.tree(tree_file)
cat("Loaded tree:", tree_file, "\n")
cat("  Tips:", length(tree$tip.label), "\n")
cat("  Internal nodes:", tree$Nnode, "\n")

# ---- 2. Extract bootstrap values ----
# In ape, node.label holds internal node labels (bootstrap values as strings)
bs_raw <- tree$node.label

# Remove empty/NA/root label
bs_vals <- suppressWarnings(as.numeric(bs_raw))
bs_vals <- bs_vals[!is.na(bs_vals)]
cat("  Bootstrap values found:", length(bs_vals), "\n")
cat("  Range:", min(bs_vals), "-", max(bs_vals), "\n")

# ---- 3. Summary statistics ----
stats <- data.frame(
  tree           = "Graph_MineGraph_BS (1000 replicates, Felsenstein)",
  n_internal     = length(bs_vals),
  mean_BS        = round(mean(bs_vals), 1),
  median_BS      = round(median(bs_vals), 1),
  sd_BS          = round(sd(bs_vals), 1),
  pct_50         = round(100 * mean(bs_vals >= 50),  1),
  pct_70         = round(100 * mean(bs_vals >= 70),  1),
  pct_90         = round(100 * mean(bs_vals >= 90),  1),
  pct_95         = round(100 * mean(bs_vals >= 95),  1),
  pct_100        = round(100 * mean(bs_vals == 100), 1)
)

cat("\n--- Bootstrap Support Summary ---\n")
print(stats)

out_csv <- file.path(TREES, "bootstrap_support_summary.csv")
write.csv(stats, out_csv, row.names = FALSE)
cat("Saved:", out_csv, "\n")

# ---- 4. Histogram ----
df <- data.frame(BS = bs_vals)

# Threshold reference lines
thresholds <- data.frame(
  x     = c(70, 90, 95),
  label = c("70%", "90%", "95%"),
  col   = c("#E69F00", "#D55E00", "#CC79A7")
)

p <- ggplot(df, aes(x = BS)) +
  geom_histogram(binwidth = 5, boundary = 0,
                 fill = "#0072B2", colour = "white", alpha = 0.85) +
  geom_vline(data = thresholds, aes(xintercept = x, colour = label),
             linetype = "dashed", linewidth = 0.9) +
  scale_colour_manual(
    name   = "Threshold",
    values = setNames(thresholds$col, thresholds$label)
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(
    x = "Bootstrap support (%)",
    y = "Number of internal nodes"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "inside", legend.position.inside = c(0.12, 0.85))

out_fig <- file.path(FIGS, "FigSX_bootstrap_histogram.png")
ggsave(out_fig, p, width = 8, height = 5, dpi = 300)
cat("Saved:", out_fig, "\n")

cat("\n=== DONE A1: Bootstrap Summary ===\n")
