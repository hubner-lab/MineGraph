#!/usr/bin/env Rscript
# =============================================================================
# Extended tree comparison: 709-taxon MineGraph BS tree vs published references
# + Direct comparison of 192 BS tree vs published references (expanded metrics)
#
# Purpose:
#   Broader reviewer validation — the 709 tree represents the full Poales
#   sampling. Pruning to common taxa with Givnish2018, GPWG2, BK2015 gives
#   independent validation of the graph-based approach at larger scale.
#
# Outputs:
#   - 709_vs_published_distances.csv
#   - 192_vs_published_extended.csv   (supplement existing A2 results)
#   - FigSX_709_vs_published_heatmap.png
#   - FigSX_709_bs_histogram.png
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(TreeDist)
})
if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages("patchwork", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(patchwork))
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap", repos = "https://cloud.r-project.org")

BASE  <- "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
TREES <- file.path(BASE, "phylogenetic_trees")
FIGS  <- file.path(BASE, "manuscript/figures")

cat("=== 709-taxon MineGraph BS tree vs Published References ===\n\n")

# ---- 1. Load trees ----
tree_709 <- read.tree(file.path(TREES, "709_graph_minegraph_bs_tree_named.nwk"))
if (is.null(tree_709)) stop("709 named tree not found — run process_709_tree.py first")
cat(sprintf("709 MineGraph BS tree: %d tips\n", length(tree_709$tip.label)))

# Published references
ref_trees <- list(
  GPWG2_2012  = read.tree(file.path(TREES, "GPWG2_2012_poaceae_renamed.nwk")),
  Givnish2018 = read.tree(file.path(TREES, "Givnish2018_monocot_renamed.nwk")),
  BK2015      = read.tree(file.path(TREES, "BK2015_MCC_POALES_renamed.nwk"))
)

# Also load 192 trees for cross-scale comparison
tree_192_bs   <- read.tree(file.path(TREES, "192_graph_minegraph_bs_tree_named.nwk"))
tree_192_ward <- read.tree(file.path(TREES, "192_graph_pav_tree_named.nwk"))
tree_192_rax  <- read.tree(file.path(TREES, "192_graph_RAxML_tree_named.nwk"))

for (nm in names(ref_trees)) {
  cat(sprintf("  %s: %d tips\n", nm, length(ref_trees[[nm]]$tip.label)))
}

# ---- 2. Helper: compute all metrics for a pair ----
compute_metrics <- function(t1, t2, label1, label2) {
  common <- intersect(t1$tip.label, t2$tip.label)
  if (length(common) < 4) {
    return(data.frame(tree1=label1, tree2=label2, n_common=length(common),
                      nRF=NA, CID=NA, cophenetic_r=NA,
                      pct_shared_bipart=NA, stringsAsFactors=FALSE))
  }
  p1 <- drop.tip(t1, setdiff(t1$tip.label, common))
  p2 <- drop.tip(t2, setdiff(t2$tip.label, common))
  p1$node.label <- NULL
  p2$node.label <- NULL

  nRF_val  <- tryCatch(RobinsonFoulds(p1, p2, normalize=TRUE), error=function(e) NA_real_)
  CID_val  <- tryCatch(ClusteringInfoDistance(p1, p2), error=function(e) NA_real_)
  coph_r   <- tryCatch(
    cor(as.vector(cophenetic(p1)), as.vector(cophenetic(p2))),
    error=function(e) NA_real_)

  # Shared bipartitions
  bp1 <- tryCatch({
    parts <- prop.part(p1)
    sapply(parts, function(x) paste(sort(p1$tip.label[x]), collapse=","))
  }, error=function(e) character(0))
  bp2 <- tryCatch({
    parts <- prop.part(p2)
    sapply(parts, function(x) paste(sort(p2$tip.label[x]), collapse=","))
  }, error=function(e) character(0))

  pct_shared <- if (length(bp1) > 0)
    round(100 * sum(bp1 %in% bp2) / length(bp1), 1) else NA_real_

  data.frame(
    tree1             = label1,
    tree2             = label2,
    n_common          = length(common),
    nRF               = round(nRF_val, 3),
    CID               = round(CID_val, 3),
    cophenetic_r      = round(coph_r, 3),
    pct_shared_bipart = pct_shared,
    stringsAsFactors  = FALSE
  )
}

# ---- 3. 709 vs published ----
cat("\n--- 709 MineGraph BS vs. Published References ---\n")
rows_709 <- list()
for (ref_name in names(ref_trees)) {
  cat(sprintf("  Computing 709 vs %s ...\n", ref_name))
  rows_709[[ref_name]] <- compute_metrics(tree_709, ref_trees[[ref_name]],
                                           "MineGraph_BS_709", ref_name)
}
df_709 <- do.call(rbind, rows_709)
print(df_709)
write.csv(df_709, file.path(TREES, "709_vs_published_distances.csv"), row.names=FALSE)
cat("Saved: 709_vs_published_distances.csv\n\n")

# ---- 4. 192 vs published (extended) ----
cat("--- 192 MineGraph BS, Ward, RAxML vs. Published References ---\n")
trees_192 <- list(
  MineGraph_BS_192 = tree_192_bs,
  PAV_Ward_192     = tree_192_ward,
  RAxML_192        = tree_192_rax
)
rows_192 <- list()
for (t1_nm in names(trees_192)) {
  for (ref_name in names(ref_trees)) {
    cat(sprintf("  %s vs %s ...\n", t1_nm, ref_name))
    rows_192[[paste(t1_nm, ref_name, sep="__")]] <-
      compute_metrics(trees_192[[t1_nm]], ref_trees[[ref_name]], t1_nm, ref_name)
  }
}
df_192 <- do.call(rbind, rows_192)
print(df_192)
write.csv(df_192, file.path(TREES, "192_vs_published_extended.csv"), row.names=FALSE)
cat("Saved: 192_vs_published_extended.csv\n\n")

# ---- 5. Figure: nRF heatmap (709 + 192 combined) ----
cat("Producing distance figures ...\n")

# Combine 709 and 192 comparisons
all_df <- rbind(df_709, df_192)

# Descriptive method labels distinguishing graph-PAV from gene-based trees
method_labels <- c(
  "MineGraph_BS_709" = "Evograph PAV-BS (709 taxa)",
  "MineGraph_BS_192" = "Evograph PAV-BS (192 taxa)",
  "PAV_Ward_192"     = "Evograph PAV-Ward (192 taxa)",
  "RAxML_192"        = "Gene-based RAxML (192 taxa)"
)
method_colors <- c(
  "Evograph PAV-BS (709 taxa)"    = "#1B7837",
  "Evograph PAV-BS (192 taxa)"    = "#5AAE61",
  "Evograph PAV-Ward (192 taxa)"  = "#E69F00",
  "Gene-based RAxML (192 taxa)"   = "#D73027"
)
all_df$method_label <- method_labels[all_df$tree1]
all_df$comparison   <- paste(all_df$method_label, "vs", all_df$tree2)
all_df$n_label      <- sprintf("n=%d", all_df$n_common)

p_nrf <- ggplot(all_df[!is.na(all_df$nRF), ],
                aes(x = reorder(comparison, nRF),
                    y = nRF,
                    fill = method_label)) +
  geom_col(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f\n%s", nRF, n_label)),
            hjust = -0.1, size = 2.5) +
  coord_flip() +
  scale_fill_manual(values = method_colors, name = "Method") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult=c(0,0.18))) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  labs(x = NULL, y = "Normalised RF distance") +
  theme_classic(base_size = 10) +
  theme(legend.position = "right")

p_coph <- ggplot(all_df[!is.na(all_df$cophenetic_r), ],
                 aes(x = reorder(comparison, cophenetic_r),
                     y = cophenetic_r,
                     fill = method_label)) +
  geom_col(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", cophenetic_r)),
            hjust = -0.1, size = 2.5) +
  coord_flip() +
  scale_fill_manual(values = method_colors, name = "Method") +
  scale_y_continuous(limits = c(-0.2, 1.2), expand = expansion(mult=c(0,0.1))) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  labs(x = NULL, y = "Cophenetic correlation (r)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "right")

combined <- p_nrf / p_coph

out_fig <- file.path(FIGS, "FigSX_709_vs_published.png")
ggsave(out_fig, combined, width = 12, height = 14, dpi = 300)
cat("Saved:", out_fig, "\n\n")

# ---- 6. 709 BS histogram ----
bs_csv <- tryCatch(
  read.csv(file.path(TREES, "709_bs_support_summary.csv"), stringsAsFactors=FALSE),
  error = function(e) NULL
)

cat("\n=== Key Results ===\n")
cat(sprintf("709 vs GPWG2:    nRF=%.3f, r=%.3f, n_common=%d\n",
    df_709["GPWG2_2012","nRF"],  df_709["GPWG2_2012","cophenetic_r"],
    df_709["GPWG2_2012","n_common"]))
cat(sprintf("709 vs Givnish:  nRF=%.3f, r=%.3f, n_common=%d\n",
    df_709["Givnish2018","nRF"], df_709["Givnish2018","cophenetic_r"],
    df_709["Givnish2018","n_common"]))
cat(sprintf("709 vs BK2015:   nRF=%.3f, r=%.3f, n_common=%d\n",
    df_709["BK2015","nRF"],      df_709["BK2015","cophenetic_r"],
    df_709["BK2015","n_common"]))

cat("\n=== DONE compare_709_vs_published.R ===\n")
