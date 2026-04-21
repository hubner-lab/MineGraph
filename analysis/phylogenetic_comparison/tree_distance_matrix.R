#!/usr/bin/env Rscript
# =============================================================================
# A2: Pairwise Tree Distance Matrix
# Computes nRF, Quartet divergence, Clustering Information Distance (CID),
# and Baker's Gamma for all 15 pairwise combinations of 6 trees.
# Each pair is pruned to common taxa before comparison.
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(dendextend)
})

# Install TreeDist / Quartet if needed
for (pkg in c("TreeDist", "Quartet")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}
suppressPackageStartupMessages({
  library(TreeDist)
  library(Quartet)
})

BASE  <- "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
TREES <- file.path(BASE, "phylogenetic_trees")
FIGS  <- file.path(BASE, "manuscript/figures")

# ---- 1. Load trees ----
tree_files <- c(
  MineGraph_BS = "192_graph_minegraph_bs_tree_named.nwk",
  PAV_Ward     = "192_graph_pav_tree_named.nwk",
  RAxML        = "192_graph_RAxML_tree_named.nwk",
  GPWG2_2012   = "GPWG2_2012_poaceae_renamed.nwk",
  Givnish2018  = "Givnish2018_monocot_renamed.nwk",
  BK2015       = "BK2015_MCC_POALES_renamed.nwk"
)

trees <- lapply(tree_files, function(f) {
  t <- read.tree(file.path(TREES, f))
  # Strip bootstrap labels (internal node labels) to avoid confusion
  t$node.label <- NULL
  t
})

tree_names <- names(trees)
n <- length(trees)
cat("Loaded", n, "trees\n")
for (nm in tree_names) cat(sprintf("  %s: %d tips\n", nm, length(trees[[nm]]$tip.label)))

# ---- 2. Helper: prune pair to common taxa ----
prune_pair <- function(t1, t2) {
  common <- intersect(t1$tip.label, t2$tip.label)
  if (length(common) < 4) stop("Too few common taxa (<4)")
  t1p <- keep.tip(t1, common)
  t2p <- keep.tip(t2, common)
  # Make unrooted for RF/Quartet/CID
  t1p <- unroot(t1p)
  t2p <- unroot(t2p)
  list(t1 = t1p, t2 = t2p, n = length(common))
}

# ---- 3. Compute pairwise metrics ----
results <- data.frame()

for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    nm1 <- tree_names[i]
    nm2 <- tree_names[j]
    cat(sprintf("Computing: %s vs %s ...\n", nm1, nm2))

    pair  <- prune_pair(trees[[i]], trees[[j]])
    t1    <- pair$t1
    t2    <- pair$t2
    n_com <- pair$n

    # --- nRF (normalized Robinson-Foulds) ---
    nRF <- tryCatch({
      rf <- RF.dist(t1, t2, normalize = TRUE, rooted = FALSE)
      round(as.numeric(rf), 4)
    }, error = function(e) NA_real_)

    # --- Clustering Information Distance (CID) — info-theoretic RF ---
    CID <- tryCatch({
      round(as.numeric(ClusteringInfoDistance(t1, t2, normalize = TRUE)), 4)
    }, error = function(e) NA_real_)

    # --- Quartet divergence (1 - similarity) ---
    QD <- tryCatch({
      sts <- QuartetStatus(list(t1, t2))
      qd  <- QuartetDivergence(sts)
      round(1 - qd[["similarity"]], 4)
    }, error = function(e) { cat("  Quartet error:", conditionMessage(e), "\n"); NA_real_ })

    # --- Cophenetic correlation (Baker's Gamma proxy, no dendrogram needed) ---
    BG <- tryCatch({
      c1 <- as.vector(cophenetic(t1))
      c2 <- as.vector(cophenetic(t2))
      round(cor(c1, c2, method = "pearson"), 4)
    }, error = function(e) { cat("  Cophenetic error:", conditionMessage(e), "\n"); NA_real_ })

    results <- rbind(results, data.frame(
      tree1       = nm1,
      tree2       = nm2,
      n_common    = n_com,
      nRF         = nRF,
      CID         = CID,
      Quartet_div = QD,
      Cophenetic_r = BG,
      stringsAsFactors = FALSE
    ))
  }
}

cat("\n--- Pairwise Distance Results ---\n")
print(results, digits = 4)

out_csv <- file.path(TREES, "pairwise_tree_distances.csv")
write.csv(results, out_csv, row.names = FALSE)
cat("\nSaved:", out_csv, "\n")

# ---- 4. Build symmetric matrices and heatmaps ----
library(reshape2)

# For visualization, use only nRF and CID (most standard)
for (metric in c("nRF", "CID", "Quartet_div", "Cophenetic_r")) {

  # Build symmetric matrix (higher = more different, except Cophenetic_r where higher = more similar)
  mat <- matrix(NA, n, n, dimnames = list(tree_names, tree_names))
  diag(mat) <- if (metric == "Cophenetic_r") 1 else 0

  for (k in 1:nrow(results)) {
    i_idx <- which(tree_names == results$tree1[k])
    j_idx <- which(tree_names == results$tree2[k])
    val   <- results[[metric]][k]
    mat[i_idx, j_idx] <- val
    mat[j_idx, i_idx] <- val
  }

  # Save matrix
  mat_out <- file.path(TREES, paste0("tree_distance_matrix_", metric, ".csv"))
  write.csv(round(mat, 4), mat_out)
  cat("Saved matrix:", mat_out, "\n")
}

# ---- 5. Heatmap figure (nRF + CID side by side) ----
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(ggplot2))

# Pretty labels
pretty_names <- c(
  MineGraph_BS = "Graph-PAV\n(MineGraph BS)",
  PAV_Ward     = "Graph-PAV\n(Ward)",
  RAxML        = "Graph\nRAxML",
  GPWG2_2012   = "GPWG2\n2012",
  Givnish2018  = "Givnish\n2018",
  BK2015       = "BK\n2015"
)

make_heatmap_df <- function(metric_name) {
  mat <- matrix(NA, n, n, dimnames = list(tree_names, tree_names))
  diag(mat) <- if (metric_name == "Cophenetic_r") 1 else 0
  for (k in 1:nrow(results)) {
    i_idx <- which(tree_names == results$tree1[k])
    j_idx <- which(tree_names == results$tree2[k])
    val   <- results[[metric_name]][k]
    mat[i_idx, j_idx] <- val
    mat[j_idx, i_idx] <- val
  }
  df <- melt(mat, varnames = c("Tree1", "Tree2"), value.name = "value")
  df$Tree1 <- factor(pretty_names[as.character(df$Tree1)],
                     levels = rev(pretty_names))
  df$Tree2 <- factor(pretty_names[as.character(df$Tree2)],
                     levels = pretty_names)
  df$metric <- metric_name
  df
}

# Combine nRF and CID
metric_labels <- c(
  nRF         = "Normalized Robinson-Foulds",
  CID         = "Clustering Information Distance",
  Quartet_div = "Quartet Divergence",
  Cophenetic_r= "Baker's Gamma (similarity)"
)

plot_list <- lapply(c("nRF", "CID"), function(m) {
  df <- make_heatmap_df(m)
  # For Baker's gamma, flip colour direction (higher = more similar)
  is_sim <- (m == "Cophenetic_r")
  ggplot(df, aes(x = Tree2, y = Tree1, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(is.na(value), "", sprintf("%.2f", value))),
              size = 2.8) +
    scale_fill_gradient2(
      low      = if (is_sim) "#D73027" else "#1A9850",
      mid      = "white",
      high     = if (is_sim) "#1A9850" else "#D73027",
      midpoint = if (is_sim) 0 else 0.5,
      na.value = "grey85",
      name     = if (is_sim) "Similarity" else "Distance",
      limits   = c(0, 1)
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(face = "bold", size = 10),
          legend.key.height = unit(0.8, "cm"))
})

if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages("patchwork", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(patchwork))

combined <- wrap_plots(plot_list, ncol = 2)

out_fig <- file.path(FIGS, "FigSX_tree_distance_heatmap.png")
ggsave(out_fig, combined, width = 13, height = 7, dpi = 300)
cat("Saved figure:", out_fig, "\n")

cat("\n=== DONE A2: Tree Distance Matrix ===\n")
