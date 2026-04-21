#!/usr/bin/env Rscript
# =============================================================================
# A6: Full Pairwise Displacement Matrix
#
# Reviewer 2 specifically requested:
#   (a) Validation of the "normalized displacement" metric
#   (b) All pairwise comparisons (not just PAV vs RAxML)
#
# Normalized displacement is formally defined here:
#   For trees T1 and T2 with n common taxa:
#     1. Prune both trees to common taxa
#     2. Root (midpoint rooting if unrooted); ladderize
#     3. Extract leaf order r_k(i) for each taxon i in tree k
#     4. displacement_i = |r_1(i) - r_2(i)|           (absolute rank shift)
#     5. norm_displacement_i = displacement_i / (n-1)  (range: [0,1])
#   Tree-level summary: mean_norm_disp, max_norm_disp, pct_displaced_>20%
#   Baker's Gamma (rho): Pearson r of all pairwise cophenetic distances
#                        — measures overall divergence-order preservation
#   Entanglement: 1 - Baker's Gamma (dendextend::entanglement)
#
# Outputs:
#   - full_pairwise_displacement_matrix.csv    (all 15 pairs × metrics)
#   - taxon_displacement_summary.csv           (per-taxon displacement profile)
#   - universal_displacement_taxa.csv          (displaced in ≥4 of 5 comparisons)
#   - FigSX_displacement_matrix_heatmap.png
#   - FigSX_universal_displacement_taxa.png
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(dendextend)
})
if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages("patchwork", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(patchwork))
if (!requireNamespace("phytools", quietly = TRUE))
  install.packages("phytools", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(phytools))

BASE  <- "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
TREES <- file.path(BASE, "phylogenetic_trees")
FIGS  <- file.path(BASE, "manuscript/figures")

cat("=== A6: Full Pairwise Displacement Matrix ===\n\n")

# ---- 1. Load all trees ----
tree_files <- c(
  "MineGraph_BS"  = "192_graph_minegraph_bs_tree_named.nwk",
  "PAV_Ward"      = "192_graph_pav_tree_named.nwk",
  "RAxML"         = "192_graph_RAxML_tree_named.nwk",
  "GPWG2_2012"    = "GPWG2_2012_poaceae_renamed.nwk",
  "Givnish2018"   = "Givnish2018_monocot_renamed.nwk",
  "BK2015"        = "BK2015_MCC_POALES_renamed.nwk"
)

trees <- lapply(tree_files, function(f) {
  t <- read.tree(file.path(TREES, f))
  # Midpoint root if unrooted
  if (!is.rooted(t)) {
    t <- tryCatch(midpoint.root(t), error = function(e) {
      # Fallback: root at first tip
      root(t, outgroup = t$tip.label[1], resolve.root = TRUE)
    })
  }
  t$node.label <- NULL
  t
})

cat("Trees loaded:\n")
for (nm in names(trees)) cat(sprintf("  %-16s %d tips\n", nm, length(trees[[nm]]$tip.label)))

# ---- 2. Helper: tree → ladderized dendrogram ----
tree_to_dend <- function(t) {
  d   <- cophenetic(t)
  hcl <- hclust(as.dist(d), method = "average")
  dnd <- as.dendrogram(hcl)
  dnd <- dendextend::ladderize(dnd)
  dnd
}

# ---- 3. Leaf order from dendrogram ----
dend_leaf_order <- function(dnd) {
  labs <- labels(dnd)
  order <- seq_along(labs)
  names(order) <- labs
  order
}

# ---- 4. Core displacement function ----
# Returns per-taxon normalized displacement + tree-level metrics
compute_displacement <- function(t1, t2, label1, label2) {
  common <- intersect(t1$tip.label, t2$tip.label)
  n      <- length(common)
  if (n < 4) {
    return(list(
      summary = data.frame(
        tree1=label1, tree2=label2, n_common=n,
        mean_norm_disp=NA, median_norm_disp=NA, max_norm_disp=NA,
        pct_above_20pct=NA, bakers_gamma=NA, entanglement=NA,
        stringsAsFactors=FALSE
      ),
      per_taxon = NULL
    ))
  }

  # Prune BOTH trees to common taxa BEFORE constructing dendrograms
  # This ensures leaf ranks span 1:n_common, making normalization valid
  p1 <- drop.tip(t1, setdiff(t1$tip.label, common))
  p2 <- drop.tip(t2, setdiff(t2$tip.label, common))

  d1 <- tree_to_dend(p1)
  d2 <- tree_to_dend(p2)

  ord1 <- dend_leaf_order(d1)
  ord2 <- dend_leaf_order(d2)

  # Intersect common taxa with those actually in both pruned dendrograms
  taxa <- intersect(common, intersect(names(ord1), names(ord2)))
  if (length(taxa) < 4) {
    return(list(
      summary = data.frame(
        tree1=label1, tree2=label2, n_common=n,
        mean_norm_disp=NA, median_norm_disp=NA, max_norm_disp=NA,
        pct_above_20pct=NA, bakers_gamma=NA, entanglement=NA,
        stringsAsFactors=FALSE
      ),
      per_taxon = NULL
    ))
  }

  rank1 <- ord1[taxa]
  rank2 <- ord2[taxa]
  ok    <- !is.na(rank1) & !is.na(rank2)
  taxa  <- taxa[ok];  rank1 <- rank1[ok];  rank2 <- rank2[ok]
  n_eff <- length(taxa)

  # Re-rank both within 1:n_eff to guarantee valid normalization
  # This handles any residual rank-mismatch from cophenetic rounding
  rank1 <- rank(rank1, ties.method = "first")
  rank2 <- rank(rank2, ties.method = "first")

  # Absolute displacement (now guaranteed to be in 0 .. n_eff-1)
  disp      <- abs(rank1 - rank2)
  # Normalized displacement: / (n_eff - 1) → range [0, 1]
  norm_disp <- disp / max(1L, n_eff - 1L)

  # Summaries
  mean_nd   <- mean(norm_disp, na.rm = TRUE)
  median_nd <- median(norm_disp, na.rm = TRUE)
  max_nd    <- max(norm_disp, na.rm = TRUE)
  pct_above <- 100 * mean(norm_disp > 0.20, na.rm = TRUE)

  # Baker's Gamma via cophenetic correlation
  bakers_g <- tryCatch({
    coph1 <- as.vector(cophenetic(p1))
    coph2 <- as.vector(cophenetic(p2))
    cor(coph1, coph2, method = "pearson")
  }, error = function(e) NA_real_)

  # Entanglement (dendextend): 1 = perfect disorder, 0 = perfect agreement
  entangl <- tryCatch({
    dl <- dendextend::dendlist(d1, d2)
    dendextend::entanglement(dl)
  }, error = function(e) NA_real_)

  per_taxon <- data.frame(
    taxon        = taxa,
    rank_tree1   = as.integer(rank1),
    rank_tree2   = as.integer(rank2),
    displacement = as.integer(disp),
    norm_disp    = round(norm_disp, 4),
    tree1        = label1,
    tree2        = label2,
    n_common     = n,
    stringsAsFactors = FALSE
  )

  list(
    summary = data.frame(
      tree1             = label1,
      tree2             = label2,
      n_common          = n,
      mean_norm_disp    = round(mean_nd,   4),
      median_norm_disp  = round(median_nd, 4),
      max_norm_disp     = round(max_nd,    4),
      pct_above_20pct   = round(pct_above, 1),
      bakers_gamma      = round(bakers_g,  3),
      entanglement      = round(entangl,   3),
      stringsAsFactors  = FALSE
    ),
    per_taxon = per_taxon
  )
}

# ---- 5. Compute all 15 pairwise combinations ----
tree_names <- names(trees)
pairs      <- combn(tree_names, 2, simplify = FALSE)

cat(sprintf("\nComputing %d pairwise comparisons ...\n", length(pairs)))

all_summary  <- list()
all_per_taxon <- list()

for (pr in pairs) {
  t1_nm <- pr[1]; t2_nm <- pr[2]
  cat(sprintf("  %s vs %s\n", t1_nm, t2_nm))
  res <- compute_displacement(trees[[t1_nm]], trees[[t2_nm]], t1_nm, t2_nm)
  all_summary[[paste(t1_nm, t2_nm, sep="__")]]   <- res$summary
  if (!is.null(res$per_taxon))
    all_per_taxon[[paste(t1_nm, t2_nm, sep="__")]] <- res$per_taxon
}

summary_df  <- do.call(rbind, all_summary)
per_taxon_df <- do.call(rbind, all_per_taxon)

cat("\n--- Pairwise Displacement Summary ---\n")
print(summary_df[, c("tree1","tree2","n_common","mean_norm_disp","bakers_gamma","entanglement","pct_above_20pct")])

# Save
write.csv(summary_df, file.path(TREES, "full_pairwise_displacement_matrix.csv"), row.names=FALSE)
cat("Saved: full_pairwise_displacement_matrix.csv\n\n")

# ---- 6. Per-taxon profile: identify universally displaced taxa ----
# A taxon is "universally displaced" if it appears among the top 20% displaced
# taxa in at least 4 out of 5 comparisons that include graph trees (vs. published refs)
# Focus comparisons: graph trees (MineGraph_BS, PAV_Ward, RAxML) vs published (GPWG2, Givnish, BK2015)

focus_pairs <- summary_df[
  (summary_df$tree1 %in% c("MineGraph_BS","PAV_Ward","RAxML") &
   summary_df$tree2 %in% c("GPWG2_2012","Givnish2018","BK2015")) |
  (summary_df$tree2 %in% c("MineGraph_BS","PAV_Ward","RAxML") &
   summary_df$tree1 %in% c("GPWG2_2012","Givnish2018","BK2015")),
  "tree1"
]
# Get per_taxon rows for those pairs
focus_keys <- unique(c(
  paste(c("MineGraph_BS","PAV_Ward","RAxML"), "GPWG2_2012", sep="__"),
  paste(c("MineGraph_BS","PAV_Ward","RAxML"), "Givnish2018", sep="__"),
  paste(c("MineGraph_BS","PAV_Ward","RAxML"), "BK2015", sep="__"),
  paste("GPWG2_2012", c("MineGraph_BS","PAV_Ward","RAxML"), sep="__"),
  paste("Givnish2018", c("MineGraph_BS","PAV_Ward","RAxML"), sep="__"),
  paste("BK2015", c("MineGraph_BS","PAV_Ward","RAxML"), sep="__")
))
focus_keys <- focus_keys[focus_keys %in% names(all_per_taxon)]

focus_pt <- do.call(rbind, all_per_taxon[focus_keys])

# For each comparison, flag taxa in top 20% by norm_disp
flag_top20 <- function(df) {
  thresh <- quantile(df$norm_disp, 0.80, na.rm=TRUE)
  df$top20 <- df$norm_disp >= thresh
  df
}

focus_pt <- do.call(rbind, lapply(
  split(focus_pt, paste(focus_pt$tree1, focus_pt$tree2)),
  flag_top20
))

# Count how many comparisons each taxon is in the top-20%
taxon_top20_count <- tapply(focus_pt$top20, focus_pt$taxon, sum, na.rm=TRUE)
taxon_n_appear    <- tapply(focus_pt$taxon, focus_pt$taxon, length)
taxon_mean_disp   <- tapply(focus_pt$norm_disp, focus_pt$taxon, mean, na.rm=TRUE)

# Also compute SD of displacement across comparisons
taxon_sd_disp     <- tapply(focus_pt$norm_disp, focus_pt$taxon, sd, na.rm=TRUE)

taxon_summary <- data.frame(
  taxon              = names(taxon_top20_count),
  n_comparisons      = as.integer(taxon_n_appear[names(taxon_top20_count)]),
  n_top20_displaced  = as.integer(taxon_top20_count),
  mean_norm_disp     = round(taxon_mean_disp[names(taxon_top20_count)], 4),
  sd_norm_disp       = round(taxon_sd_disp[names(taxon_top20_count)], 4),
  stringsAsFactors   = FALSE
)
taxon_summary <- taxon_summary[order(-taxon_summary$n_top20_displaced, -taxon_summary$mean_norm_disp), ]

# Universal displacement: top-20% in ≥4 comparisons
universal <- taxon_summary[taxon_summary$n_top20_displaced >= 4, ]

cat(sprintf("Universal displacement taxa (top-20%% in ≥4 comparisons): %d\n", nrow(universal)))
if (nrow(universal) > 0) print(head(universal, 20))

# Selective displacement: high mean but few comparisons
selective <- taxon_summary[taxon_summary$n_comparisons >= 2 &
                            taxon_summary$mean_norm_disp > 0.30 &
                            taxon_summary$n_top20_displaced < 4, ]
cat(sprintf("Selective displacement taxa (high mean but < 4 comparisons): %d\n", nrow(selective)))

write.csv(taxon_summary, file.path(TREES, "taxon_displacement_summary.csv"), row.names=FALSE)
write.csv(universal,     file.path(TREES, "universal_displacement_taxa.csv"),   row.names=FALSE)
cat("Saved: taxon_displacement_summary.csv, universal_displacement_taxa.csv\n\n")

# ---- 7. Figure A: Displacement matrix heatmap ----
cat("Producing displacement heatmap ...\n")

# Build symmetric matrix for mean_norm_disp
all_labels <- unique(c(summary_df$tree1, summary_df$tree2))
mat_disp   <- matrix(NA, length(all_labels), length(all_labels),
                     dimnames = list(all_labels, all_labels))
mat_gamma  <- matrix(NA, length(all_labels), length(all_labels),
                     dimnames = list(all_labels, all_labels))

for (i in seq_len(nrow(summary_df))) {
  t1 <- summary_df$tree1[i]; t2 <- summary_df$tree2[i]
  mat_disp[t1, t2]  <- summary_df$mean_norm_disp[i]
  mat_disp[t2, t1]  <- summary_df$mean_norm_disp[i]
  mat_gamma[t1, t2] <- summary_df$bakers_gamma[i]
  mat_gamma[t2, t1] <- summary_df$bakers_gamma[i]
  diag(mat_disp)    <- 0
  diag(mat_gamma)   <- 1
}

# Melt for ggplot
melt_mat <- function(mat, value_name) {
  df <- as.data.frame(as.table(mat))
  names(df) <- c("tree1","tree2",value_name)
  df
}
hm_disp  <- melt_mat(mat_disp,  "mean_norm_disp")
hm_gamma <- melt_mat(mat_gamma, "bakers_gamma")

# Label order: graph trees first, then published
label_order <- c("MineGraph_BS","PAV_Ward","RAxML","GPWG2_2012","Givnish2018","BK2015")
hm_disp$tree1  <- factor(hm_disp$tree1,  levels = rev(label_order))
hm_disp$tree2  <- factor(hm_disp$tree2,  levels = label_order)
hm_gamma$tree1 <- factor(hm_gamma$tree1, levels = rev(label_order))
hm_gamma$tree2 <- factor(hm_gamma$tree2, levels = label_order)

p_hm_disp <- ggplot(hm_disp[!is.na(hm_disp$mean_norm_disp), ],
                    aes(x = tree2, y = tree1, fill = mean_norm_disp)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", mean_norm_disp)), size = 3) +
  scale_fill_gradient2(low = "#1A9850", mid = "#FFFFBF", high = "#D73027",
                       midpoint = 0.25, na.value = "grey80",
                       name = "Mean norm.\ndisplacement") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle=30, hjust=1))

p_hm_gamma <- ggplot(hm_gamma[!is.na(hm_gamma$bakers_gamma), ],
                     aes(x = tree2, y = tree1, fill = bakers_gamma)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", bakers_gamma)), size = 3) +
  scale_fill_gradient2(low = "#D73027", mid = "#FFFFBF", high = "#1A9850",
                       midpoint = 0.5, na.value = "grey80",
                       name = "Baker's\nGamma (r)") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle=30, hjust=1))

combined_hm <- p_hm_disp / p_hm_gamma

out_hm <- file.path(FIGS, "FigSX_displacement_matrix_heatmap.png")
ggsave(out_hm, combined_hm, width = 10, height = 12, dpi = 300)
cat("Saved:", out_hm, "\n\n")

# ---- 8. Figure B: Universal displacement taxa ----
if (nrow(universal) > 0) {
  top_taxa <- head(universal, 30)

  # Per-taxon displacement per comparison
  focus_top <- focus_pt[focus_pt$taxon %in% top_taxa$taxon, ]
  focus_top$comparison <- paste(focus_top$tree1, "vs", focus_top$tree2)

  p_univ <- ggplot(focus_top,
                   aes(x = reorder(taxon, -norm_disp, FUN=mean),
                       y = norm_disp,
                       fill = top20)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(colour = comparison), position = position_jitter(width=0.15),
               size = 1.2, alpha = 0.6) +
    scale_fill_manual(values = c("TRUE"="#D73027", "FALSE"="#91BFDB"),
                      name = "Top 20%\ndisplaced") +
    geom_hline(yintercept = 0.20, linetype="dashed", colour="grey40") +
    coord_flip() +
    labs(x = NULL, y = "Normalised displacement") +
    theme_classic(base_size = 9) +
    theme(legend.position = "right", legend.text = element_text(size=6))

  out_univ <- file.path(FIGS, "FigSX_universal_displacement_taxa.png")
  ggsave(out_univ, p_univ, width=12, height=max(5, 0.35*nrow(top_taxa)+2), dpi=300)
  cat("Saved:", out_univ, "\n\n")
} else {
  cat("No universally displaced taxa found (threshold: top-20% in ≥4 comparisons).\n\n")
  # Try relaxed threshold
  modest <- taxon_summary[taxon_summary$n_top20_displaced >= 2 &
                           taxon_summary$mean_norm_disp > 0.25, ]
  cat(sprintf("Relaxed threshold (≥2 comparisons, mean_nd>0.25): %d taxa\n", nrow(modest)))
  if (nrow(modest) > 0) {
    cat("Top taxa:\n"); print(head(modest, 10))
    write.csv(modest, file.path(TREES, "universal_displacement_taxa.csv"), row.names=FALSE)
    cat("Saved: universal_displacement_taxa.csv (relaxed threshold)\n\n")
  }
}

# ---- 9. Print formal definition for Methods section ----
cat("=== FORMAL DEFINITION (for manuscript Methods) ===\n")
cat("
Normalised displacement for taxon i between trees T1 and T2 is:

  nd_i = |r_1(i) - r_2(i)| / (n - 1)

where r_k(i) is the rank (1-indexed leaf order position) of taxon i in
the ladderized dendrogram of tree Tk (derived from cophenetic distance
matrix via average-linkage UPGMA), and n is the number of taxa common
to both trees. nd_i ∈ [0, 1]: nd_i = 0 means identical position in
both trees; nd_i = 1 means maximum possible rank reversal.

Summary metrics per pair:
  mean_norm_disp   = mean(nd_i) over all common taxa
  max_norm_disp    = max(nd_i)
  pct_above_20pct  = % taxa with nd_i > 0.20

Baker's Gamma (rho): Pearson correlation of all pairwise cophenetic
distances between T1 and T2 on common taxa. rho = 1 indicates
identical relative divergence order; rho = 0 indicates no relationship.

Universal displacement: a taxon is universally displaced if it ranks
in the top 20% highest-nd taxa in ≥4 of 9 graph-vs-published
pairwise comparisons. These represent taxa where structural variation
signal (graph topology) systematically disagrees with published
sequence-based phylogenies — candidate regions for phylogenetic
inconsistency driven by structural rather than sequence evolution.
")

cat("=== DONE A6: pairwise_displacement_matrix.R ===\n")
