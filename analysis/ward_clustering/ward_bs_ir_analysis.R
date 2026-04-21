#!/usr/bin/env Rscript
# =============================================================================
# Ward vs Bootstrap Tree Comparison + Bootstrap × IR Family Analysis
#
# Two linked questions:
# Q1. Does the original Ward-linkage PAV tree agree with the new Felsenstein
#     bootstrap UPGMA tree (MineGraph, Jaccard distances, 1,000 BS replicates)?
#     (If yes: Ward topology was correct, just lacked stats.)
# Q2. Does bootstrap support at family-level nodes correlate with IR size?
#     (IR-expanded families: does higher IR content inflate apparent BS?)
#
# Outputs:
#   - ward_vs_bs_comparison.csv   (nRF, CID, cophenetic r: Ward vs BS)
#   - family_bs_ir_table.csv      (family, BS, sCF, IR_size_bp)
#   - FigSX_ward_vs_bs_tanglegram.png
#   - FigSX_bs_scf_ir_scatter.png
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(TreeDist)
})
if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages("patchwork", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(patchwork))

BASE  <- "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
TREES <- file.path(BASE, "phylogenetic_trees")
FIGS  <- file.path(BASE, "manuscript/figures")

cat("=== Ward vs Bootstrap + IR Analysis ===\n\n")

# ---- 1. Load trees ----
bs_tree   <- read.tree(file.path(TREES, "192_graph_minegraph_bs_tree_named.nwk"))
ward_tree <- read.tree(file.path(TREES, "192_graph_pav_tree_named.nwk"))

cat(sprintf("MineGraph BS tree:   %d tips\n", length(bs_tree$tip.label)))
cat(sprintf("PAV Ward tree:       %d tips\n", length(ward_tree$tip.label)))

# Both should have the same 192 tips
common_tips <- intersect(bs_tree$tip.label, ward_tree$tip.label)
cat(sprintf("Common tips:         %d\n\n", length(common_tips)))

# Prune to common tips
t_bs   <- drop.tip(bs_tree,   setdiff(bs_tree$tip.label,   common_tips))
t_ward <- drop.tip(ward_tree, setdiff(ward_tree$tip.label, common_tips))

# Strip node labels for distance functions
t_bs_stripped   <- t_bs;   t_bs_stripped$node.label   <- NULL
t_ward_stripped <- t_ward; t_ward_stripped$node.label <- NULL

# ---- 2. Pairwise metrics: Ward vs BS ----
cat("Computing Ward vs BS distance metrics ...\n")

nRF_val  <- RobinsonFoulds(t_bs_stripped, t_ward_stripped, normalize = TRUE)
CID_val  <- ClusteringInfoDistance(t_bs_stripped, t_ward_stripped)

# Cophenetic correlation
coph_bs   <- as.vector(cophenetic(t_bs_stripped))
coph_ward <- as.vector(cophenetic(t_ward_stripped))
coph_r    <- cor(coph_bs, coph_ward, method = "pearson")

# Topological concordance: what fraction of BS-tree bipartitions appear in Ward?
bs_parts   <- prop.part(t_bs_stripped)
ward_parts <- prop.part(t_ward_stripped)

# Convert to sorted character strings for matching
parts_to_str <- function(parts, tips) {
  lapply(parts, function(p) paste(sort(tips[p]), collapse=","))
}
bs_str   <- parts_to_str(bs_parts,   t_bs_stripped$tip.label)
ward_str <- parts_to_str(ward_parts, t_ward_stripped$tip.label)

shared_bipartitions  <- sum(bs_str %in% ward_str)
pct_bs_in_ward       <- 100 * shared_bipartitions / length(bs_str)

cat(sprintf("  nRF (Ward vs BS):          %.3f\n", nRF_val))
cat(sprintf("  CID (Ward vs BS):          %.3f\n", CID_val))
cat(sprintf("  Cophenetic r (Ward vs BS): %.3f\n", coph_r))
cat(sprintf("  BS bipartitions in Ward:   %d / %d (%.1f%%)\n",
            shared_bipartitions, length(bs_str), pct_bs_in_ward))
cat("\n")

# Save comparison table
comp_df <- data.frame(
  comparison     = "Ward (PAV) vs MineGraph_BS",
  n_common_tips  = length(common_tips),
  nRF            = round(nRF_val, 3),
  CID            = round(CID_val, 3),
  cophenetic_r   = round(coph_r, 3),
  pct_bs_in_ward = round(pct_bs_in_ward, 1),
  stringsAsFactors = FALSE
)
out_comp <- file.path(TREES, "ward_vs_bs_comparison.csv")
write.csv(comp_df, out_comp, row.names = FALSE)
cat("Saved:", out_comp, "\n\n")

# ---- 3. Bootstrap support at family MRCA nodes ----
cat("Extracting bootstrap support at family-level MRCA nodes ...\n")

# Family definitions — genera present in the 192 dataset
families_genera <- list(
  Poaceae      = c("Achnatherum","Acidosasa","Aegilops","Aeluropus","Agrostis",
                   "Alloteropsis","Alopecurus","Andropogon","Aristida","Avena",
                   "Axonopus","Bambusa","Bouteloua","Brachiaria","Brachypodium",
                   "Cenchrus","Chikusichloa","Chusquea","Coix","Dendrocalamus",
                   "Deschampsia","Digitaria","Diheteropogon","Echinochloa",
                   "Eragrostis","Eriachne","Fargesia","Festuca","Guadua",
                   "Hordeum","Leersia","Leptaspis","Merostachys","Miscanthus",
                   "Neomicrocalamus","Neyraudia","Oryza","Panicum","Phragmites",
                   "Phyllostachys","Poa","Puelia","Saccharum","Sarga","Schizachyrium",
                   "Schizostachyum","Secale","Sorghum","Sporobolus","Stipa",
                   "Streptochaeta","Streptogyna","Taeniatherum","Temochloa",
                   "Triticum","Zea","Gynerium","Lygeum","Nardus","Amphipogon",
                   "Agropyron","Agenium","Elymus","Leymus","Alloeochaete",
                   "Anomochloa","Sacciolepis","Stipagrostis","Phleum","Holcus",
                   "Beckmannia","Glyceria","Melocalamus","Rhipidocladum","Guaduella",
                   "Sinobambusa","Ampelocalamus","Garnotia","Joinvillea"),
  Cyperaceae   = c("Bolboschoenus","Carex","Cyperus","Eleocharis","Gahnia",
                   "Schoenoplectus","Schoenus","Scleria","Uncinia"),
  Bromeliaceae = c("Alcantarea","Ananas","Brocchinia","Guzmania","Hechtia",
                   "Lindmania","Navia","Neoregelia","Pitcairnia","Puya",
                   "Racinaea","Stigmatodon","Tillandsia","TPA_asm_Racinaea",
                   "Viridantha","Vriesea","Werauhia","Zizkaea","Pseudalcantarea"),
  Typhaceae    = c("Sparganium","Typha"),
  Juncaceae    = c("Juncus","Luzula"),
  Eriocaulaceae= c("Eriocaulon","Paepalanthus"),
  Xyridaceae   = c("Xyris"),
  Flagellariaceae = c("Flagellaria"),
  Joinvilleaceae  = c("Joinvillea"),
  Restionaceae    = c("Anarthria","Elegia","Restio","Ecdeiocolea"),
  Rapateaceae     = c("Rapatea","Navia")
)

# IR size data per family (from IR_GRAPH_ANALYSIS.md, mean of observed range)
ir_size_bp <- c(
  Poaceae         = 21614,   # mean of 20506–22722 bp
  Cyperaceae      = 37385,   # mean of 37276–37494 bp
  Bromeliaceae    = 27112,   # mean of 26708–27516 bp
  Typhaceae       = 26925,
  Juncaceae       = 23000,   # estimated from literature
  Eriocaulaceae   = 26437,
  Xyridaceae      = 17405,
  Flagellariaceae = 23500,   # estimated
  Joinvilleaceae  = 24048,
  Restionaceae    = 33752,
  Rapateaceae     = 25000    # estimated
)

# Helper: get bootstrap at MRCA
get_bs_at_mrca <- function(tree, tip_list) {
  tips_in_tree <- intersect(tip_list, tree$tip.label)
  if (length(tips_in_tree) < 2) return(list(bs = NA_real_, n = length(tips_in_tree)))
  node <- getMRCA(tree, tips_in_tree)
  if (is.null(node)) return(list(bs = NA_real_, n = length(tips_in_tree)))
  idx <- node - length(tree$tip.label)
  if (idx < 1 || idx > length(tree$node.label))
    return(list(bs = NA_real_, n = length(tips_in_tree)))
  val <- suppressWarnings(as.numeric(tree$node.label[idx]))
  list(bs = val, n = length(tips_in_tree))
}

# Get all tips matching each family's genera
get_family_tips <- function(tree, genera_list) {
  all_tips <- tree$tip.label
  matched  <- c()
  for (gen in genera_list) {
    pat <- paste0("^", gen, "_")
    matched <- c(matched, all_tips[grepl(pat, all_tips)])
  }
  unique(matched)
}

# Also load sCF branch table to get sCF at family MRCA
scf_branch <- tryCatch(
  read.csv(file.path(TREES, "scf_branch_table.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL
)
# Load sCF tree (with branch IDs) for ID lookup
scf_tree <- tryCatch(
  read.tree(file.path(TREES, "pav_scf_minegraph_bs.cf.tree")),
  error = function(e) NULL
)

results <- list()
for (fam_name in names(families_genera)) {
  genera    <- families_genera[[fam_name]]
  fam_tips  <- get_family_tips(bs_tree, genera)
  bs_result <- get_bs_at_mrca(bs_tree, fam_tips)
  ward_tips <- get_family_tips(ward_tree, genera)
  mono_ward <- if (length(ward_tips) >= 2)
    is.monophyletic(ward_tree, ward_tips) else NA

  # sCF at same MRCA node — find branch ID via sCF tree
  scf_val <- NA_real_
  if (!is.null(scf_tree) && !is.null(scf_branch) && length(fam_tips) >= 2) {
    tryCatch({
      # Map tips to scf_tree labels (same tip labels as bs_tree)
      scf_tips <- intersect(fam_tips, scf_tree$tip.label)
      if (length(scf_tips) >= 2) {
        scf_node <- getMRCA(scf_tree, scf_tips)
        if (!is.null(scf_node)) {
          # Node ID in scf_tree corresponds to ID column in branch table
          scf_idx <- scf_node - length(scf_tree$tip.label)
          branch_row <- scf_branch[scf_branch$ID == scf_node, ]
          if (nrow(branch_row) > 0) scf_val <- branch_row$sCF_MG[1]
        }
      }
    }, error = function(e) {})
  }

  results[[fam_name]] <- data.frame(
    family       = fam_name,
    n_tips       = bs_result$n,
    BS_support   = bs_result$bs,
    sCF_MG       = scf_val,
    ward_monophyletic = mono_ward,
    ir_size_bp   = ir_size_bp[fam_name],
    stringsAsFactors = FALSE
  )
}

family_df <- do.call(rbind, results)
family_df <- family_df[!is.na(family_df$BS_support), ]
cat(sprintf("  Families with BS values: %d\n", nrow(family_df)))
print(family_df[, c("family","n_tips","BS_support","sCF_MG","ward_monophyletic","ir_size_bp")])

out_fam <- file.path(TREES, "family_bs_ir_table.csv")
write.csv(family_df, out_fam, row.names = FALSE)
cat("Saved:", out_fam, "\n\n")

# ---- 4. Figure A: Cophenetic distance scatter (Ward vs BS) ----
cat("Producing Ward vs BS cophenetic distance scatter ...\n")

coph_df <- data.frame(
  ward = coph_ward,
  bs   = coph_bs
)
# Subsample to 5000 pairs for plotting speed
set.seed(42)
if (nrow(coph_df) > 5000) coph_df <- coph_df[sample(nrow(coph_df), 5000), ]

p_coph <- ggplot(coph_df, aes(x = ward, y = bs)) +
  geom_point(alpha = 0.08, size = 0.8, colour = "#2166AC") +
  geom_smooth(method = "lm", colour = "#D73027", linewidth = 0.9, se = FALSE) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = sprintf("r = %.3f\nnRF = %.3f", coph_r, nRF_val),
           size = 3.5, fontface = "italic") +
  labs(
    x = "Cophenetic distance (Ward-linkage PAV tree, Euclidean)",
    y = "Cophenetic distance (MineGraph UPGMA tree, Jaccard, 1000 Felsenstein BS)"
  ) +
  theme_classic(base_size = 11)

# ---- 5. Figure B: BS and sCF vs IR size per family ----
cat("Producing BS / sCF vs IR size scatter ...\n")

fam_plot <- family_df[!is.na(family_df$ir_size_bp), ]

# Long format for dual-panel
bs_long  <- fam_plot[!is.na(fam_plot$BS_support),  ]
bs_long$metric <- "Bootstrap support (BS)"
bs_long$value  <- bs_long$BS_support

scf_long <- fam_plot[!is.na(fam_plot$sCF_MG), ]
scf_long$metric <- "PAV site concordance (sCF)"
scf_long$value  <- scf_long$sCF_MG

plot_df <- rbind(
  bs_long[,  c("family","ir_size_bp","metric","value","ward_monophyletic","n_tips")],
  scf_long[, c("family","ir_size_bp","metric","value","ward_monophyletic","n_tips")]
)

p_ir <- ggplot(plot_df, aes(x = ir_size_bp / 1000, y = value,
                             colour = metric, label = family)) +
  geom_point(aes(shape = ward_monophyletic), size = 3) +
  geom_smooth(aes(group = metric), method = "lm", se = TRUE, alpha = 0.12,
              linewidth = 0.8) +
  geom_text(size = 2.5, vjust = -0.8, show.legend = FALSE) +
  scale_colour_manual(values = c("Bootstrap support (BS)"    = "#0072B2",
                                  "PAV site concordance (sCF)" = "#E69F00"),
                      name = "Support metric") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 17, "NA" = 15),
                     na.value = 15, name = "Ward monophyletic") +
  geom_hline(yintercept = 33.3, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  annotate("text", x = Inf, y = 35, label = "sCF random threshold (33.3%)",
           hjust = 1.1, size = 2.8, colour = "grey40") +
  labs(
    x = "Mean inverted repeat (IR) size (kb)",
    y = "Support value (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# ---- 6. Save figures ----
combined <- p_coph / p_ir

out_fig <- file.path(FIGS, "FigSX_ward_bs_ir_analysis.png")
ggsave(out_fig, combined, width = 10, height = 14, dpi = 300)
cat("Saved figure:", out_fig, "\n\n")

# ---- 7. Print interpretation summary ----
cat("=== INTERPRETATION SUMMARY ===\n")
cat(sprintf("Ward vs BS nRF = %.3f — %s topology overlap\n",
    nRF_val,
    ifelse(nRF_val < 0.2, "HIGH", ifelse(nRF_val < 0.4, "MODERATE", "LOW"))))
cat(sprintf("Cophenetic r = %.3f — divergence order is %s preserved\n",
    coph_r,
    ifelse(coph_r > 0.95, "very well", ifelse(coph_r > 0.85, "well", "moderately"))))
cat(sprintf("%.1f%% of BS-tree bipartitions present in Ward tree\n", pct_bs_in_ward))
cat("\nFamily BS support range:\n")
cat(sprintf("  Max: %s (BS=%.0f, IR=%.0f kb)\n",
    family_df$family[which.max(family_df$BS_support)],
    max(family_df$BS_support, na.rm=TRUE),
    family_df$ir_size_bp[which.max(family_df$BS_support)] / 1000))
cat(sprintf("  Min: %s (BS=%.0f, IR=%.0f kb)\n",
    family_df$family[which.min(family_df$BS_support)],
    min(family_df$BS_support, na.rm=TRUE),
    family_df$ir_size_bp[which.min(family_df$BS_support)] / 1000))

bs_ir_cor <- cor(family_df$BS_support, family_df$ir_size_bp, use = "complete.obs")
cat(sprintf("\nCorrelation BS vs IR size: r = %.3f\n", bs_ir_cor))

cat("\n=== DONE ward_bs_ir_analysis.R ===\n")
