#!/usr/bin/env Rscript
# =============================================================================
# A3: Clade Recovery Analysis
# Tests monophyly of 20 established Poales clades across all graph trees and
# published reference trees. Outputs recovery table and summary figure.
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(reshape2)
})

BASE  <- "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
TREES <- file.path(BASE, "phylogenetic_trees")
FIGS  <- file.path(BASE, "manuscript/figures")

# ============================================================
# 1. Load trees (strip node labels for clean monophyly test)
# ============================================================
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
  t$node.label <- NULL
  t
})

get_tips <- function(t) t$tip.label

# Helper: get bootstrap value at MRCA of a clade in MineGraph_BS tree
bs_tree <- read.tree(file.path(TREES, "192_graph_minegraph_bs_tree_named.nwk"))

get_bs_at_clade <- function(tree_with_bs, tips_present) {
  tryCatch({
    node <- getMRCA(tree_with_bs, tips_present)
    if (is.null(node)) return(NA_real_)
    idx  <- node - length(tree_with_bs$tip.label)  # offset into node.label
    val  <- suppressWarnings(as.numeric(tree_with_bs$node.label[idx]))
    val
  }, error = function(e) NA_real_)
}

# ============================================================
# 2. Define 20 target clades
#    Each entry: name, genera (genus prefix) or direct tip match
# ============================================================
# Genus-to-clade mapping based on APG IV / GPWG II classification
genus_to_clade <- list(
  # ---- Families ----
  Poaceae = c(
    "Achnatherum","Acidosasa","Acrachne","Aegilops","Aeluropus","Agenium",
    "Agropyron","Agrostis","Alcantarea","Alloeochaete","Alloteropsis","Alopecurus",
    "Amphipogon","Anadelphia","Anatherum","Andropogon","Andropterum","Anomochloa",
    "Apluda","Aristida","Arthraxon","Arundinella","Astrebla","Avena","Axonopus",
    "Bambusa","Barfussia","Beckmannia","Bolboschoenus","Bothriochloa","Bouteloua",
    "Brachiaria","Brachypodium","Briza","Bromus","Cenchrus","Chasmanthium",
    "Chikusichloa","Chionachne","Chrysopogon","Chusquea","Cleistogenes",
    "Coelachne","Coix","Coleanthus","Crinipes","Ctenium","Cynosurus",
    "Dactylis","Danthoniopsis","Dendrocalamus","Deschampsia","Desmostachya",
    "Dichanthelium","Digitaria","Diheteropogon","Echinochloa","Ehrharta",
    "Eleocharis","Elionurus","Elymus","Enneapogon","Eragrostis","Eriachne",
    "Eremochloa","Fargesia","Festuca","Froesiochloa","Garnotia","Germainia",
    "Glyceria","Guadua","Guaduella","Gynerium","Hakonechloa","Heteropogon",
    "Hitchcockella","Holcus","Homolepis","Hordeum","Humbertochloa",
    "Hyparrhenia","Joinvillea","Kampochloa","Lecomtella","Leersia","Leptagrostis",
    "Leymus","Littledalea","Loxodera","Lygeum","Melocalamus","Merostachys",
    "Merxmuellera","Miscanthus","Muhlenbergia","Nardus","Neomicrocalamus",
    "Neurachne","Neyraudia","Oryza","Otachyrium","Panicum","Perotis",
    "Phaenosperma","Phragmites","Phyllostachys","Phleum","Poa","Puelia",
    "Reynaudia","Rhipidocladum","Saccharum","Sacciolepis","Sarga","Schizachyrium",
    "Schizostachyum","Secale","Sinobambusa","Sorghum","Sporobolus","Stipa",
    "Stipagrostis","Streptochaeta","Streptogyna","Taeniatherum","Temochloa",
    "Themeda","Tragus","Trikeraia","Trichoneura","Triodia","Tripsacum",
    "Triticum","Uniola","Urochondra","Vaseyochloa","Yushania","Zea","Zeugites",
    "x_Triticosecale","Leptaspis"
  ),
  Bromeliaceae = c(
    "Aechmea","Alcantarea","Ananas","Bakerantha","Barfussia","Brocchinia",
    "Catopsis","Goudaea","Lindmania","Navia","Neoregelia","Pitcairnia",
    "Puya","Tillandsia","Vriesea"
  ),
  Cyperaceae = c(
    "Bolboschoenus","Carex","Cyperus","Eleocharis","Gahnia","Hypolytrum",
    "Isolepis","Scleria"
  ),
  Typhaceae = c("Rohrbachia","Sparganium","Typha"),
  Joinvilleaceae = c("Joinvillea"),
  Eriocaulaceae = c("Eriocaulon"),
  # ---- Poaceae subfamilies ----
  Bambusoideae = c(
    "Acidosasa","Bambusa","Chusquea","Dendrocalamus","Fargesia","Froesiochloa",
    "Guadua","Guaduella","Hitchcockella","Humbertochloa","Melocalamus",
    "Merostachys","Neomicrocalamus","Phyllostachys","Puelia","Rhipidocladum",
    "Schizostachyum","Sinobambusa","Streptogyna","Temochloa","Yushania"
  ),
  Oryzoideae = c(
    "Chikusichloa","Ehrharta","Leersia","Oryza"
  ),
  Pooideae = c(
    "Achnatherum","Agropyron","Agrostis","Alopecurus","Avena","Beckmannia",
    "Brachypodium","Briza","Bromus","Catapodium","Coleanthus","Cynosurus",
    "Dactylis","Deschampsia","Elymus","Festuca","Glyceria","Holcus",
    "Hordeum","Kampochloa","Leptagrostis","Leymus","Littledalea","Lygeum",
    "Nardus","Phaenosperma","Phleum","Poa","Secale","Stipa","Stipagrostis",
    "Taeniatherum","Trikeraia","Triticum","x_Triticosecale"
  ),
  BOP_clade = c(  # Bambusoideae + Oryzoideae + Pooideae
    "Acidosasa","Bambusa","Chusquea","Dendrocalamus","Fargesia","Froesiochloa",
    "Guadua","Guaduella","Hitchcockella","Humbertochloa","Melocalamus",
    "Merostachys","Neomicrocalamus","Phyllostachys","Puelia","Rhipidocladum",
    "Schizostachyum","Sinobambusa","Streptogyna","Temochloa","Yushania",
    "Chikusichloa","Ehrharta","Leersia","Oryza",
    "Achnatherum","Agropyron","Agrostis","Alopecurus","Avena","Beckmannia",
    "Brachypodium","Briza","Bromus","Catapodium","Coleanthus","Cynosurus",
    "Dactylis","Deschampsia","Elymus","Festuca","Glyceria","Holcus",
    "Hordeum","Kampochloa","Leptagrostis","Leymus","Littledalea","Lygeum",
    "Nardus","Phaenosperma","Phleum","Poa","Secale","Stipa","Stipagrostis",
    "Taeniatherum","Trikeraia","Triticum","x_Triticosecale"
  ),
  Panicoideae = c(
    "Alloeochaete","Alloteropsis","Andropogon","Andropterum","Anadelphia",
    "Anatherum","Apluda","Arthraxon","Arundinella","Bothriochloa","Brachiaria",
    "Cenchrus","Chasmanthium","Chionachne","Chrysopogon","Coix","Dichanthelium",
    "Diheteropogon","Digitaria","Echinochloa","Elionurus","Eremochloa",
    "Garnotia","Germainia","Gynerium","Heteropogon","Homolepis","Hyparrhenia",
    "Lecomtella","Loxodera","Miscanthus","Neurachne","Otachyrium","Panicum",
    "Reynaudia","Saccharum","Sacciolepis","Sarga","Schizachyrium","Sorghum",
    "Themeda","Tragus","Tripsacum","Uniola","Zeugites","Zea","Axonopus"
  ),
  Chloridoideae = c(
    "Acrachne","Aeluropus","Astrebla","Bouteloua","Cleistogenes","Ctenium",
    "Dactylis","Desmostachya","Eleocharis","Enneapogon","Eragrostis",
    "Muhlenbergia","Neyraudia","Perotis","Sporobolus","Trichoneura","Triodia",
    "Urochondra","Vaseyochloa"
  ),
  Aristidoideae = c("Aristida","Stipagrostis"),
  Danthonioideae = c("Amphipogon","Danthoniopsis","Merxmuellera",
                     "Crinipes","Leptagrostis","Hakonechloa","Phragmites"),
  Micrairoideae = c("Coelachne","Eriachne"),
  PACMAD_clade = c(  # all PACMAD subfamilies
    "Acrachne","Aeluropus","Alloeochaete","Alloteropsis","Amphipogon",
    "Andropogon","Andropterum","Anadelphia","Anatherum","Apluda","Aristida",
    "Arthraxon","Arundinella","Astrebla","Axonopus","Bothriochloa","Bouteloua",
    "Brachiaria","Cenchrus","Chasmanthium","Chionachne","Chrysopogon","Coelachne",
    "Coix","Crinipes","Ctenium","Danthoniopsis","Desmostachya","Dichanthelium",
    "Diheteropogon","Digitaria","Echinochloa","Elionurus","Enneapogon","Eragrostis",
    "Eriachne","Eremochloa","Garnotia","Germainia","Gynerium","Hakonechloa",
    "Heteropogon","Homolepis","Hyparrhenia","Lecomtella","Leptagrostis","Loxodera",
    "Merxmuellera","Miscanthus","Muhlenbergia","Neurachne","Neyraudia","Otachyrium",
    "Panicum","Perotis","Phragmites","Reynaudia","Saccharum","Sacciolepis","Sarga",
    "Schizachyrium","Sorghum","Sporobolus","Stipagrostis","Themeda","Trichoneura",
    "Triodia","Tripsacum","Uniola","Urochondra","Vaseyochloa","Zeugites","Zea",
    "Aeluropus","Anadelphia","Amphipogon","Alloeochaete"
  ),
  # ---- Broader Poales groupings ----
  Commelinids_core = c(  # non-grass Poales: Bromeliaceae+Typhaceae+Eriocaulaceae+Xyridaceae+Flagellariaceae+etc
    "Aechmea","Alcantarea","Ananas","Bakerantha","Barfussia","Brocchinia",
    "Catopsis","Eriocaulon","Flagellaria","Goudaea","Lindmania","Navia",
    "Neoregelia","Pitcairnia","Puya","Rapatea","Rohrbachia","Sparganium",
    "Tillandsia","Typha","Vriesea","Xyris"
  ),
  Anomochlooideae_Pharoideae = c("Anomochloa","Leptaspis","Streptochaeta")
)

# ============================================================
# 3. Match genera to actual tips in each tree
# ============================================================
match_clade_to_tips <- function(tree, genera_list) {
  all_tips <- tree$tip.label
  matched  <- unlist(lapply(genera_list, function(g) {
    all_tips[startsWith(all_tips, paste0(g, "_"))]
  }))
  unique(matched)
}

# ============================================================
# 4. Test monophyly per clade per tree
# ============================================================
# Returns: "Y" (monophyletic), "P" (paraphyletic/polyphyletic, ≥50% members in one clade),
#          "N" (not recovered), "X" (insufficient taxa, <2 in tree)
test_monophyly <- function(tree, tips_in_clade) {
  avail <- intersect(tips_in_clade, tree$tip.label)
  n_avail <- length(avail)
  n_total <- length(tips_in_clade)
  if (n_avail < 2) return(list(status = "X", n = n_avail, n_total = n_total))
  mono <- is.monophyletic(tree, avail)
  if (mono) {
    status <- "Y"
  } else {
    # Check if ≥50% clustered (partial recovery)
    # Find largest monophyletic subset using MRCA node
    node <- getMRCA(tree, avail)
    clade_tips <- extract.clade(tree, node)$tip.label
    overlap <- length(intersect(avail, clade_tips)) / n_avail
    status <- if (overlap >= 0.75) "P+" else if (overlap >= 0.50) "P" else "N"
  }
  list(status = status, n = n_avail, n_total = n_total)
}

# ============================================================
# 5. Run all tests
# ============================================================
results_list <- list()

for (clade_nm in names(genus_to_clade)) {
  genera <- genus_to_clade[[clade_nm]]
  row    <- list(Clade = clade_nm)

  # BS value for this clade in MineGraph_BS tree
  bs_tips <- match_clade_to_tips(bs_tree, genera)
  bs_val  <- if (length(bs_tips) >= 2) get_bs_at_clade(bs_tree, bs_tips) else NA_real_

  for (tr_nm in names(trees)) {
    clade_tips <- match_clade_to_tips(trees[[tr_nm]], genera)
    res        <- test_monophyly(trees[[tr_nm]], clade_tips)
    row[[tr_nm]] <- res$status
    if (tr_nm == "MineGraph_BS") {
      row[["n_in_MineGraph"]]  <- res$n
      row[["n_total_defined"]] <- res$n_total
      row[["BS_support"]]      <- round(bs_val, 1)
    }
  }
  results_list[[clade_nm]] <- as.data.frame(row, stringsAsFactors = FALSE)
}

results <- do.call(rbind, results_list)
rownames(results) <- NULL

# Reorder columns
col_order <- c("Clade","n_total_defined","n_in_MineGraph","BS_support",
               "MineGraph_BS","PAV_Ward","RAxML","GPWG2_2012","Givnish2018","BK2015")
results <- results[, intersect(col_order, colnames(results))]

cat("\n--- Clade Recovery Table ---\n")
print(results)

out_csv <- file.path(TREES, "clade_recovery_table.csv")
write.csv(results, out_csv, row.names = FALSE)
cat("\nSaved:", out_csv, "\n")

# ============================================================
# 6. Summary statistics
# ============================================================
status_to_score <- function(s) {
  # Y=full, P+=high partial, P=partial, N=not recovered, X=not tested
  c("Y"=1, "P+"=0.75, "P"=0.5, "N"=0, "X"=NA)[s]
}

cat("\n--- Recovery rates (excluding X=insufficient taxa) ---\n")
for (tr in c("MineGraph_BS","PAV_Ward","RAxML","GPWG2_2012","Givnish2018","BK2015")) {
  if (!(tr %in% colnames(results))) next
  vals   <- results[[tr]]
  tested <- vals[vals != "X"]
  full_Y <- mean(tested == "Y") * 100
  any_Y_P<- mean(tested %in% c("Y","P+","P")) * 100
  cat(sprintf("  %-20s  Y: %4.1f%%  Y+P: %4.1f%%  (n tested = %d)\n",
              tr, full_Y, any_Y_P, length(tested)))
}

# ============================================================
# 7. Bar chart figure
# ============================================================
tree_labels <- c(
  MineGraph_BS = "Graph-PAV\n(MineGraph BS)",
  PAV_Ward     = "Graph-PAV\n(Ward)",
  RAxML        = "Graph\nRAxML",
  GPWG2_2012   = "GPWG2 2012",
  Givnish2018  = "Givnish 2018",
  BK2015       = "BK 2015"
)

tree_cols <- c("MineGraph_BS","PAV_Ward","RAxML","GPWG2_2012","Givnish2018","BK2015")
plot_data <- data.frame()

for (tr in tree_cols) {
  if (!(tr %in% colnames(results))) next
  vals   <- results[[tr]]
  tested <- vals[vals != "X"]
  n      <- length(tested)
  plot_data <- rbind(plot_data, data.frame(
    Tree     = tree_labels[tr],
    Status   = c("Full (Y)","Partial (P+/P)","Not recovered (N)"),
    Count    = c(sum(tested == "Y"),
                 sum(tested %in% c("P+","P")),
                 sum(tested == "N")),
    Pct      = c(sum(tested == "Y") / n * 100,
                 sum(tested %in% c("P+","P")) / n * 100,
                 sum(tested == "N") / n * 100),
    stringsAsFactors = FALSE
  ))
}

plot_data$Tree   <- factor(plot_data$Tree,   levels = tree_labels)
plot_data$Status <- factor(plot_data$Status,
                            levels = c("Full (Y)","Partial (P+/P)","Not recovered (N)"))

pal <- c("Full (Y)"="#1A9850","Partial (P+/P)"="#FEE08B","Not recovered (N)"="#D73027")

p <- ggplot(plot_data, aes(x = Tree, y = Pct, fill = Status)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")),
            position = position_stack(vjust = 0.5),
            size = 3.2, colour = "black") +
  scale_fill_manual(values = pal, name = "Monophyly") +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  labs(
    x = NULL, y = "% of testable clades"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 9))

out_fig <- file.path(FIGS, "FigSX_clade_recovery_bars.png")
ggsave(out_fig, p, width = 10, height = 5.5, dpi = 300)
cat("Saved figure:", out_fig, "\n")

cat("\n=== DONE A3: Clade Recovery ===\n")
