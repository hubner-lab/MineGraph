#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(topGO)
  library(GO.db)
  library(AnnotationDbi)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
})

# ======================
# ====== CONFIG ========
# ======================
OBO_FILE     <- "go-basic.obo"          # your downloaded GO
ASSOC_FILE   <- "gene2go.assoc"         # "gene<TAB>GO:...;GO:..."
POP_FILE     <- "population.ids"        # background genes (IDs must match ASSOC_FILE)
CONS_FILE    <- "study_consensus.ids"   # consensus study list (subset of POP_FILE)
NONCONS_FILE <- "study_nonconsensus.ids" # non-consensus study list (subset of POP_FILE)

# Optional: provide a text file with GO IDs (one per line) to restrict to a GO-Slim
# For Plant slim, put its GO IDs into a file and set SLIM_FILE to that path; else "" to disable.
SLIM_FILE    <- ""  # e.g., "goslim_plant_ids.txt" or leave "" to use full GO

# Plot settings
TOP_N        <- 20       # top terms to plot by p-value (per namespace, per set, per algorithm)
ALPHA_FDR    <- 0.1      # show dashed line at this FDR in plots (visual aid)

# Output prefix
PREFIX       <- "topgo"

# ======================
# === LOAD MAPPINGS ====
# ======================
message("Reading inputs...")
stopifnot(file.exists(OBO_FILE), file.exists(ASSOC_FILE), file.exists(POP_FILE))


# gene2GO: list(gene -> character vector of GO IDs)
read_assoc <- function(path) {
  # Expect "gene<TAB>GO:...;GO:..." (as built earlier)
  df <- read_tsv(path, col_names = c("gene","gos"), show_col_types = FALSE)
  df <- df %>% mutate(
    gene = gene %>% as.character(),
    gos  = ifelse(is.na(gos), "", gos)
  )
  # split semicolon list
  split_list <- strsplit(df$gos, ";")
  names(split_list) <- df$gene
  # keep only valid GO IDs
  split_list <- lapply(split_list, function(v) unique(grep("^GO:\\d{7}$", v, value = TRUE)))
  # drop empty
  split_list[sapply(split_list, length) > 0]
}

gene2GO_all <- read_assoc(ASSOC_FILE)
pop_ids <- readLines(POP_FILE)
pop_ids <- unique(pop_ids[pop_ids %in% names(gene2GO_all)])
if (length(pop_ids) == 0) stop("No overlap between population.ids and gene2go.assoc")

has_file <- function(p) !is.null(p) && nzchar(p) && file.exists(p)

cons_ids <- if (has_file(CONS_FILE)) readLines(CONS_FILE) else character()
cons_ids <- unique(cons_ids[cons_ids %in% pop_ids])

noncons_ids <- if (has_file(NONCONS_FILE)) readLines(NONCONS_FILE) else character()
noncons_ids <- unique(noncons_ids[noncons_ids %in% pop_ids])

if (length(cons_ids) == 0 && length(noncons_ids) == 0) {
  stop("Both study sets are empty after intersecting with population.")
}

slim_set <- if (has_file(SLIM_FILE)) unique(readLines(SLIM_FILE)) else character()

# ======================
# === topGO helpers ====
# ======================
make_topgo <- function(ontology = c("BP","MF","CC"), study_ids, universe_ids, gene2GO, use_slim = FALSE, slim_terms = character()) {
  ontology <- match.arg(ontology)

  # Build 'geneList': named numeric vector over universe; 1 = study gene
  geneList <- factor(as.integer(universe_ids %in% study_ids))
  names(geneList) <- universe_ids

  # topGO requires mapping for the universe only
  geneID2GO <- gene2GO[universe_ids]
  geneID2GO <- geneID2GO[!sapply(geneID2GO, is.null)]

  if (length(geneID2GO) == 0) return(NULL)

  obj <- new("topGOdata",
             ontology = ontology,
             allGenes = geneList,
             annot = annFUN.gene2GO,
             gene2GO = geneID2GO,
             nodeSize = 1)  # allow small terms; we’ll filter later

  if (use_slim && length(slim_terms) > 0) {
    # We can't “restrict” the graph in topGO directly, but we can post-filter results to slim.
    slot(obj, "description") <- paste(slot(obj, "description"), "(Slim view)", sep = " ")
  }
  obj
}

run_topgo <- function(tgd) {
  if (is.null(tgd)) return(NULL)
  # 3 scoring strategies
  res_classic  <- runTest(tgd, algorithm = "classic",  statistic = "fisher")
  res_elim     <- runTest(tgd, algorithm = "elim",     statistic = "fisher")
  res_weight01 <- runTest(tgd, algorithm = "weight01", statistic = "fisher")
  
  allGO <- usedGO(tgd)
  tab <- GenTable(tgd,
                  classic  = res_classic,
                  elim     = res_elim,
                  weight01 = res_weight01,
                  orderBy  = "weight01",
                  ranksOf  = "weight01",
                  topNodes = length(allGO))
  
  # numeric counts from the universe/study
  # tgd@allGenes is a factor with levels "0","1"
  all_vec <- as.integer(as.character(tgd@allGenes))
  study_n <- sum(all_vec == 1, na.rm = TRUE)
  pop_n   <- length(all_vec)
  
  # clean and attach counts
  num_cols <- c("classic","elim","weight01")
  for (c in num_cols) tab[[c]] <- suppressWarnings(as.numeric(tab[[c]]))
  colnames(tab)[colnames(tab)=="Significant"] <- "study_count"
  colnames(tab)[colnames(tab)=="Annotated"]   <- "pop_count"
  tab$study_n <- study_n
  tab$pop_n   <- pop_n
  tab$NS      <- switch(tgd@ontology, BP="biological_process", MF="molecular_function", CC="cellular_component")
  colnames(tab)[colnames(tab)=="GO.ID"] <- "GO"
  colnames(tab)[colnames(tab)=="Term"]  <- "name"
  
  # long format with a BH-FDR per algorithm
  long <- tab %>%
    dplyr::select(GO, name, NS, study_count, study_n, pop_count, pop_n, classic, elim, weight01) %>%
    tidyr::pivot_longer(cols = c(classic, elim, weight01), names_to = "method", values_to = "p_raw") %>%
    dplyr::mutate(p_fdr = p.adjust(p_raw, method = "BH"))
  
  long
}


restrict_to_slim <- function(df, slim_terms) {
  if (length(slim_terms) == 0) return(df)
  df %>% filter(GO %in% slim_terms)
}

tidy_and_filter <- function(df, min_study = 1, min_pop = 3) {
  df %>%
    filter(!is.na(p_raw), study_count >= min_study, pop_count >= min_pop)
}

# ======================
# ==== RUN ANALYSIS ====
# ======================
sets <- list(
  Consensus     = cons_ids,
  `Non-consensus` = noncons_ids
)
sets <- sets[sapply(sets, length) > 0]

all_outputs <- list()
for (set_name in names(sets)) {
  study <- sets[[set_name]]
  message(sprintf("[%s] study=%d, population=%d", set_name, length(study), length(pop_ids)))

  for (onto in c("BP","MF","CC")) {
    tgd <- make_topgo(onto, study, pop_ids, gene2GO_all, use_slim = nzchar(SLIM_FILE), slim_terms = slim_set)
    if (is.null(tgd)) next
    res <- run_topgo(tgd)
    if (is.null(res) || nrow(res)==0) next

    res <- tidy_and_filter(res, min_study = 1, min_pop = 3)
    if (nzchar(SLIM_FILE)) res <- restrict_to_slim(res, slim_set)

    # save TSV
    out_tsv <- sprintf("%s_%s_%s.tsv", PREFIX, tolower(gsub("[^A-Za-z]+", "", set_name)), tolower(onto))
    write_tsv(res, out_tsv)
    message("Wrote: ", out_tsv)

    # store for plotting
    res$Set <- set_name
    res$Onto <- onto
    all_outputs[[paste(set_name, onto, sep="::")]] <- res
  }
}

if (length(all_outputs) == 0) {
  stop("No results to plot. Check that your study & population overlap with gene2go, and that GO coverage exists.")
}

all_df <- bind_rows(all_outputs)

# ======================
# ======= PLOTS ========
# ======================
plot_one <- function(d, set_name, onto, top_n = TOP_N) {
  # top by p_raw within each method (so you can compare classic/elim/weight01)
  top <- d %>%
    filter(Set == set_name, Onto == onto) %>%
    group_by(method) %>%
    arrange(p_raw, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()

  if (nrow(top) == 0) return(NULL)

  top <- top %>%
    mutate(
      neglogp = -log10(pmax(p_raw, .Machine$double.xmin)),
      gene_ratio = study_count / pmax(study_n, 1),
      method = factor(method, levels = c("weight01","elim","classic")),
      name_fac = factor(name, levels = rev(unique(name[order(p_raw, decreasing = FALSE)])))
    )

  p <- ggplot(top, aes(x = gene_ratio, y = name_fac)) +
    geom_point(aes(size = study_count, fill = neglogp, shape = method),
               color = "grey15", alpha = 0.95) +
    scale_size_continuous(name = "Study count") +
    scale_fill_viridis(name = expression(-log[10](p)), option = "C") +
    scale_shape_manual(values = c(21, 22, 24), name = "Algorithm") +
    labs(x = "Gene ratio (study_count / study_n)", y = "GO term",
         title = paste0(set_name, " — ", onto, " (top ", top_n, " per method)")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.title.y = element_text(margin = margin(r = 10)))
  p
}

for (set_name in unique(all_df$Set)) {
  for (onto in c("BP","MF","CC")) {
    g <- plot_one(all_df, set_name, onto, TOP_N)
    if (!is.null(g)) {
      png(sprintf("%s_%s_%s_bubbles.png", PREFIX, tolower(gsub("[^A-Za-z]+","",set_name)), tolower(onto)),
          width = 1200, height = max(800, 40*TOP_N), res = 120)
      print(g)
      dev.off()

      ggsave(sprintf("%s_%s_%s_bubbles.pdf", PREFIX, tolower(gsub("[^A-Za-z]+","",set_name)), tolower(onto)),
             g, width = 10, height = max(6, 0.4*TOP_N))
    }
  }
}

# Also write combined table with an FDR flag for convenience
all_df <- all_df %>%
  mutate(sig_fdr = p_fdr < ALPHA_FDR)
write_tsv(all_df, sprintf("%s_all_results.tsv", PREFIX))

message("\nDone. Outputs include:")
message("  - ", paste0(PREFIX, "_consensus_bp/mf/cc.tsv (if consensus provided)"))
message("  - ", paste0(PREFIX, "_nonconsensus_bp/mf/cc.tsv (if non-consensus provided)"))
message("  - ", paste0(PREFIX, "_*_bubbles.png/.pdf bubble plots (per set × namespace)"))
message("  - ", paste0(PREFIX, "_all_results.tsv (combined)"))
if (nzchar(SLIM_FILE)) message("  - Results restricted to terms in GO-Slim: ", SLIM_FILE)

