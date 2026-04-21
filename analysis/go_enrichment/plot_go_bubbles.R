#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(viridis)
  library(tidyr)
})

# =========================
# ======= CONFIG ==========
# =========================
# Put your file paths here. Leave non-consensus as "" if you only want consensus.
files_consensus <- c(
  "enrich_consensus_bi.tsv",  # BP
  "enrich_consensus_mo.tsv",  # MF
  "enrich_consensus_ce.tsv"   # CC
)

files_nonconsensus <- c(
  "enrich_nonconsensus_bi.tsv",  # BP
  "enrich_nonconsensus_mo.tsv",  # MF
  "enrich_nonconsensus_ce.tsv"   # CC
)

# Plot parameters
PREFIX       <- "go_debug"     # output prefix
TOP_N        <- 15             # top terms per namespace&set
RANK_COL     <- "p_raw"        # "p_raw" or "p_fdr" to rank/score terms
X_MODE       <- "ratio"        # "ratio" or "neglogp"
ADD_JITTER   <- FALSE          # set TRUE to add tiny horizontal jitter to separate identical x-values

# =========================
# ======= HELPERS =========
# =========================
need_cols <- c("GO","name","NS","p_raw","p_fdr","study_count","study_n","pop_count","pop_n")

read_ns <- function(fpath) {
  if (!file.exists(fpath)) stop(paste("File not found:", fpath))
  df <- suppressMessages(read_tsv(fpath, show_col_types = FALSE))
  miss <- setdiff(need_cols, names(df))
  if (length(miss) > 0) stop(paste("Missing columns in", fpath, ":", paste(miss, collapse=", ")))

  # coerce p columns to numeric in case they got read as character
  df <- df %>%
    mutate(
      p_raw = suppressWarnings(as.numeric(p_raw)),
      p_fdr = suppressWarnings(as.numeric(p_fdr))
    )

  # normalize NS labels
  df <- df %>%
    mutate(
      NS = recode(NS,
                  biological_process   = "BP",
                  molecular_function   = "MF",
                  cellular_component   = "CC",
                  .default = NS)
    )
  df$source_file <- basename(fpath)
  df
}

infer_set_label <- function(fps, fallback) {
  s <- tolower(paste(fps, collapse = " "))
  if (str_detect(s, "noncon")) return("Non-consensus")
  if (str_detect(s, "consensus")) return("Consensus")
  return(fallback)
}

load_set <- function(files3, label_fallback) {
  if (all(files3 != "") && all(file.exists(files3))) {
    d <- bind_rows(lapply(files3, read_ns))
    d$Set <- infer_set_label(files3, label_fallback)
    return(d)
  }
  return(NULL)
}

summarize_ns <- function(df, label) {
  if (nrow(df) == 0) return(invisible())
  cat("\n==== Diagnostics:", label, "====\n")
  # Overall
  cat("Rows:", nrow(df), "\n")
  # Per namespace summary
  s <- df %>%
    group_by(NS) %>%
    summarise(
      n_terms = n(),
      uniq_p_raw   = n_distinct(p_raw, na.rm = TRUE),
      p_raw_min    = min(p_raw, na.rm = TRUE),
      p_raw_max    = max(p_raw, na.rm = TRUE),
      uniq_p_fdr   = n_distinct(p_fdr, na.rm = TRUE),
      uniq_ratio   = n_distinct(study_count / pmax(study_n, 1)),
      uniq_study_n = n_distinct(study_n),
      uniq_study_c = n_distinct(study_count),
      .groups = "drop"
    )
  print(s, n = Inf)

  # Warn on low variation
  s %>%
    mutate(flag_one_color = uniq_p_raw <= 1,
           flag_one_x     = uniq_ratio <= 1) %>%
    rowwise() %>%
    mutate(msg = paste0(
      if (flag_one_color) "[WARN one color (neglogp constant)] " else "",
      if (flag_one_x)     "[WARN x collapses (ratio constant)] " else ""
    )) %>%
    ungroup() %>%
    filter(msg != "") %>%
    { if (nrow(.)>0) print(.) }

  # Show a tiny sample for each NS (first 5 rows by p_raw)
  sample_tbl <- df %>%
    arrange(p_raw) %>%
    group_by(NS) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    select(NS, name, p_raw, p_fdr, study_count, study_n, pop_count, pop_n, source_file)
  cat("\n-- Small sample (first 5 lowest p_raw per NS):\n")
  print(sample_tbl, n = Inf)
}

# =========================
# ======== LOAD =========== 
# =========================
d1 <- load_set(files_consensus, "Consensus")
d2 <- load_set(files_nonconsensus, "Non-consensus")

if (is.null(d1) && is.null(d2)) {
  stop("No valid input files found. Check the CONFIG paths at the top of the script.")
}

if (!is.null(d1) && is.null(d2)) {
  df_all <- d1
} else if (is.null(d1) && !is.null(d2)) {
  df_all <- d2
} else {
  df_all <- bind_rows(d1, d2)
}

# Compute derived fields
# choose which p to use for ranking/neglogp
RANK_COL <- if (RANK_COL %in% c("p_raw","p_fdr")) RANK_COL else "p_raw"

df_all <- df_all %>%
  mutate(
    rank_p   = ifelse(!is.na(.data[[RANK_COL]]), .data[[RANK_COL]], p_raw),
    neglogp  = -log10(pmax(rank_p, .Machine$double.xmin)),
    gene_ratio = study_count / pmax(study_n, 1)
  )

# Diagnostics before trimming
if (!is.null(d1)) summarize_ns(d1, "Consensus (raw input)")
if (!is.null(d2)) summarize_ns(d2, "Non-consensus (raw input)")

# Select top-N per (Set, NS) by RANK_COL
df_top <- df_all %>%
  group_by(Set, NS) %>%
  arrange(rank_p, .by_group = TRUE) %>%
  slice_head(n = TOP_N) %>%
  ungroup()

# Order y axis by rank within each facet
df_top <- df_top %>%
  group_by(Set, NS) %>%
  mutate(name_fac = factor(name, levels = rev(unique(name[order(rank_p, decreasing = FALSE)])))) %>%
  ungroup()

# More diagnostics on the trimmed set
cat("\n==== Diagnostics: Trimmed top-N per (Set, NS) ====\n")
diag2 <- df_top %>%
  group_by(Set, NS) %>%
  summarise(
    n_terms     = n(),
    uniq_neglogp = n_distinct(neglogp),
    uniq_ratio   = n_distinct(gene_ratio),
    uniq_study_c = n_distinct(study_count),
    uniq_study_n = n_distinct(study_n),
    .groups = "drop"
  )
print(diag2, n = Inf)

# =========================
# ========= PLOT ==========
# =========================
xvar <- if (tolower(X_MODE) == "neglogp") "neglogp" else "gene_ratio"
xlab <- if (xvar == "neglogp") expression(-log[10](p)) else "Gene ratio (study_count / study_n)"

base <- ggplot(df_top, aes(x = .data[[xvar]], y = name_fac))

if (ADD_JITTER) {
  base <- base + geom_point(
    aes(size = study_count, fill = neglogp, color = Set),
    shape = 21, alpha = 0.95, position = position_jitter(width = ifelse(xvar=="gene_ratio", 0.002, 0.02), height = 0)
  )
} else {
  base <- base + geom_point(
    aes(size = study_count, fill = neglogp, color = Set),
    shape = 21, alpha = 0.95
  )
}

plt <- base +
  scale_size_continuous(name = "Study count") +
  scale_fill_viridis_c(name = expression(-log[10](p)), option = "C") +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(NS ~ Set, scales = "free_y", switch = "y") +
  labs(x = xlab, y = "GO term",
       title = paste0(PREFIX, " — top ", TOP_N, " by ", RANK_COL, " | x = ", xvar)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.placement = "outside",
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.key.height = unit(0.6, "cm")
  )

if (xvar == "neglogp") {
  plt <- plt + geom_vline(xintercept = -log10(0.05), linetype = "dashed")
}

# =========================
# ======== SAVE ===========
# =========================
png_out <- paste0(PREFIX, "_bubble_new.png")
pdf_out <- paste0(PREFIX, "_bubble_new.pdf")
height_rows <- max(6, 0.35 * (nrow(df_top) / length(unique(df_top$Set))))  # rough auto height
ggsave(png_out, plt, width = 10, height = height_rows, dpi = 300)
ggsave(pdf_out, plt, width = 10, height = height_rows)

cat("\nWrote plots to:\n  ", png_out, "\n  ", pdf_out, "\n")

# =========================
# ====== EXTRA DEBUG ======
# =========================
# Write a small debug CSV per facet with the key numeric columns so you can inspect values causing overlaps.
debug_csv <- paste0(PREFIX, "_facet_debug.csv")
df_top %>%
  arrange(Set, NS, rank_p, name) %>%
  select(Set, NS, name, GO, p_raw, p_fdr, rank_p, neglogp, gene_ratio, study_count, study_n, pop_count, pop_n, source_file) %>%
  write_csv(debug_csv)

cat("Wrote debug table:\n  ", debug_csv, "\n\n")
cat("TIP: If you still see one color only,\n",
    " - check the diagnostics above: if uniq_neglogp == 1 for a facet → all points have same p (same color).\n",
    " - try X_MODE = 'neglogp' to spread points and verify p variation visually.\n",
    " - consider re-running GOEA with propagate_counts=FALSE to reduce ancestor duplication.\n")

