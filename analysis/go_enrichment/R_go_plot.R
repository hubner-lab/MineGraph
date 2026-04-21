#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste0(
    "Usage:\n",
    "  Rscript plot_go_bubbles.R <BP.tsv> <MF.tsv> <CC.tsv> [<BP2.tsv> <MF2.tsv> <CC2.tsv>] ",
    "--prefix <name> --top <N> --rank <p_raw|p_fdr> --x <ratio|neglogp>\n",
    "\nExamples:\n",
    "  Rscript plot_go_bubbles.R enrich_consensus_bi.tsv enrich_consensus_mo.tsv enrich_consensus_ce.tsv ",
    "--prefix go_consensus --top 15 --rank p_raw --x ratio\n",
    "  Rscript plot_go_bubbles.R enrich_consensus_bi.tsv enrich_consensus_mo.tsv enrich_consensus_ce.tsv ",
    "enrich_nonconsensus_bi.tsv enrich_nonconsensus_mo.tsv enrich_nonconsensus_ce.tsv ",
    "--prefix go_cmp --top 15 --rank p_raw --x neglogp\n"
  ))
}

# ---- parse positional inputs ----
# We expect 3 or 6 TSVs in total. Anything extra are flags.
pos <- args
# find where flags start
flag_idx <- which(grepl("^--", args))
if (length(flag_idx) > 0) {
  pos <- args[seq_len(min(flag_idx) - 1)]
}
if (!(length(pos) %in% c(3,6))) stop("Provide either 3 (one set) or 6 (two sets) TSV files.")

files <- pos

# ---- defaults & flags ----
prefix  <- "go_out"
top_n   <- 15L
rankcol <- "p_raw"   # or p_fdr
xmode   <- "ratio"   # or neglogp

if (length(flag_idx) > 0) {
  flags <- args[min(flag_idx):length(args)]
  for (i in seq(1, length(flags), by = 2)) {
    k <- flags[i]; v <- ifelse(i + 1 <= length(flags), flags[i + 1], NA)
    if (k %in% c("--prefix","-p")) prefix  <- v
    if (k %in% c("--top","-t"))    top_n   <- as.integer(v)
    if (k %in% c("--rank","-r"))   rankcol <- v
    if (k %in% c("--x","-x"))      xmode   <- v
  }
}

# ---- helpers ----
read_ns <- function(fpath) {
  if (!file.exists(fpath)) stop(paste("File not found:", fpath))
  df <- suppressMessages(read_tsv(fpath, show_col_types = FALSE))
  need <- c("GO","name","NS","p_raw","p_fdr","study_count","study_n","pop_count","pop_n")
  miss <- setdiff(need, names(df))
  if (length(miss) > 0) stop(paste("Missing columns in", fpath, ":", paste(miss, collapse=", ")))
  df$source_file <- basename(fpath)
  df
}

infer_set_label <- function(fn) {
  s <- tolower(fn)
  if (str_detect(s, "noncon")) return("Non-consensus")
  if (str_detect(s, "consensus")) return("Consensus")
  # fallback from trailing group (1 or 2)
  return(NA_character_)
}

# ---- load data ----
if (length(files) == 3) {
  d1 <- bind_rows(lapply(files, read_ns))
  d1$Set <- "Consensus"  # or just "Set1"
  df <- d1
} else { # 6 files: combine two sets
  set1 <- files[1:3]; set2 <- files[4:6]
  d1 <- bind_rows(lapply(set1, read_ns)); d1$Set <- infer_set_label(set1[1])
  if (is.na(d1$Set[1])) d1$Set <- "Set1"
  d2 <- bind_rows(lapply(set2, read_ns)); d2$Set <- infer_set_label(set2[1])
  if (is.na(d2$Set[1])) d2$Set <- "Set2"
  df <- bind_rows(d1, d2)
}

# ---- prepare ----
# normalize namespace short labels
df <- df %>%
  mutate(
    NS = recode(NS,
                biological_process   = "BP",
                molecular_function   = "MF",
                cellular_component   = "CC",
                .default = NS),
    neglogp = -log10(pmax(ifelse(rankcol %in% names(.), .data[[rankcol]], p_raw), .Machine$double.xmin)),
    gene_ratio = study_count / pmax(study_n, 1)
  )

# guard on rank column
if (!rankcol %in% names(df)) {
  warning("Requested rank column not found; using p_raw")
  rankcol <- "p_raw"
}

# select top-N per (NS, Set)
df_top <- df %>%
  group_by(Set, NS) %>%
  arrange(.data[[rankcol]], .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()

# order terms within facet
df_top <- df_top %>%
  group_by(Set, NS) %>%
  mutate(name_fac = factor(name, levels = rev(unique(name[order(.data[[rankcol]], decreasing = FALSE)])))) %>%
  ungroup()

# choose x variable
if (tolower(xmode) == "neglogp") {
  xvar <- "neglogp";  xlab <- expression(-log[10](p))
} else {
  xvar <- "gene_ratio"; xlab <- "Gene ratio (study_count / study_n)"
}

# ---- plot ----
plt <- ggplot(df_top, aes(x = .data[[xvar]], y = name_fac)) +
  geom_point(aes(size = study_count, fill = neglogp),
             shape = 21, color = "grey25", alpha = 0.95) +
  scale_size_continuous(name = "Study count") +
  scale_fill_viridis_c(name = expression(-log[10](p)), option = "C") +
  facet_grid(NS ~ Set, scales = "free_y", switch = "y") +
  labs(
    x = xlab,
    y = "GO term",
    title = paste0(prefix, " — top ", top_n, " by ", rankcol, " | x = ", xvar)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.placement = "outside",
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.key.height = unit(0.6, "cm")
  )

# dashed reference line for FDR=0.05 if plotting neglogp
if (tolower(xmode) == "neglogp") {
  plt <- plt + geom_vline(xintercept = -log10(0.05), linetype = "dashed")
}

# ---- save ----
png_out <- paste0(prefix, "_bubble.png")
pdf_out <- paste0(prefix, "_bubble.pdf")
height_rows <- max(6, 0.35 * (nrow(df_top) / length(unique(df_top$Set))))
ggsave(png_out, plt, width = 10, height = height_rows, dpi = 300)
ggsave(pdf_out, plt, width = 10, height = height_rows)
message("Wrote: ", png_out, " and ", pdf_out)

