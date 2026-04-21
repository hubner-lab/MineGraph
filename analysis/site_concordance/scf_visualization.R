#!/usr/bin/env Rscript
# =============================================================================
# A4+A5: sCF Visualization — Bootstrap vs. Site Concordance Factor
# Reads IQ-TREE .cf.stat files for two reference trees:
#   1. pav_scf_minegraph_bs  (MineGraph_BS as reference)
#   2. pav_scf_pav_ward      (PAV_Ward as reference)
# Produces:
#   - Scatter plot BS vs sCF per internal branch (both references)
#   - Summary statistics table
#   - Histogram of sCF distribution per reference
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
})
if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages("patchwork", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(patchwork))

BASE  <- "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
TREES <- file.path(BASE, "phylogenetic_trees")
FIGS  <- file.path(BASE, "manuscript/figures")

# ---- Helper: read cf.stat (skip comment lines starting with #) ----
read_cfstat <- function(path) {
  lines <- readLines(path)
  data_lines <- lines[!grepl("^#", lines)]
  con <- textConnection(paste(data_lines, collapse = "\n"))
  df  <- read.table(con, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  close(con)
  df
}

# ---- 1. Load stat files ----
cat("Loading sCF stat files ...\n")
scf_mg  <- read_cfstat(file.path(TREES, "pav_scf_minegraph_bs.cf.stat"))
scf_pav <- read_cfstat(file.path(TREES, "pav_scf_pav_ward.cf.stat"))

cat(sprintf("  MineGraph_BS ref: %d branches\n", nrow(scf_mg)))
cat(sprintf("  PAV_Ward ref:     %d branches\n", nrow(scf_pav)))

# ---- 2. Filter to internal branches only (non-NA sCF) ----
# IQ-TREE assigns NA sCF to root and tips; Label = bootstrap (may be empty)
filter_internal <- function(df) {
  df <- df[!is.na(df$sCF), ]
  # Convert Label to numeric bootstrap (NA if empty / not numeric)
  df$BS <- suppressWarnings(as.numeric(df$Label))
  df
}

mg_int  <- filter_internal(scf_mg)
pav_int <- filter_internal(scf_pav)

cat(sprintf("  Internal branches (MineGraph_BS): %d\n", nrow(mg_int)))
cat(sprintf("  Internal branches (PAV_Ward):     %d\n", nrow(pav_int)))

# ---- 3. Summary statistics ----
summarise_scf <- function(df, label) {
  data.frame(
    reference     = label,
    n_branches    = nrow(df),
    mean_sCF      = round(mean(df$sCF), 1),
    median_sCF    = round(median(df$sCF), 1),
    sd_sCF        = round(sd(df$sCF), 1),
    pct_above_33  = round(100 * mean(df$sCF > 33.3), 1),  # above random (1/3)
    pct_above_50  = round(100 * mean(df$sCF > 50), 1),
    pct_above_75  = round(100 * mean(df$sCF > 75), 1),
    cor_BS_sCF    = if (sum(!is.na(df$BS)) > 2) round(cor(df$BS, df$sCF, use = "complete.obs"), 3) else NA_real_,
    stringsAsFactors = FALSE
  )
}

stats_mg  <- summarise_scf(mg_int,  "MineGraph_BS")
stats_pav <- summarise_scf(pav_int, "PAV_Ward")
stats_all <- rbind(stats_mg, stats_pav)

cat("\n--- sCF Summary Statistics ---\n")
print(stats_all)

out_csv <- file.path(TREES, "scf_summary_statistics.csv")
write.csv(stats_all, out_csv, row.names = FALSE)
cat("Saved:", out_csv, "\n")

# ---- 4. Annotate branches by sCF support tier ----
annotate_tier <- function(df) {
  df$tier <- cut(df$sCF,
    breaks = c(-Inf, 33.3, 50, 75, Inf),
    labels = c("<33% (below random)", "33–50%", "50–75%", ">75% (strong)"),
    right  = TRUE
  )
  df
}

mg_int  <- annotate_tier(mg_int)
pav_int <- annotate_tier(pav_int)

tier_cols <- c(
  "<33% (below random)" = "#D73027",
  "33–50%"              = "#FC8D59",
  "50–75%"              = "#91BFDB",
  ">75% (strong)"       = "#1A9850"
)

# ---- 5. Scatter plot: BS vs sCF ----
make_scatter <- function(df, cor_val) {
  ggplot(df, aes(x = BS, y = sCF, colour = tier)) +
    geom_hline(yintercept = 33.3, linetype = "dashed", colour = "grey50", linewidth = 0.7) +
    geom_point(alpha = 0.65, size = 1.8) +
    scale_colour_manual(values = tier_cols, name = "sCF tier") +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    annotate("text", x = 5, y = 35.5, label = "1/3 random", hjust = 0,
             colour = "grey40", size = 3) +
    annotate("text", x = 95, y = 5,
             label = sprintf("r = %.3f", cor_val),
             hjust = 1, vjust = 0, size = 3.5, fontface = "italic") +
    labs(
      x        = "Bootstrap support (%)",
      y        = "Site concordance factor sCF (%)"
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "right",
          legend.key.size = unit(0.45, "cm"),
          legend.text     = element_text(size = 8))
}

# BS scatter only for MineGraph_BS (which embeds bootstrap values)
# PAV_Ward is a Ward-linkage tree with no bootstrap → no BS axis
mg_bs <- mg_int[!is.na(mg_int$BS), ]
p1 <- make_scatter(mg_bs, stats_mg$cor_BS_sCF)

# For PAV_Ward: sCF-only plot (histogram panel, skip scatter)
p2 <- ggplot(pav_int, aes(x = sCF, fill = tier)) +
  geom_histogram(binwidth = 5, boundary = 0, colour = "white") +
  geom_vline(xintercept = 33.3, linetype = "dashed", colour = "grey50", linewidth = 0.8) +
  scale_fill_manual(values = tier_cols, name = "sCF tier") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  annotate("text", x = 35, y = Inf, label = "1/3 random", hjust = 0, vjust = 1.5,
           size = 3, colour = "grey40") +
  labs(x = "sCF (%)", y = "Number of internal branches") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right", legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 8))

scatter_combined <- p1 / p2

out_scatter <- file.path(FIGS, "FigSX_scf_vs_bootstrap.png")
ggsave(out_scatter, scatter_combined, width = 9, height = 11, dpi = 300)
cat("Saved scatter:", out_scatter, "\n")

# ---- 6. sCF histogram (both references overlaid) ----
hist_df <- rbind(
  data.frame(sCF = mg_int$sCF,  reference = "MineGraph_BS"),
  data.frame(sCF = pav_int$sCF, reference = "PAV_Ward")
)

p_hist <- ggplot(hist_df, aes(x = sCF, fill = reference)) +
  geom_histogram(binwidth = 5, boundary = 0, alpha = 0.7,
                 position = "identity", colour = "white") +
  geom_vline(xintercept = 33.3, linetype = "dashed", colour = "grey40", linewidth = 0.9) +
  scale_fill_manual(values = c("MineGraph_BS" = "#0072B2", "PAV_Ward" = "#E69F00"),
                    name = "Reference tree") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  annotate("text", x = 35, y = Inf, label = "1/3 chance\n(random)", hjust = 0,
           vjust = 1.4, size = 3, colour = "grey40") +
  labs(
    x = "sCF (%)",
    y = "Number of internal branches"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.15, 0.85))

out_hist <- file.path(FIGS, "FigSX_scf_histogram.png")
ggsave(out_hist, p_hist, width = 8, height = 5, dpi = 300)
cat("Saved histogram:", out_hist, "\n")

# ---- 7. Branch-level table: ID, BS, sCF, tier ----
# Merge both references on branch ID for joint inspection
mg_tbl  <- mg_int[,  c("ID", "BS", "sCF", "sDF1", "sDF2", "Length", "tier")]
pav_tbl <- pav_int[, c("ID",       "sCF", "sDF1", "sDF2", "Length", "tier")]
names(mg_tbl)[3:7]  <- c("sCF_MG",  "sDF1_MG",  "sDF2_MG",  "Length_MG",  "tier_MG")
names(pav_tbl)[2:6] <- c("sCF_PAV", "sDF1_PAV", "sDF2_PAV", "Length_PAV", "tier_PAV")

branch_table <- merge(mg_tbl, pav_tbl, by = "ID", all = TRUE)
branch_table <- branch_table[order(branch_table$ID), ]

out_branch <- file.path(TREES, "scf_branch_table.csv")
write.csv(branch_table, out_branch, row.names = FALSE)
cat("Saved branch table:", out_branch, "\n")

cat("\n=== DONE A4+A5: sCF Visualization ===\n")
cat(sprintf("Key result (MineGraph_BS ref): mean sCF = %.1f%%, %.1f%% branches > 1/3, r(BS,sCF) = %.3f\n",
            stats_mg$mean_sCF, stats_mg$pct_above_33, stats_mg$cor_BS_sCF))
cat(sprintf("Key result (PAV_Ward ref):     mean sCF = %.1f%%, %.1f%% branches > 1/3, r(BS,sCF) = %.3f\n",
            stats_pav$mean_sCF, stats_pav$pct_above_33, stats_pav$cor_BS_sCF))
