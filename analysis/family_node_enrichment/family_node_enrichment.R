#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Family-level node analysis from odgi paths -H matrix:
#   - nodes "unique/enriched" for each family
#   - nodes shared across families
#   - per-family output folders + final summary TSV
#
# Inputs:
#   odgi_matrix.tsv  (from: odgi paths -H)
#   Supplementary_Table_S3_Full_Poales_Data.xlsx  (sheet: "Full List")
#
# Outputs:
#   results_families/
#     <Family>/
#       family_nodes_enriched.tsv
#       family_nodes_unique.tsv
#       family_nodes_core_within_family.tsv
#       family_summary.tsv
#     ALL_FAMILIES_node_summary.tsv
#     ALL_FAMILIES_family_summary.tsv
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(stringr)
  library(fs)
})

# ================== USER PARAMETERS =========================

matrix_file <- "odgi_matrix.tsv"
xlsx_s3     <- "Supplementary_Table_S3_Full_Poales_Data.xlsx"
sheet_name  <- "Full List"

out_root    <- "results_families"

# Thresholds (tune if needed)
# "present in family" means fraction of species in family having node=1 is >= in_family_min
in_family_min        <- 0.80

# "absent outside" means fraction outside family having node=1 is <= out_family_max
out_family_max       <- 0.10

# "unique to family" means node never appears outside family (outside_fraction == 0)
unique_outside_exact <- TRUE

# Optional filters for stability
min_family_size      <- 1     # skip tiny families (still reported in summary)
min_total_presence   <- 2     # node must appear in at least this many species overall

# ============================================================
# Helper: make tip labels compatible with your previous mapping
# (same logic as your PAV heatmap script)
# ============================================================

make_raxml_tip <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\[|\\]", "")
  x <- str_replace_all(x, "-", "")
  x <- str_replace_all(x, "\\(|\\)", "")
  x <- str_replace_all(x, "[^A-Za-z0-9]", "_")
  substr(x, 1, 30)
}

# ============================================================
# 0) Create output root
# ============================================================
dir_create(out_root)

# ============================================================
# 1) Read metadata (Species -> Family)
# ============================================================

meta_raw <- read_excel(xlsx_s3, sheet = sheet_name) %>%
  mutate(across(everything(), as.character))

# drop blank column names if exist
if (any(names(meta_raw) == "")) {
  meta_raw <- meta_raw %>% select(-which(names(meta_raw) == ""))
}

if (!all(c("Species", "Family") %in% names(meta_raw))) {
  stop("Excel sheet must contain columns named 'Species' and 'Family'.")
}

meta_family <- meta_raw %>%
  transmute(
    Species   = Species,
    tip_label = make_raxml_tip(Species),
    Family    = if_else(is.na(Family) | Family == "", "UnassignedFamily", Family)
  ) %>%
  distinct(tip_label, Family, .keep_all = TRUE)

# ============================================================
# 2) Read odgi_matrix.tsv and build binary node matrix
# ============================================================

dt <- fread(matrix_file)

if (!("path.name" %in% names(dt))) {
  stop("odgi_matrix.tsv must include a 'path.name' column.")
}

node_cols <- grep("^node\\.", names(dt), value = TRUE)
if (length(node_cols) == 0) stop("No node.* columns found in odgi_matrix.tsv")

# Clean path names (remove #1#1 etc.) to match RAxML tip style
dt[, path.name := sub("#\\d+#\\d+$", "", path.name)]

# Convert node columns to numeric matrix and binarize
mat <- as.matrix(dt[, ..node_cols])
suppressWarnings(storage.mode(mat) <- "numeric")
mat <- (mat > 0) * 1L
rownames(mat) <- dt$path.name

# Filter nodes that are never present or too rare globally (optional)
total_presence <- colSums(mat)
keep_cols <- total_presence >= min_total_presence
mat <- mat[, keep_cols, drop = FALSE]
node_cols_kept <- colnames(mat)

cat("Total species:", nrow(mat), "\n")
cat("Node columns kept:", ncol(mat), " (after min_total_presence =", min_total_presence, ")\n")

# ============================================================
# 3) Map each path to Family
# ============================================================

map_df <- tibble(path.name = rownames(mat)) %>%
  mutate(tip_label = make_raxml_tip(path.name)) %>%
  left_join(meta_family %>% select(tip_label, Family), by = "tip_label") %>%
  mutate(Family = if_else(is.na(Family), "UnassignedFamily", Family))

# Sanity check mapping
n_mapped <- sum(map_df$Family != "UnassignedFamily")
cat("Mapped to known Family:", n_mapped, "/", nrow(map_df), "\n")

# ============================================================
# 4) Core computations: per-family node frequencies
# ============================================================

families <- sort(unique(map_df$Family))
family_sizes <- map_df %>%
  count(Family, name = "n_in_family") %>%
  arrange(desc(n_in_family))

# Precompute global outside-family convenience
all_rows <- rownames(mat)

# We'll build:
# - per-family summary (counts)
# - per-family node tables
all_family_summaries <- list()
all_family_node_summary_rows <- list()

for (fam in families) {
  
  fam_dir <- path(out_root, fam)
  dir_create(fam_dir)
  
  fam_rows <- map_df %>% filter(Family == fam) %>% pull(path.name)
  out_rows <- setdiff(all_rows, fam_rows)
  
  n_in  <- length(fam_rows)
  n_out <- length(out_rows)
  
  fam_summary <- tibble(
    Family = fam,
    n_in_family = n_in,
    n_outside   = n_out
  )
  
  # Skip small families for “enriched/unique” calls (still record summary)
  if (n_in < min_family_size) {
    fam_summary <- fam_summary %>%
      mutate(
        status = paste0("SKIPPED_small_family(<", min_family_size, ")"),
        enriched_nodes = 0,
        unique_nodes   = 0,
        core_within_family_nodes = 0
      )
    fwrite(fam_summary, file = path(fam_dir, "family_summary.tsv"), sep = "\t")
    all_family_summaries[[fam]] <- fam_summary
    next
  }
  
  # Fractions in family and outside for each node
  in_frac  <- colMeans(mat[fam_rows, , drop = FALSE])
  out_frac <- if (n_out > 0) colMeans(mat[out_rows, , drop = FALSE]) else rep(0, length(in_frac))
  names(out_frac) <- names(in_frac)
  
  # Counts (more interpretable)
  in_count  <- colSums(mat[fam_rows, , drop = FALSE])
  out_count <- if (n_out > 0) colSums(mat[out_rows, , drop = FALSE]) else rep(0L, length(in_frac))
  names(out_count) <- names(in_frac)
  
  # Definitions
  # Core-within-family: present in (almost) all family members
  core_within <- (in_frac >= 0.95)
  
  # Enriched in family + mostly absent outside
  enriched <- (in_frac >= in_family_min) & (out_frac <= out_family_max)
  
  # Unique to family: enriched AND absent outside (exactly)
  if (unique_outside_exact) {
    unique <- enriched & (out_count == 0L)
  } else {
    unique <- enriched & (out_frac <= out_family_max)
  }
  
  # Build per-node table (for this family)
  node_tbl <- tibble(
    node = names(in_frac),
    in_count  = as.integer(in_count),
    in_frac   = as.numeric(in_frac),
    out_count = as.integer(out_count),
    out_frac  = as.numeric(out_frac),
    is_core_within_family = core_within,
    is_enriched_family    = enriched,
    is_unique_family      = unique
  ) %>%
    arrange(desc(is_unique_family), desc(is_enriched_family), desc(in_frac), out_frac)
  
  # Write family node outputs
  fwrite(node_tbl %>% filter(is_enriched_family),
         file = path(fam_dir, "family_nodes_enriched.tsv"), sep = "\t")
  
  fwrite(node_tbl %>% filter(is_unique_family),
         file = path(fam_dir, "family_nodes_unique.tsv"), sep = "\t")
  
  fwrite(node_tbl %>% filter(is_core_within_family),
         file = path(fam_dir, "family_nodes_core_within_family.tsv"), sep = "\t")
  
  # Family summary counts
  fam_summary <- fam_summary %>%
    mutate(
      status = "OK",
      enriched_nodes = sum(enriched),
      unique_nodes   = sum(unique),
      core_within_family_nodes = sum(core_within)
    )
  
  fwrite(fam_summary, file = path(fam_dir, "family_summary.tsv"), sep = "\t")
  
  all_family_summaries[[fam]] <- fam_summary
  
  # Add to global node summary: one row per (family, node) for enriched/unique
  node_summary_small <- node_tbl %>%
    filter(is_enriched_family | is_unique_family) %>%
    mutate(Family = fam) %>%
    select(Family, node, in_count, in_frac, out_count, out_frac,
           is_enriched_family, is_unique_family)
  
  all_family_node_summary_rows[[fam]] <- node_summary_small
  
  cat("Family:", fam, "n=", n_in,
      " enriched=", sum(enriched),
      " unique=", sum(unique),
      " core_within=", sum(core_within), "\n")
}

# ============================================================
# 5) Final global outputs
# ============================================================

family_summary_all <- bind_rows(all_family_summaries) %>%
  arrange(desc(n_in_family))

fwrite(family_summary_all,
       file = path(out_root, "ALL_FAMILIES_family_summary.tsv"),
       sep = "\t")

node_summary_all <- bind_rows(all_family_node_summary_rows) %>%
  arrange(Family, desc(is_unique_family), desc(in_frac), out_frac)

fwrite(node_summary_all,
       file = path(out_root, "ALL_FAMILIES_node_summary.tsv"),
       sep = "\t")

cat("\nDONE.\n")
cat("Wrote per-family folders in:", out_root, "\n")
cat("Global summaries:\n  -", path(out_root, "ALL_FAMILIES_family_summary.tsv"), "\n")
cat("  -", path(out_root, "ALL_FAMILIES_node_summary.tsv"), "\n")
