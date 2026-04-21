#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Reference-free PAV-like heatmap from:
#   odgi paths -H  -> odgi_matrix.tsv
#
# Uses graph node sort order (node.1..node.n), bins nodes for plotting,
# and selects stratified representatives by Family (medoids).
#
# Output: PAV_from_odgi_matrix.png
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(vegan)
  library(readxl)
  library(stringr)
  library(ggnewscale)
})

# ================== USER PARAMETERS =========================

matrix_file       <- "odgi_matrix.tsv"                 # from: odgi paths -H
xlsx_s3           <- "Supplementary_Table_S3_Full_Poales_Data.xlsx"
output_png        <- "PAV_from_odgi_matrix2212.png"

total_reps        <- 25
per_family_target <- 10

# Bin size in "node columns" to make the heatmap tractable.
# If your matrix has many node.* columns (e.g., 50k), use 200–1000.
node_bin_size     <- 500

# Label placement/margins (same approach you liked)
label_hjust       <- 1.15
left_margin_pt    <- 260

# ============================================================
# FAMILY MAPPING (same as your t-SNE / other scripts)
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

# Excel may have an empty first column; ignore it safely
meta_raw <- read_excel(xlsx_s3, sheet = "Full List") %>%
  mutate(across(everything(), as.character))

# If there is a blank column name (common from Excel), drop it
if (any(names(meta_raw) == "")) {
  meta_raw <- meta_raw %>% select(-which(names(meta_raw) == ""))
}

if (!all(c("Species", "Family") %in% names(meta_raw))) {
  stop("Excel sheet 'Full List' must contain columns named 'Species' and 'Family'.")
}

meta_family <- meta_raw %>%
  mutate(
    tip_label = make_raxml_tip(Species),
    Family    = if_else(is.na(Family) | Family == "", "UnassignedFamily", Family)
  ) %>%
  select(tip_label, Family) %>%
  distinct()

# ============================================================
# 1) READ ODGI -H MATRIX
# ============================================================

dt <- fread(matrix_file)

if (!("path.name" %in% names(dt))) {
  stop("odgi_matrix.tsv must include a 'path.name' column (first column).")
}

node_cols <- grep("^node\\.", names(dt), value = TRUE)
if (length(node_cols) == 0) stop("No node.* columns found in odgi_matrix.tsv")

# Clean path names (remove #1#1, #2#2, etc.)
dt[, path.name := sub("#\\d+#\\d+$", "", path.name)]

# Convert node columns -> numeric matrix, then binarize (FAST)
mat <- as.matrix(dt[, ..node_cols])

# Ensure numeric (fread usually gives integer, but be safe)
suppressWarnings(storage.mode(mat) <- "numeric")

# Binary presence (vectorized, VERY FAST)
mat <- (mat > 0) * 1L

# Attach row names
rownames(mat) <- dt$path.name

# ============================================================
# 2) DROP NODES WITH NO VARIANCE (FAST FOR BINARY)
#    SD>0  <=>  colSum in (0, nrow)
# ============================================================

col_sums <- colSums(mat)
var_keep <- (col_sums > 0) & (col_sums < nrow(mat))
mat <- mat[, var_keep, drop = FALSE]

if (ncol(mat) < 2) stop("Too few variable node columns after filtering.")

# ============================================================
# 3) MAP ALL PATHS -> FAMILY
# ============================================================

all_family <- tibble(group = rownames(mat)) %>%
  mutate(tip_label = make_raxml_tip(group)) %>%
  left_join(meta_family, by = "tip_label") %>%
  mutate(Family = if_else(is.na(Family), "UnassignedFamily", Family))

# ============================================================
# 4) DISTANCE, CLUSTERING
# ============================================================

dist_mat <- vegdist(mat, method = "jaccard")
hc <- hclust(dist_mat, method = "ward.D2")
D <- as.matrix(dist_mat)

# Helper: choose k medoids within a set using the distance matrix
pick_medoids <- function(members, k, D) {
  members <- unique(members)
  n <- length(members)
  if (n == 0) return(character(0))
  if (n <= k) return(members)
  
  subD <- D[members, members, drop = FALSE]
  
  if (requireNamespace("cluster", quietly = TRUE)) {
    pam_fit <- cluster::pam(as.dist(subD), k = k, diss = TRUE)
    return(rownames(subD)[pam_fit$id.med])
  } else {
    # Greedy fallback
    chosen <- character(0)
    remaining <- members
    for (i in seq_len(k)) {
      subR <- subD[remaining, remaining, drop = FALSE]
      avgD <- rowMeans(subR)
      med  <- names(which.min(avgD))
      chosen <- c(chosen, med)
      remaining <- setdiff(remaining, med)
      if (length(remaining) == 0) break
    }
    return(chosen)
  }
}

# ============================================================
# 5) STRATIFIED REPRESENTATIVE SELECTION (per-family medoids)
# ============================================================

families <- all_family %>%
  count(Family, name = "n") %>%
  arrange(desc(n)) %>%
  pull(Family)

selected <- character(0)

# First pass: up to per_family_target per family
for (fam in families) {
  members <- all_family %>% filter(Family == fam) %>% pull(group)
  k <- min(per_family_target, length(members))
  selected <- c(selected, pick_medoids(members, k, D))
}
selected <- unique(selected)

# Fill to total_reps by taking extra medoids from families with remaining members
remaining_needed <- total_reps - length(selected)
if (remaining_needed > 0) {
  remaining_map <- all_family %>%
    group_by(Family) %>%
    summarise(remaining = list(setdiff(group, selected)), .groups = "drop")
  
  while (remaining_needed > 0) {
    remaining_sizes <- remaining_map %>% mutate(n = lengths(remaining))
    if (all(remaining_sizes$n == 0)) break
    
    fam_pick <- remaining_sizes %>% arrange(desc(n)) %>% slice(1) %>% pull(Family)
    pool <- remaining_map %>% filter(Family == fam_pick) %>% pull(remaining) %>% .[[1]]
    if (length(pool) == 0) break
    
    add_one <- pick_medoids(pool, 1, D)
    selected <- unique(c(selected, add_one))
    
    remaining_map <- remaining_map %>%
      mutate(remaining = if_else(Family == fam_pick,
                                 list(setdiff(pool, add_one)),
                                 remaining))
    
    remaining_needed <- total_reps - length(selected)
  }
}

cat("Selected representative paths:\n")
print(selected)

# Order representatives by dendrogram order
dend_order <- rownames(mat)[hc$order]
rep_order  <- intersect(dend_order, selected)

# ============================================================
# 6) MAP REPRESENTATIVES -> FAMILY + PALETTE
# ============================================================

rep_family <- tibble(group = rep_order) %>%
  mutate(tip_label = make_raxml_tip(group)) %>%
  left_join(meta_family, by = "tip_label") %>%
  mutate(Family = if_else(is.na(Family), "UnassignedFamily", Family))

fam_levels <- sort(unique(rep_family$Family))
n_fam <- length(fam_levels)

if (requireNamespace("Polychrome", quietly = TRUE) && n_fam > 3) {
  seed_cols <- grDevices::hcl.colors(min(8, n_fam), palette = "Dark 3")
  fam_pal <- Polychrome::createPalette(n_fam, seedcolors = seed_cols, M = 1000)
} else {
  fam_pal <- grDevices::hcl.colors(n_fam, palette = "Dark 3")
}
names(fam_pal) <- fam_levels

# ============================================================
# 7) BIN NODE COLUMNS + BUILD LONG DATA FOR PLOT
# ============================================================

# Subset matrix to representatives (rows)
mat_rep <- mat[rep_order, , drop = FALSE]

# Recover node indices from column names node.X
node_idx <- as.integer(sub("^node\\.", "", colnames(mat_rep)))

# Bin assignment per column
bin_id <- floor((node_idx - 1) / node_bin_size) + 1
n_bins <- max(bin_id, na.rm = TRUE)

# Compute per-path per-bin mean presence (0..1)
# Faster than pivot_longer for huge matrices
binned_mat <- sapply(seq_len(n_bins), function(b) {
  cols <- which(bin_id == b)
  if (length(cols) == 1) {
    mat_rep[, cols]
  } else {
    rowMeans(mat_rep[, cols, drop = FALSE])
  }
})

# Ensure matrix shape
binned_mat <- as.matrix(binned_mat)
colnames(binned_mat) <- as.character(seq_len(n_bins))
rownames(binned_mat) <- rep_order

binned_df <- as_tibble(binned_mat, rownames = "path.name") %>%
  pivot_longer(
    cols = -path.name,
    names_to = "bin",
    values_to = "pav"
  ) %>%
  mutate(
    bin = as.integer(bin),
    path.name = factor(path.name, levels = rev(rep_order))
  )

# Labels (colored, in left margin)
labels_df <- rep_family %>%
  mutate(group = factor(group, levels = rev(rep_order))) %>%
  transmute(
    group,
    Family,
    x     = 1,
    label = as.character(group)
  )

# ============================================================
# 8) PLOT
# ============================================================

p <- ggplot() +
  geom_tile(
    data = binned_df,
    aes(x = bin, y = path.name, fill = pav),
    color = "black",
    size  = 0.05
  ) +
  scale_fill_gradient(low = "white", high = "brown", name = "PAV") +
  
  ggnewscale::new_scale_color() +
  geom_text(
    data = labels_df,
    aes(x = x, y = group, label = label, color = Family),
    hjust     = label_hjust,
    fontface  = "bold",
    size      = 5
  ) +
  scale_color_manual(values = fam_pal, guide = "none") +
  
  coord_cartesian(
    xlim = c(1, max(binned_df$bin, na.rm = TRUE)),
    clip = "off"
  ) +
  
  labs(
    x = paste0("Graph node bins (", node_bin_size, " nodes/bin; sort order)"),
    y = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    plot.margin  = margin(t = 10, r = 10, b = 10, l = left_margin_pt)
  )

ggsave(output_png, p, width = 14, height = 8, dpi = 300)
cat("\nSaved heatmap to:", output_png, "\n")
