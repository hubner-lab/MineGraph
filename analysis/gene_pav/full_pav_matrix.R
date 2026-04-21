#!/usr/bin/env Rscript

# ------------------------------------------------------------
# PAV heatmap from:
#   odgi paths -H  -> odgi_matrix.tsv
#
# Uses graph node sort order (node.1..node.n), bins nodes for plotting,
# and colors labels by Family.
#
# Output: PAV_from_odgi_matrix_ALL2.png
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

matrix_file       <- "odgi_matrix_192.tsv"                 # from: odgi paths -H
xlsx_s3           <- "Supplementary_Table_S3_Full_Poales_Data.xlsx"
output_png        <- "PAV_from_odgi_matrix_ALL_final.png"

# Bin size in "node columns" to make the heatmap tractable.
node_bin_size     <- 500

# Label placement/margins
label_hjust       <- 1.15
left_margin_pt    <- 260

# ---- Heatmap color scheme (absence clarity) ----
# PAV is a ratio 0..1 (after binning). 0 should look "empty" but still visible.
pav_low_color     <- "#F7FBFF"  # very light blue (instead of grey)
pav_high_color    <- "#08306B"  # dark navy
tile_border_color <- NA         # set to "black" if you want borders
tile_border_size  <- 0.00       # if borders enabled, try 0.05

# ============================================================
# FAMILY MAPPING
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

# Convert node columns -> numeric matrix, then binarize
mat <- as.matrix(dt[, ..node_cols])
suppressWarnings(storage.mode(mat) <- "numeric")
mat <- (mat > 0) * 1L
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
# 3.5) SUBSAMPLE POACEAE ONLY (KEEP OTHERS FULL)
# ============================================================

set.seed(1)                 # change if you want a different random 10
poaceae_keep_n <- 10        # <-- your request

poaceae_paths <- all_family %>%
  filter(Family == "Poaceae") %>%
  pull(group)

if (length(poaceae_paths) > poaceae_keep_n) {
  keep_poaceae <- sample(poaceae_paths, poaceae_keep_n)
} else {
  keep_poaceae <- poaceae_paths
}

keep_all <- c(
  keep_poaceae,
  all_family %>% filter(Family != "Poaceae") %>% pull(group)
)

# subset matrix + family table (keep order for now; clustering later)
mat <- mat[rownames(mat) %in% keep_all, , drop = FALSE]
all_family <- all_family %>% filter(group %in% keep_all)

cat("Kept", nrow(mat), "paths total; Poaceae kept:", length(keep_poaceae), "\n")

# ============================================================
# 4) DISTANCE, CLUSTERING (ALL SPECIES)
# ============================================================

dist_mat <- vegdist(mat, method = "jaccard")
hc <- hclust(dist_mat, method = "ward.D2")
all_order <- rownames(mat)[hc$order]

# ============================================================
# 5) FAMILY PALETTE
# ============================================================

fam_levels <- sort(unique(all_family$Family))
n_fam <- length(fam_levels)

if (requireNamespace("Polychrome", quietly = TRUE) && n_fam > 3) {
  seed_cols <- grDevices::hcl.colors(min(8, n_fam), palette = "Dark 3")
  fam_pal <- Polychrome::createPalette(n_fam, seedcolors = seed_cols, M = 1000)
} else {
  fam_pal <- grDevices::hcl.colors(n_fam, palette = "Dark 3")
}
names(fam_pal) <- fam_levels

# ============================================================
# 6) BIN NODE COLUMNS + BUILD LONG DATA FOR PLOT
# ============================================================

mat_all <- mat[all_order, , drop = FALSE]

node_idx <- as.integer(sub("^node\\.", "", colnames(mat_all)))
bin_id   <- floor((node_idx - 1) / node_bin_size) + 1
n_bins   <- max(bin_id, na.rm = TRUE)

binned_mat <- sapply(seq_len(n_bins), function(b) {
  cols <- which(bin_id == b)
  if (length(cols) == 1) mat_all[, cols] else rowMeans(mat_all[, cols, drop = FALSE])
})

binned_mat <- as.matrix(binned_mat)
colnames(binned_mat) <- as.character(seq_len(n_bins))
rownames(binned_mat) <- all_order

binned_df <- as_tibble(binned_mat, rownames = "path.name") %>%
  pivot_longer(cols = -path.name, names_to = "bin", values_to = "pav") %>%
  mutate(
    bin = as.integer(bin),
    path.name = factor(path.name, levels = rev(all_order))
  )

labels_df <- all_family %>%
  filter(group %in% all_order) %>%
  mutate(group = factor(group, levels = rev(all_order))) %>%
  transmute(
    group,
    Family,
    x     = 1,
    label = as.character(group)
  )

# ============================================================
# 7) PLOT (WHITE BACKGROUND + CLEAR ABSENCE)
# ============================================================

p <- ggplot() +
  geom_tile(
    data = binned_df,
    aes(x = bin, y = path.name, fill = pav),
    color = tile_border_color,
    size  = tile_border_size
  ) +
  scale_fill_gradient(
    low = pav_low_color,
    high = pav_high_color,
    name = "PAV"
  ) +
  
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
    # Critical: remove the default grey panel that makes absence look grey
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    plot.margin  = margin(t = 10, r = 10, b = 10, l = left_margin_pt)
  )

ggsave(output_png, p, width = 14, height = 10, dpi = 300)
cat("\nSaved heatmap to:", output_png, "\n")
