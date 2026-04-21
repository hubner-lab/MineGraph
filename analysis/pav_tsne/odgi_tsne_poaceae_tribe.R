
# Load required libraries (unique and grouped logically)
library(data.table)      # Fast data reading
library(tidyverse)       # Data wrangling
library(ape)             # Phylogenetics
library(dendextend)      # Dendrogram manipulation
library(phytools)        # Tree rooting
library(phylobase)       # Tree conversions
library(methods)         # S4 method support
library(FactoMineR)      # PCA
library(factoextra)      # PCA visualization
library(ggrepel)         # Label repulsion
library(Rtsne)           # t-SNEe

# 1. Load and preprocess data
df <- fread("odgi_matrix.tsv") %>%
  as_tibble() %>%
  select(path.name, starts_with("node.")) %>%
  mutate(path.name = sub("#1#1$", "", path.name))

# 2. Filter nodes with variance > 0
df_clean <- df %>%
  select(-path.name) %>%
  select(where(~ sd(.) > 0))

# 3. Add back path names
df_with_names <- df_clean %>%
  mutate(path.name = df$path.name)

# 4. Prepare numeric-only matrix for PCA
df_for_pca <- df_with_names %>%
  column_to_rownames("path.name") %>%
  as.data.frame()

# 5. PCA
pca_result <- PCA(df_for_pca, graph = FALSE, scale.unit = TRUE)

### PCA PLOT ####
library(ggplot2)
library(ggrepel)
library(tibble)

coords <- as.data.frame(pca_result$ind$coord) %>%
  rownames_to_column("path.name")

p <- ggplot(coords, aes(Dim.1, Dim.2)) +
  geom_point(size = 2) +
  geom_text_repel(
    aes(label = path.name),
    size = 2.5,
    max.overlaps = Inf,   # don't drop labels
    box.padding = 0.4,
    point.padding = 0.2
  ) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA of ODGI paths", x = "PC1", y = "PC2")

print(p)

ggsave("pca_labeled_five_families.png", p, width = 16, height = 12, dpi = 300)

#########

### PCA PLOT TRY 2 #####

library(tidyverse)
library(ggrepel)

# PCA coordinates
coords <- as.data.frame(pca_result$ind$coord) %>%
  rownames_to_column("path.name")

# % variance explained by each PC (FactoMineR: pca_result$eig)
pc1_var <- pca_result$eig[1, 2]
pc2_var <- pca_result$eig[2, 2]

# "Top-left -> bottom-right" score:
# top-left means PC1 low, PC2 high  => score should be HIGH there
coords <- coords %>%
  mutate(
    grad_score = as.numeric(scale((-Dim.1) + (Dim.2)))  # top-left high, bottom-right low
  )

p <- ggplot(coords, aes(Dim.1, Dim.2, color = grad_score)) +
  geom_point(size = 2) +
  scale_color_gradient(low = "navy", high = "gold") +
  labs(
    title = "PCA of ODGI paths",
    x = paste0("PC1 (", round(pc1_var, 1), "%)"),
    y = paste0("PC2 (", round(pc2_var, 1), "%)"),
    color = "Top-left →\nBottom-right"
  ) +
  theme_minimal()

print(p)
ggsave("pca_gradient_PC1_PC2.png", p, width = 10, height = 8, dpi = 300)

#################




library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(scales)

# PCA coordinates
coords <- as.data.frame(pca_result$ind$coord) %>%
  rownames_to_column("path.name") %>%
  rename(PC1 = Dim.1, PC2 = Dim.2)

# % variance explained (FactoMineR)
pc1_pct <- pca_result$eig[1, 2]
pc2_pct <- pca_result$eig[2, 2]

# Gradient score: top-left -> bottom-right (PC1 up, PC2 down)
coords <- coords %>%
  mutate(grad = as.numeric(scale(PC1)) - as.numeric(scale(PC2)))

# Always label these (pattern match to be safe)
targets_idx <- grepl("^Cyperus_mutica", coords$path.name) | grepl("^Ananas_comosus", coords$path.name)

# Label bottom extremes (most negative PC2), not top extremes
n_bottom <- 30  # adjust
bottom_idx <- order(coords$PC2)[1:n_bottom]

label_idx <- unique(c(which(targets_idx), bottom_idx))
coords <- coords %>%
  mutate(label = ifelse(row_number() %in% label_idx, path.name, NA))

# Heuristic "Poaceae bulb" center: right-side, near PC2 ~ 0 (dense cluster)
bulb <- coords %>%
  filter(PC1 >= quantile(PC1, 0.75), abs(PC2) <= quantile(abs(PC2), 0.65))
bulb_center <- data.frame(PC1 = mean(bulb$PC1), PC2 = mean(bulb$PC2))

# Plot
p <- ggplot(coords, aes(PC1, PC2)) +
  geom_point(aes(color = grad), size = 2, alpha = 0.9) +
  # Arrow from targets to bulb center (to show "far from Poaceae bulb")
  geom_segment(
    data = coords %>% filter(targets_idx),
    aes(x = PC1, y = PC2, xend = bulb_center$PC1, yend = bulb_center$PC2),
    inherit.aes = FALSE,
    linewidth = 0.4,
    arrow = arrow(length = unit(0.14, "inches"))
  ) +
  geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  scale_color_viridis_c(option = "plasma", name = "Top-left →\nBottom-right") +
  labs(
    title = "PCA of ODGI paths (bottom extremes labeled + key taxa)",
    x = sprintf("PC1 (%.1f%%)", pc1_pct),
    y = sprintf("PC2 (%.1f%%)", pc2_pct)
  ) +
  theme_minimal(base_size = 14)

print(p)
ggsave("pca_gradient_bottom_extremes.png", p, width = 13, height = 9, dpi = 400)


#################
# 6. t-SNE using PCA components
set.seed(42)

# PCA coordinates
pca_matrix <- pca_result$ind$coord

tsne_result <- Rtsne(
  pca_matrix,
  dims        = 2,
  perplexity  = 10,
  verbose     = TRUE,
  max_iter    = 1000
)

# 7. Prepare t-SNE data frame
tsne_df <- as_tibble(tsne_result$Y) %>%
  mutate(path.name = rownames(df_for_pca))



##################### FAMILY MAPPING (RAxML-tip-style) + t-SNE PLOT #####################

library(readxl)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)

xlsx_s3 <- "Supplementary_Table_S3_Full_Poales_Data.xlsx"

# ---- Make labels match RAxML_fixed_labels_192.tree tip formatting ----
# Rules:
#  - remove brackets []
#  - remove hyphens (NOT replace with "_")
#  - remove parentheses characters only (keep the inside text)
#  - replace remaining non-alphanumeric with "_"
#  - truncate to 30 characters
make_raxml_tip <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  
  x <- str_replace_all(x, "\\[|\\]", "")     # remove [ ]
  x <- str_replace_all(x, "-", "")           # remove hyphens entirely
  x <- str_replace_all(x, "\\(|\\)", "")     # remove ( )
  x <- str_replace_all(x, "[^A-Za-z0-9]", "_")# non-alnum -> _
  substr(x, 1, 30)                           # truncate to 30 chars
}

# 1) Read metadata and build lookup: raxml_tip_label -> Family
meta_family <- read_excel(xlsx_s3, sheet = "Full List") %>%
  mutate(across(everything(), as.character)) %>%
  mutate(
    # In your Excel, the full binomial string is in column "Species"
    tip_label = make_raxml_tip(Species),
    Family    = if_else(is.na(Family) | Family == "", "UnassignedFamily", Family)
  ) %>%
  select(tip_label, Family) %>%
  distinct()

# 2) Map t-SNE labels to the same RAxML-style tip label and join
tsne_df2 <- tsne_df %>%
  mutate(
    tip_label = make_raxml_tip(path.name)
  ) %>%
  left_join(meta_family, by = "tip_label") %>%
  mutate(Family = if_else(is.na(Family), "UnassignedFamily", Family))

# 3) Diagnostics (highly recommended once)
cat("\nFamily assignment summary:\n")
print(sort(table(tsne_df2$Family), decreasing = TRUE))

cat("\nExamples of unmatched labels (first 20):\n")
print(
  tsne_df2 %>%
    filter(Family == "UnassignedFamily") %>%
    distinct(path.name, tip_label) %>%
    head(20)
)

# 4) Palette (Polychrome optional)
fam_levels <- sort(unique(tsne_df2$Family))
n_fam <- length(fam_levels)

if (requireNamespace("Polychrome", quietly = TRUE) && n_fam > 3) {
  seed_cols <- grDevices::hcl.colors(min(8, n_fam), palette = "Dark 3")
  fam_pal <- Polychrome::createPalette(n_fam, seedcolors = seed_cols, M = 1000)
} else {
  fam_pal <- grDevices::hcl.colors(n_fam, palette = "Dark 3")
}
names(fam_pal) <- fam_levels

# 5) Plot + save
# Compute tight limits with a small padding
pad_x <- 0.5
pad_y <- 0.5
x_lim <- range(tsne_df2$V1, na.rm = TRUE) + c(-pad_x, pad_x)
y_lim <- range(tsne_df2$V2, na.rm = TRUE) + c(-pad_y, pad_y)

p_tsne <- ggplot(tsne_df2, aes(x = V1, y = V2, color = Family)) +
  geom_point(size = 3.2, alpha = 0.9) +  # you can still enlarge points
  scale_color_manual(values = fam_pal) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  coord_cartesian(xlim = x_lim, ylim = y_lim) +
  labs(x = "", y = "") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 13, face = "bold")
  )

print(p_tsne)
ggsave("tsne_family_new_small.png", p_tsne, width = 8, height = 5, dpi = 300)

########################################################################################
#### NICE PCA ### 
# Add Family to PCA coordinates using the SAME mapping as t-SNE
coords2 <- coords %>%
  left_join(
    tsne_df2 %>% select(path.name, Family),
    by = "path.name"
  ) %>%
  mutate(Family = if_else(is.na(Family), "UnassignedFamily", Family))


coords2 <- coords2 %>%
  rename(grad_score = grad)


coords_lab <- coords2 %>%
  arrange(desc(grad_score)) %>%
  slice_head(n = 20) %>%
  bind_rows(coords2 %>% arrange(grad_score) %>% slice_head(n = 20)) %>%
  distinct(path.name, .keep_all = TRUE)


p2 <- ggplot(coords2, aes(PC1, PC2, color = Family)) +
  geom_point(size = 2, alpha = 0.9) +
  
  geom_text_repel(
    data = coords_lab,
    aes(label = path.name),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = fam_pal) +
  
  labs(
    title = "PCA of ODGI paths (colored by Family)",
    x = paste0("PC1 (", round(pc1_var, 1), "%)"),
    y = paste0("PC2 (", round(pc2_var, 1), "%)"),
    color = "Family"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11, face = "bold")
  )

print(p2)
ggsave("pca_family_new.png", p2, width = 12, height = 9, dpi = 300)





#####



##DEC 22 END ##
################################

ggplot(tsne_df, aes(x = V1, y = V2)) +
  geom_point(color = "steelblue", size = 2) +
  geom_text_repel(aes(label = path.name), size = 2.5, max.overlaps = 100) +
  labs(title = "t-SNE of ODGI Paths", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()


ggsave("tsne_plot_2323.png", width = 10, height = 8, dpi = 300)

# ===== Replace everything from here (XLSX load) =====
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(fuzzyjoin)
library(tidyr)

# tsne_df must already exist with columns: V1, V2, path.name

# ---------------- Params for legend size control ----------------
top_k      <- 40   # keep the top tribe-groups in legend (diagnostic)
min_count  <- NULL # e.g., 5 to keep any group with >=5 samples (overrides top_k if set)
min_prop   <- NULL # e.g., 0.01 to keep groups with >=1% of points (overrides top_k if set)
fam_top    <- 8    # how many families get unique shapes in the legend

# ---------------- Load Excel + normalized keys (with tribe) ----------------
xlsx_path <- "Supplementary_Table_S3_Full_Poales_Data.xlsx"
full_list <- read_excel(xlsx_path, sheet = "Full List") %>%
  mutate(
    binomial     = str_squish(paste(Genus, Species)),
    binom_clean  = binomial |>
      str_to_lower() |>
      str_replace_all("_+", " ") |>
      str_replace_all("\\s+", " ") |>
      str_trim(),
    key_strict   = binom_clean |> str_replace_all("[^a-z0-9]", ""),
    genus_clean  = str_to_lower(Genus) |> str_replace_all("[^a-z0-9]", "")
  )

# ---------------- Clean t-SNE labels (underscores → spaces) ----------------
tsne_clean <- tsne_df %>%
  mutate(
    label_clean = path.name |>
      str_replace_all("_+", " ") |>
      str_replace_all("\\s+", " ") |>
      str_trim() |>
      str_to_lower(),
    key_strict  = label_clean |> str_replace_all("[^a-z0-9]", ""),
    genus_token = str_split_fixed(label_clean, "\\s+", 2)[,1] |> str_replace_all("[^a-z0-9]", "")
  )

# ---------------- (A) strict join (fast, exact) ----------------
m <- tsne_clean %>%
  left_join(full_list %>% select(Family, tribe, Genus, Species, binomial, key_strict),
            by = "key_strict") %>%
  mutate(match_source = if_else(!is.na(Family), "exact", NA_character_))

# ---------------- (B) fuzzy full-binomial rescue (Jaro–Winkler) ----------------
leftover <- m %>% filter(is.na(Family)) %>% select(path.name, label_clean) %>% distinct()
if (nrow(leftover) > 0) {
  fuzzy_full <- stringdist_inner_join(
    leftover,
    full_list %>% select(Family, tribe, Genus, Species, binom_clean),
    by = c("label_clean" = "binom_clean"),
    method = "jw",
    max_dist = 0.18,
    distance_col = "dist"
  ) %>%
    group_by(path.name) %>%
    slice_min(dist, with_ties = FALSE) %>%
    ungroup()
  
  m <- m %>%
    left_join(fuzzy_full %>% select(path.name, Family_f = Family, Tribe_f = tribe), by = "path.name") %>%
    mutate(
      Family = coalesce(Family, Family_f),
      tribe  = coalesce(tribe,  Tribe_f),
      match_source = if_else(is.na(match_source) & !is.na(Family), "fuzzy_full", match_source)
    ) %>%
    select(-ends_with("_f"))
}

# ---------------- (C) genus-only fuzzy rescue ----------------
leftover <- m %>% filter(is.na(Family)) %>% select(path.name, genus_token) %>% distinct()
if (nrow(leftover) > 0) {
  fuzzy_genus <- stringdist_inner_join(
    leftover,
    full_list %>% select(Family, tribe, genus_clean),
    by = c("genus_token" = "genus_clean"),
    method = "jw",
    max_dist = 0.12,
    distance_col = "dist_g"
  ) %>%
    group_by(path.name) %>%
    slice_min(dist_g, with_ties = FALSE) %>%
    ungroup()
  
  m <- m %>%
    left_join(fuzzy_genus %>% select(path.name, Family_g = Family, Tribe_g = tribe), by = "path.name") %>%
    mutate(
      Family = coalesce(Family, Family_g),
      tribe  = coalesce(tribe,  Tribe_g),
      match_source = if_else(is.na(match_source) & !is.na(Family), "fuzzy_genus", match_source)
    ) %>%
    select(-Family_g, -Tribe_g)
}

# ---------------- Diagnostics ----------------
n_total   <- nrow(m)
n_matched <- sum(!is.na(m$Family))
cat("Matched:", n_matched, "/", n_total, "\n")

# ---------------- Keep all points; build cross-family tribe-group ----------------
m <- m %>%
  mutate(
    Family = replace_na(Family, "UnassignedFamily"),
    tribe  = replace_na(tribe,  "Unknown tribe"),
    Group  = ifelse(Family == "UnassignedFamily", "Unassigned",
                    paste(Family, tribe, sep = ": "))
  )

# ---------------- Collapse legend to Top-K (or threshold) for tribes (diagnostic only) ----------------
tribe_counts <- m %>%
  filter(Group != "Unassigned") %>%
  count(Group, sort = TRUE) %>%
  mutate(prop = n / sum(n))

if (!is.null(min_count)) {
  keep_groups <- tribe_counts %>% filter(n >= min_count) %>% pull(Group)
} else if (!is.null(min_prop)) {
  keep_groups <- tribe_counts %>% filter(prop >= min_prop) %>% pull(Group)
} else {
  keep_groups <- head(tribe_counts$Group, top_k)
}

m_top <- m %>%
  mutate(
    GroupTop = case_when(
      Group == "Unassigned"      ~ "Unassigned",
      Group %in% keep_groups     ~ Group,
      TRUE                       ~ "Other tribes"
    )
  )

# ---------------- Shapes by Family (collapse to Top-F for legend readability) ----------------
fam_counts <- m %>% count(Family, sort = TRUE)
keep_fams  <- head(fam_counts$Family, fam_top)

m_top <- m_top %>%
  mutate(
    FamilyTop = case_when(
      Family == "UnassignedFamily" ~ "Unassigned",
      Family %in% keep_fams        ~ Family,
      TRUE                         ~ "Other families"
    ),
    # helpful ordering for prettier legends
    GroupTop  = factor(GroupTop,
                       levels = c(sort(setdiff(unique(GroupTop), c("Other tribes","Unassigned"))),
                                  "Other tribes","Unassigned")),
    FamilyTop = factor(FamilyTop,
                       levels = c(sort(setdiff(unique(FamilyTop), c("Other families","Unassigned"))),
                                  "Other families","Unassigned"))
  )

# ---------------- Family palette (for non-Poaceae) ----------------
# Build a family color palette (prefer RColorBrewer Set3; fallback to hue)
fam_levels <- levels(m_top$FamilyTop)
n_fams     <- length(fam_levels)

if (requireNamespace("RColorBrewer", quietly = TRUE)) {
  max_set3 <- min(max(n_fams, 3), 12)
  fam_cols <- RColorBrewer::brewer.pal(max_set3, "Set3")
  if (n_fams > max_set3) {
    extra <- scales::hue_pal()(n_fams - max_set3)
    fam_pal_vec <- c(fam_cols, extra)
  } else {
    fam_pal_vec <- fam_cols[seq_len(n_fams)]
  }
} else {
  fam_pal_vec <- scales::hue_pal()(n_fams)
}
names(fam_pal_vec) <- fam_levels

# ---------------- Poaceae tribe spectrum + per-tribe labels ----------------
# identify Poaceae rows using the original Family column (not collapsed)
poa_idx <- which(m_top$Family == "Poaceae")
poa_tribes <- m_top$tribe[poa_idx]
poa_tribes <- poa_tribes[poa_tribes != "Unknown tribe"]
poa_tribe_levels <- sort(unique(poa_tribes))

# spectrum for tribes (prefer viridisLite::turbo; fallback to hue)
if (requireNamespace("viridisLite", quietly = TRUE)) {
  poa_cols <- viridisLite::turbo(max(length(poa_tribe_levels), 3))[seq_along(poa_tribe_levels)]
} else {
  poa_cols <- scales::hue_pal()(max(length(poa_tribe_levels), 3))[seq_along(poa_tribe_levels)]
}
tribe_col_map <- setNames(poa_cols, poa_tribe_levels)

# Construct final color per point:
# - default to FamilyTop color
# - override to tribe color for Poaceae rows
base_colors <- fam_pal_vec[as.character(m_top$FamilyTop)]
if (length(poa_idx) > 0) {
  base_colors[poa_idx] <- tribe_col_map[ m_top$tribe[poa_idx] ]
}
m_top$ColorFinal <- unname(base_colors)

# Build tribe centroid labels for Poaceae (avoid clutter with a minimum size)
min_points_per_tribe <- 3
label_df <- m_top %>%
  filter(Family == "Poaceae", tribe %in% names(tribe_col_map)) %>%
  group_by(tribe) %>%
  summarize(
    x = median(V1, na.rm = TRUE),
    y = median(V2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= min_points_per_tribe) %>%
  mutate(label_col = unname(tribe_col_map[tribe]))

# ---------------- Shapes for families ----------------
base_shapes <- c(16,17,15,18,7,8,3,4,0,2,5,6,9,10,11)  # plenty for ~10 families
if (length(base_shapes) < n_fams) {
  extra_pool  <- setdiff(0:25, base_shapes)
  base_shapes <- c(base_shapes, extra_pool[seq_len(n_fams - length(base_shapes))])
}
fam_shapes <- base_shapes[seq_len(n_fams)]
names(fam_shapes) <- fam_levels

# ---------------- Final plot ----------------
# - Points: color = ColorFinal (identity), shape = FamilyTop
# - Tribe labels for Poaceae: placed at centroids, colored by tribe (identity), no legend
p <- ggplot(m_top, aes(V1, V2)) +
  geom_point(aes(color = ColorFinal, shape = FamilyTop),
             size = 2.2, alpha = 0.95, stroke = 0.6, show.legend = TRUE) +
  # Poaceae tribe labels (no legend; uses the same identity color scale)
  ggrepel::geom_text_repel(
    data = label_df,
    aes(x = x, y = y, label = tribe, color = label_col),
    inherit.aes = FALSE,
    size = 3.2, fontface = "bold",
    box.padding = 0.3, point.padding = 0.2,
    segment.size = 0.2, min.segment.length = 0.1,
    show.legend = FALSE, max.overlaps = Inf
  ) +
  scale_color_identity(guide = "none") +  # suppress huge color legend
  scale_shape_manual(values = fam_shapes, name = "Family") +
  labs(
    title = "t-SNE of ODGI Paths — Poaceae tribes (spectrum) & Families (shape)",
    x = "t-SNE 1", y = "t-SNE 2"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.text  = element_text(size = 9),
    legend.title = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 4), ncol = 1)
  )

print(p)
ggsave("tsne_poaceae_tribes_spectrum_with_family_shapes.png", p, width = 11, height = 8.5, dpi = 300)

# ===== End replacement =====


# 8. Hierarchical clustering (no k-means)
hclust_ward <- function(x) hclust(x, method = "ward.D2")
Z <- hclust_ward(dist(df_for_pca))
####################################
# UMAP (uses PCA as input)
library(uwot)

set.seed(42)
umap_result <- umap(pca_matrix, init = "spectral", n_neighbors = 15, min_dist = 0.1)
umap_df <- as_tibble(umap_result) %>%
  mutate(path.name = rownames(df_for_pca))

ggplot(umap_df, aes(x = V1, y = V2)) +
  geom_point(color = "forestgreen", size = 2) +
  geom_text_repel(aes(label = path.name), size = 2.5, max.overlaps = 100) +
  labs(title = "UMAP (PCA-based)", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()

ggsave("umap_plot.png", width = 10, height = 8, dpi = 300)
###################################

###############



##############

# Bootstrap support
n_bootstraps <- 1
bootstrap_support <- numeric(length(Z$height))

for (i in 1:n_bootstraps) {
  boot_idx <- sample(nrow(df_clean), replace = TRUE)
  boot_data <- df_clean[boot_idx, ]
  Z_boot <- hclust_ward(dist(boot_data))
  bootstrap_support <- bootstrap_support + (Z$height <= quantile(Z_boot$height, 0.95, na.rm = TRUE))
}
bootstrap_support <- (bootstrap_support / n_bootstraps) * 100

# Dendrogram
dend <- as.dendrogram(Z)
dend <- dendextend::set(dend, "branches_lwd", 1.5)
dend <- dendextend::set(dend, "labels_cex", 0.6)

png("bootstrapped_dendrogram_222.png", width = 1200, height = 800, res = 150)
par(mar = c(5, 12, 2, 2))
plot(dend,
     main = paste("Bootstrapped Dendrogram (n =", n_bootstraps, ")"),
     xlab = "", ylab = "Distance", cex.axis = 0.8)
dev.off()

# Convert to phylo tree
tree_phylo <- as.phylo(Z)

#####
# 1) Convert your graph clustering (hclust) to phylo
graph_phy <- as.phylo(Z)

# 2) Ensure labels match your cleaned path names (IMPORTANT)
#    df_for_pca rownames are the names used in PCA/t-SNE, already stripped of "#1#1"
graph_phy$tip.label <- rownames(df_for_pca)

# (Optional) If you want RAxML-style formatting for SplitsTree matching:
# graph_phy$tip.label <- make_raxml_tip(graph_phy$tip.label)

# 3) Root it (SplitsTree works best with rooted trees)
if (!ape::is.rooted(graph_phy)) {
  graph_phy <- phytools::midpoint.root(graph_phy)
}

# 4) Write Newick
ape::write.tree(graph_phy, file = "50_Graph_tree_from_ODGI.nwk")

cat("Saved graph tree to: Graph_tree_from_ODGI.nwk\n")
#####

plot(tree_phylo, main = paste("Phylogenetic Tree with Bootstrap (n =", n_bootstraps, ")"),
     cex = 0.6)
nodelabels(text = round(bootstrap_support), cex = 0.5, frame = "none")
#################
#deal with names
normalize_labels <- function(labels) {
  labels %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]", "")  # keep only alphanumerics
}

# 8. Tree Comparison Function (graph tree vs RAxML) — with faint/transparent low-gradient lines
# 8. Tree Comparison Function (graph tree vs RAxML) — with faint/transparent low-gradient lines
# 8. Tree Comparison Function (graph tree vs RAxML) — with faint/transparent low-gradient lines
safe_tree_comparison <- function(hclust_tree, raxml_file,
                                 k_branches = 5,
                                 alpha_power = 2.0,
                                 alpha_floor = 0.02,
                                 alpha_threshold = 0,
                                 lwd_min = 0.3,
                                 lwd_max = 4.0) {
  require(ape)
  require(phytools)
  require(dendextend)
  require(grDevices)
  require(readxl)
  require(stringr)
  
  phylo_to_dend_safe <- function(phy) {
    d <- cophenetic(phy)
    hcl <- hclust(as.dist(d), method = "average")
    as.dendrogram(hcl)
  }
  
  if (!file.exists(raxml_file)) stop("RAxML tree file not found: ", raxml_file)
  
  # --- Graph-based tree (hclust) ---
  dend_hclust <- as.dendrogram(hclust_tree)
  dend_hclust <- dendextend::set(dend_hclust, "branches_lwd", 1.5)
  dend_hclust <- dendextend::set(dend_hclust, "labels_cex", 0.6)
  
  # --- RAxML tree ---
  dend_raxml <- tryCatch({
    raxml_tree <- ape::read.tree(raxml_file)
    if (!ape::is.rooted(raxml_tree)) {
      message("RAxML tree is unrooted - applying midpoint rooting")
      raxml_tree <- phytools::midpoint.root(raxml_tree)
    }
    dend_tmp <- phylo_to_dend_safe(raxml_tree)
    dend_tmp <- dendextend::set(dend_tmp, "branches_lwd", 1.5)
    dend_tmp <- dendextend::set(dend_tmp, "labels_cex", 0.6)
    dend_tmp
  }, error = function(e) {
    warning("RAxML tree processing failed: ", e$message)
    return(NULL)
  })
  if (is.null(dend_raxml)) return(NULL)
  
  # ============================================================
  # NAME MAPPING (robust)  ——  APPLY WITH dendextend::set()
  # ============================================================
  
  norm_name <- function(x) {
    x <- as.character(x)
    x <- stringr::str_trim(x)
    x <- stringr::str_replace_all(x, "\\[|\\]", "")
    x <- stringr::str_replace_all(x, "\\(|\\)", "")
    x <- stringr::str_replace_all(x, "-", "")
    x <- stringr::str_replace_all(x, "[^A-Za-z0-9]", "_")
    x <- stringr::str_replace_all(x, "_{2,}", "_")
    x <- stringr::str_replace_all(x, "^_|_$", "")
    x
  }
  
  base_key <- function(x) {
    x <- norm_name(x)
    toks <- unlist(strsplit(x, "_", fixed = TRUE))
    toks <- toks[toks != ""]
    if (length(toks) >= 2) paste(toks[1], toks[2], sep = "_") else x
  }
  
  str_dist <- function(a, b) {
    if (requireNamespace("stringdist", quietly = TRUE)) {
      stringdist::stringdist(a, b, method = "jw")
    } else {
      d <- utils::adist(a, b)
      as.numeric(d) / max(1, max(nchar(a), nchar(b)))
    }
  }
  
  build_one_to_one_map <- function(h_labs, r_labs) {
    h_labs <- as.character(h_labs)
    r_labs <- as.character(r_labs)
    
    h_norm <- norm_name(h_labs)
    r_norm <- norm_name(r_labs)
    
    h_key <- vapply(h_norm, base_key, character(1))
    r_key <- vapply(r_norm, base_key, character(1))
    
    map <- setNames(rep(NA_character_, length(h_labs)), h_labs)
    used_r <- setNames(rep(FALSE, length(r_labs)), r_labs)
    
    keys <- intersect(unique(h_key), unique(r_key))
    
    for (k in keys) {
      hi <- which(h_key == k)
      ri <- which(r_key == k)
      if (!length(hi) || !length(ri)) next
      
      pairs <- expand.grid(hi = hi, ri = ri, KEEP.OUT.ATTRS = FALSE)
      pairs$h_lab  <- h_labs[pairs$hi]
      pairs$r_lab  <- r_labs[pairs$ri]
      pairs$h_norm <- h_norm[pairs$hi]
      pairs$r_norm <- r_norm[pairs$ri]
      pairs$dist   <- mapply(str_dist, pairs$h_norm, pairs$r_norm)
      
      pairs <- pairs[order(pairs$dist), ]
      taken_h <- setNames(rep(FALSE, length(hi)), h_labs[hi])
      taken_r <- setNames(rep(FALSE, length(ri)), r_labs[ri])
      
      for (i in seq_len(nrow(pairs))) {
        hname <- pairs$h_lab[i]
        rname <- pairs$r_lab[i]
        if (isTRUE(taken_h[hname])) next
        if (isTRUE(taken_r[rname])) next
        if (isTRUE(used_r[rname]))  next
        
        map[hname] <- rname
        taken_h[hname] <- TRUE
        taken_r[rname] <- TRUE
        used_r[rname]  <- TRUE
      }
    }
    
    map
  }
  
  h0 <- as.character(labels(dend_hclust))
  r0 <- as.character(labels(dend_raxml))
  
  h_to_r <- build_one_to_one_map(h0, r0)
  matched <- sum(!is.na(h_to_r))
  message("Name mapping: matched ", matched, " / ", length(h0), " labels")
  
  # APPLY AUTO RENAMING
  new_h <- h0
  idx_m <- which(!is.na(h_to_r))
  new_h[idx_m] <- unname(h_to_r[h0[idx_m]])
  
  # ---------- MANUAL FALLBACK (only for remaining unmatched) ----------
  manual_map <- c(
    "Coix_lacrymajobi"               = "Coix_lacryma_jobi",
    "Echinochloa_crusgalli_var__pra" = "Echinochloa_crus_galli_var_praticola",
    "Taeniatherum_caputmedusae_subs" = "Taeniatherum_caput_medusae_subsp_caput_medusae"
  )
  
  # apply only if target exists in RAxML labels
  unmatched_now <- setdiff(new_h, r0)
  for (u in intersect(unmatched_now, names(manual_map))) {
    tgt <- manual_map[[u]]
    if (tgt %in% r0) new_h[new_h == u] <- tgt
  }
  # -------------------------------------------------------------------
  
  dend_hclust <- dendextend::set(dend_hclust, "labels", as.character(new_h))
  
  exact_after <- sum(labels(dend_hclust) %in% labels(dend_raxml))
  message("Exact label matches after renaming: ", exact_after, " / ", length(labels(dend_hclust)))
  
  # Now strict intersection works
  common_labels <- intersect(as.character(labels(dend_hclust)), as.character(labels(dend_raxml)))
  if (length(common_labels) < 3) {
    warning("Only ", length(common_labels), " common labels between trees after mapping")
    return(NULL)
  }
  
  dend1 <- dendextend::prune(dend_hclust, setdiff(labels(dend_hclust), common_labels))
  dend2 <- dendextend::prune(dend_raxml,  setdiff(labels(dend_raxml),  common_labels))
  
  common_labels <- unique(common_labels[!is.na(common_labels) & common_labels != ""])
  dend1 <- dendextend::rotate(dend1, common_labels)
  dend2 <- dendextend::rotate(dend2, common_labels)
  
  assign("dend1", dend1, envir = .GlobalEnv)
  assign("dend2", dend2, envir = .GlobalEnv)
  
  #################################################################
  # FAMILY-BASED LABEL COLORS — keep your logic unchanged
  #################################################################
  
  make_raxml_tip <- function(x) {
    x <- as.character(x)
    x <- stringr::str_trim(x)
    x <- stringr::str_replace_all(x, "\\[|\\]", "")
    x <- stringr::str_replace_all(x, "-", "")
    x <- stringr::str_replace_all(x, "\\(|\\)", "")
    x <- stringr::str_replace_all(x, "[^A-Za-z0-9]", "_")
    x <- stringr::str_replace_all(x, "_{2,}", "_")
    x <- stringr::str_replace_all(x, "^_|_$", "")
    substr(x, 1, 30)
  }
  
  xlsx_s3 <- "Supplementary_Table_S3_Full_Poales_Data.xlsx"
  if (!file.exists(xlsx_s3)) {
    warning("Family file not found — skipping family coloring")
  } else {
    taxa_meta <- suppressMessages(
      readxl::read_excel(xlsx_s3, sheet = "Full List")
    ) %>%
      dplyr::mutate(across(everything(), as.character)) %>%
      dplyr::mutate(
        tip_label = make_raxml_tip(Species),
        Family    = dplyr::if_else(is.na(Family) | Family == "", "UnassignedFamily", Family)
      ) %>%
      dplyr::select(tip_label, Family) %>%
      dplyr::distinct()
    
    if (is.null(fam_pal)) {
      fam_levels <- sort(unique(taxa_meta$Family))
      n_fam <- length(fam_levels)
      if (requireNamespace("Polychrome", quietly = TRUE) && n_fam > 3) {
        seed_cols <- grDevices::hcl.colors(min(8, n_fam), palette = "Dark 3")
        fam_pal <- Polychrome::createPalette(n_fam, seedcolors = seed_cols, M = 1000)
      } else {
        fam_pal <- grDevices::hcl.colors(n_fam, palette = "Dark 3")
      }
      names(fam_pal) <- fam_levels
    }
    
    labs1 <- labels(dend1)
    labs2 <- labels(dend2)
    
    label_df <- tibble::tibble(label = unique(c(labs1, labs2))) %>%
      dplyr::left_join(taxa_meta, by = c("label" = "tip_label")) %>%
      dplyr::mutate(
        Family = dplyr::if_else(is.na(Family), "UnassignedFamily", Family),
        color  = fam_pal[Family]
      )
    
    dend1 <- dendextend::set(dend1, "labels_col", label_df$color[match(labs1, label_df$label)])
    dend2 <- dendextend::set(dend2, "labels_col", label_df$color[match(labs2, label_df$label)])
  }
  
  # --- Displacement → line color & width (same as before) ---
  labels1 <- labels(dend1)
  labels2 <- labels(dend2)
  match_idx <- match(labels1, labels2)
  displacement <- abs(seq_along(labels1) - match_idx)
  max_disp <- max(displacement)
  if (max_disp == 0) max_disp <- 1
  
  base_cols <- colorRampPalette(c("blue", "red"))(max_disp + 1)
  base_cols <- base_cols[displacement + 1]
  disp_norm <- displacement / max_disp
  
  show_mask   <- disp_norm > alpha_threshold
  disp_scaled <- ifelse(show_mask, (disp_norm - alpha_threshold) / (1 - alpha_threshold), 0)
  alpha_vals  <- disp_scaled ^ alpha_power
  alpha_vals  <- ifelse(show_mask, pmax(alpha_vals, alpha_floor), 0)
  
  line_colors <- mapply(function(col, a) adjustcolor(col, alpha.f = a), base_cols, alpha_vals)
  line_widths <- lwd_min + (lwd_max - lwd_min) * alpha_vals
  
  tg_file <- "50_tanglegram_comparison_RAxML_gradient21223.png"
  ok <- TRUE
  tryCatch({
    png(tg_file, width = 12, height = 20, units = "in", res = 300)
    tanglegram(
      dend1, dend2,
      highlight_distinct_edges = TRUE,
      common_subtrees_color_lines = TRUE,
      color_lines = line_colors,
      lwd = line_widths,
      lab.cex = 0.8,
      margin_inner = 12,
      columns_width = c(5, 2, 5),
      edge.lwd = 0.5
    )
    legend(
      x = 0.5, y = -0.02,
      xpd = NA,
      legend = c("Low displacement", "High displacement"),
      col = c(adjustcolor("blue", alpha.f = 0.2),
              adjustcolor("red",  alpha.f = 1.0)),
      lwd = c(lwd_min, lwd_max),
      horiz = TRUE,
      bty = "n",
      xjust = 0.5
    )
    dev.off()
  }, error = function(e) {
    ok <<- FALSE
    if (dev.cur() > 1) try(dev.off(), silent = TRUE)
    warning("Failed to save tanglegram: ", e$message)
  })
  if (!ok) return(NULL)
  
  entanglement_score <- entanglement(dendlist(dend1, dend2))
  
  dend_to_phylo_safe <- function(dend) {
    # ape exports as.phylo.dendrogram, so this works if 'ape' is loaded
    phy <- ape::as.phylo(dend)
    phy$tip.label <- as.character(phy$tip.label)
    phy
  }
  
  
  save_fixed <- TRUE
  
  if (save_fixed) {
    phy1_fixed <- dend_to_phylo_safe(dend1)
    phy2_fixed <- dend_to_phylo_safe(dend2)
    
    ape::write.tree(phy1_fixed, file = "50_hclust_fixed_names_new.nwk")
    ape::write.tree(phy2_fixed, file = "50_raxml_fixed_names_new.nwk")
    
    png("hclust_fixed_names.png", width = 2800, height = 4200, res = 400)
    plot(phy1_fixed, cex = 0.35, no.margin = TRUE)
    dev.off()
    
    png("raxml_fixed_names.png", width = 2800, height = 4200, res = 400)
    plot(phy2_fixed, cex = 0.35, no.margin = TRUE)
    dev.off()
    
    message("Saved fixed-name trees: hclust_fixed_names.nwk/.png and raxml_fixed_names.nwk/.png")
  }
  
  
  return(list(
    entanglement = entanglement_score,
    total_labels_hclust = length(labels(dend_hclust)),
    total_labels_raxml  = length(labels(dend_raxml)),
    common_labels       = length(common_labels),
    unique_to_hclust    = length(setdiff(labels(dend_hclust), common_labels)),
    unique_to_raxml     = length(setdiff(labels(dend_raxml), common_labels)),
    percent_matched     = round(
      100 * length(common_labels) /
        min(length(labels(dend_hclust)), length(labels(dend_raxml))), 2
    ),
    plot_file              = tg_file,
    unique_to_hclust_names = setdiff(labels(dend_hclust), common_labels),
    unique_to_raxml_names  = setdiff(labels(dend_raxml), common_labels),
    hclust_names           = labels(dend_hclust),
    raxml_names            = labels(dend_raxml)
  ))
}


# 9. Run tree comparison
result <- safe_tree_comparison(Z, "50_supermatrix_partitions.treefile")
print(result)

if (!is.null(result)) {
  cat("Entanglement score:", result$entanglement, "\n")
  cat("Common labels:", result$common_labels, "\n")
  if (!is.null(result$plot_file)) {
    cat("Tanglegram saved as:", result$plot_file, "\n")
  }
} else {
  message("Comparison failed - check warnings for details")
}

############################## DEC 24 START ###############

#!/usr/bin/env Rscript
# ------------------------------------------------------------
# POST-PROCESSING MENU for your tanglegram result
# Run this AFTER you already ran:
#   result <- safe_tree_comparison(Z, "RAxML_fixed_labels_192.tree")
#   print(result)
#
# It will:
#  1) Reconstruct dend1/dend2 (they were assign()'d in your function)
#  2) Compute a displacement table
#  3) Save multiple figure "options" you can choose from:
#     - full (only high displacement lines)
#     - pruned/zoom panels around top displaced tips
#     - a tip list you can paste into a caption/supp table
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dendextend)
  library(dplyr)
  library(stringr)
  library(grDevices)
})

# ----------------------------
# 0) Safety: confirm dend1/dend2 exist
# ----------------------------
if (!exists("dend1", envir = .GlobalEnv) || !exists("dend2", envir = .GlobalEnv)) {
  stop("I can't find `dend1` and `dend2` in .GlobalEnv.\n",
       "Make sure you ran safe_tree_comparison() first.\n",
       "Your function must have these lines (it already does):\n",
       "  assign('dend1', dend1, envir = .GlobalEnv)\n",
       "  assign('dend2', dend2, envir = .GlobalEnv)\n")
}

dend1 <- get("dend1", envir = .GlobalEnv)
dend2 <- get("dend2", envir = .GlobalEnv)

# ----------------------------
# 1) Displacement table
# ----------------------------
labs1 <- labels(dend1)
labs2 <- labels(dend2)

match_idx <- match(labs1, labs2)
if (any(is.na(match_idx))) {
  warning("Some labels in dend1 were not found in dend2 AFTER pruning/rotation (unexpected).")
}

displacement <- abs(seq_along(labs1) - match_idx)
max_disp <- max(displacement, na.rm = TRUE)
if (max_disp == 0) max_disp <- 1

disp_df <- tibble(
  label       = labs1,
  pos_tree1   = seq_along(labs1),
  pos_tree2   = match_idx,
  displacement = displacement,
  disp_norm   = displacement / max_disp
) %>%
  arrange(desc(displacement))

write.csv(disp_df, "tanglegram_displacement_table.csv", row.names = FALSE)

cat("\nSaved: tanglegram_displacement_table.csv\n")
cat("Top 20 displaced tips:\n")
print(disp_df %>% slice_head(n = 20))

# ----------------------------
# 2) Helper functions
# ----------------------------

# Build line colors/widths for a given pair of dendrograms
make_line_style <- function(d1, d2,
                            alpha_threshold = 0.20,
                            alpha_power = 2.2,
                            alpha_floor = 0.00,
                            lwd_min = 0.2,
                            lwd_max = 4.0,
                            low_col = "blue",
                            high_col = "red") {
  l1 <- labels(d1)
  l2 <- labels(d2)
  
  mi <- match(l1, l2)
  disp <- abs(seq_along(l1) - mi)
  md <- max(disp, na.rm = TRUE)
  if (md == 0) md <- 1
  
  disp_norm <- disp / md
  
  show_mask <- disp_norm > alpha_threshold
  disp_scaled <- ifelse(show_mask, (disp_norm - alpha_threshold) / (1 - alpha_threshold), 0)
  
  alpha_vals <- disp_scaled ^ alpha_power
  alpha_vals <- ifelse(show_mask, pmax(alpha_vals, alpha_floor), 0)
  
  base_cols <- colorRampPalette(c(low_col, high_col))(md + 1)
  base_cols <- base_cols[disp + 1]
  
  line_colors <- mapply(function(col, a) adjustcolor(col, alpha.f = a), base_cols, alpha_vals)
  line_widths <- lwd_min + (lwd_max - lwd_min) * alpha_vals
  
  list(colors = line_colors, widths = line_widths, disp = disp, disp_norm = disp_norm)
}

# Choose labels for a "zoom" window: top displaced + ± neighbors on each side (in tree1 order)
zoom_subset_labels <- function(dend, top_labels, k_side = 20) {
  labs <- labels(dend)
  idx  <- match(top_labels, labs)
  idx  <- idx[!is.na(idx)]
  keep_idx <- unique(unlist(lapply(idx, function(i) {
    lo <- max(1, i - k_side)
    hi <- min(length(labs), i + k_side)
    lo:hi
  })))
  labs[keep_idx]
}

# Prune both trees to the same label set and rotate tree2 to tree1 order
prune_pair <- function(d1, d2, keep_labels) {
  d1p <- dendextend::prune(d1, setdiff(labels(d1), keep_labels))
  d2p <- dendextend::prune(d2, setdiff(labels(d2), keep_labels))
  
  # Make sure both have identical label sets
  common <- intersect(labels(d1p), labels(d2p))
  d1p <- dendextend::prune(d1p, setdiff(labels(d1p), common))
  d2p <- dendextend::prune(d2p, setdiff(labels(d2p), common))
  
  # Align ordering
  d1p <- rotate(d1p, labels(d1p))
  d2p <- rotate(d2p, labels(d1p))
  
  list(d1 = d1p, d2 = d2p)
}

save_tangle <- function(filename, d1, d2,
                        alpha_threshold = 0.20,
                        alpha_power = 2.2,
                        alpha_floor = 0.00,
                        lwd_min = 0.2,
                        lwd_max = 4.0,
                        width_in = 12,
                        height_in = 20,
                        res = 900,
                        lab_cex = 0.75,
                        margin_inner = 10,
                        columns_width = c(5, 2, 5),
                        add_legend = TRUE) {
  style <- make_line_style(
    d1, d2,
    alpha_threshold = alpha_threshold,
    alpha_power = alpha_power,
    alpha_floor = alpha_floor,
    lwd_min = lwd_min,
    lwd_max = lwd_max
  )
  
  png(filename, width = width_in, height = height_in, units = "in", res = res)
  tanglegram(
    d1, d2,
    highlight_distinct_edges = TRUE,
    common_subtrees_color_lines = TRUE,
    color_lines = style$colors,
    lwd = style$widths,
    lab.cex = lab_cex,
    margin_inner = margin_inner,
    columns_width = columns_width,
    edge.lwd = 0.5
  )
  
  if (add_legend) {
    legend(
      x = 0.5, y = -0.02,
      xpd = NA,
      legend = c("Low displacement", "High displacement"),
      col = c(adjustcolor("blue", alpha.f = 0.2),
              adjustcolor("red", alpha.f = 1.0)),
      lwd = c(lwd_min, lwd_max),
      horiz = TRUE,
      bty = "n",
      xjust = 0.5
    )
  }
  
  dev.off()
  invisible(style)
}

# ----------------------------
# 3) OUTPUT OPTIONS (you'll get several PNGs)
# ----------------------------
dir.create("tanglegram_options", showWarnings = FALSE)
setwd("tanglegram_options")

cat("\nWriting outputs into: ", normalizePath("."), "\n", sep = "")

# Option 1: FULL but only show strong discordances (good for main paper panel)
cat("\n[Option 1] Full tanglegram with only high-displacement lines visible...\n")
save_tangle(
  filename = "OPTION1_full_high_discordance.png",
  d1 = dend1, d2 = dend2,
  alpha_threshold = 0.20,   # <- raise to 0.30 for even less clutter
  alpha_power = 2.2,
  alpha_floor = 0.00,
  lab_cex = 0.70,
  width_in = 12, height_in = 20, res = 900
)

# Option 2: FULL but slightly more permissive (for Supplementary, still cleaner than all-lines)
cat("\n[Option 2] Full tanglegram with medium discordance visibility...\n")
save_tangle(
  filename = "OPTION2_full_medium_discordance.png",
  d1 = dend1, d2 = dend2,
  alpha_threshold = 0.10,
  alpha_power = 2.0,
  alpha_floor = 0.02,
  lab_cex = 0.70,
  width_in = 12, height_in = 20, res = 900
)

# Option 3: ZOOM around top displaced tips (objective zoom panel)
topN <- 15
k_side <- 20
top_labels <- disp_df$label[1:topN]
keep_labels <- zoom_subset_labels(dend1, top_labels, k_side = k_side)
pair_zoom <- prune_pair(dend1, dend2, keep_labels)

cat("\n[Option 3] Zoom panel around top ", topN,
    " displaced tips (+/- ", k_side, " neighbors)...\n", sep = "")
save_tangle(
  filename = sprintf("OPTION3_zoom_top%d_k%d.png", topN, k_side),
  d1 = pair_zoom$d1, d2 = pair_zoom$d2,
  alpha_threshold = 0.00,  # in zoom we can show all, it’s not too crowded
  alpha_power = 2.0,
  alpha_floor = 0.05,
  lab_cex = 0.85,
  width_in = 8, height_in = 10, res = 900,
  margin_inner = 8,
  columns_width = c(5, 2, 5)
)

# Option 4: ZOOM of the NEXT set (ranks 16–30) (second zoom panel)
topN2 <- 15
offset <- 15
labels2 <- disp_df$label[(offset + 1):(offset + topN2)]
labels2 <- labels2[!is.na(labels2)]
keep2 <- zoom_subset_labels(dend1, labels2, k_side = k_side)
pair_zoom2 <- prune_pair(dend1, dend2, keep2)

cat("\n[Option 4] Second zoom panel (displacement ranks ",
    offset + 1, "–", offset + topN2, ")...\n", sep = "")
save_tangle(
  filename = sprintf("OPTION4_zoom_rank%dto%d_k%d.png", offset + 1, offset + topN2, k_side),
  d1 = pair_zoom2$d1, d2 = pair_zoom2$d2,
  alpha_threshold = 0.00,
  alpha_power = 2.0,
  alpha_floor = 0.05,
  lab_cex = 0.85,
  width_in = 8, height_in = 10, res = 900,
  margin_inner = 8,
  columns_width = c(5, 2, 5)
)

# Option 5: PRUNED discordant-only subset (very clean “main figure”)
# Keep top 25 most displaced tips (no neighbors), purely discordant subset
top_only <- 25
keep_only <- disp_df$label[1:top_only]
pair_pruned <- prune_pair(dend1, dend2, keep_only)

cat("\n[Option 5] Pruned tanglegram showing ONLY the top ", top_only, " displaced taxa...\n", sep = "")
save_tangle(
  filename = sprintf("OPTION5_pruned_top%d_only.png", top_only),
  d1 = pair_pruned$d1, d2 = pair_pruned$d2,
  alpha_threshold = 0.00,
  alpha_power = 2.0,
  alpha_floor = 0.10,
  lab_cex = 0.95,
  width_in = 7, height_in = 7, res = 900,
  margin_inner = 6,
  columns_width = c(5, 2, 5)
)

# Save lists for caption/supp table
writeLines(top_labels, "top_displaced_labels_rank1to15.txt")
writeLines(labels2, "top_displaced_labels_rank16to30.txt")

cat("\nDONE.\n")
cat("Created PNG options in folder: tanglegram_options/\n")
cat("Also created:\n")
cat("  - tanglegram_displacement_table.csv (one folder above)\n")
cat("  - top_displaced_labels_rank1to15.txt\n")
cat("  - top_displaced_labels_rank16to30.txt\n\n")

# Quick summary of what you now have
cat("Figure options you can pick:\n")
cat("  OPTION1_full_high_discordance.png  (main paper candidate)\n")
cat("  OPTION2_full_medium_discordance.png (supp or backup)\n")
cat(sprintf("  OPTION3_zoom_top%d_k%d.png (zoom panel 1)\n", topN, k_side))
cat(sprintf("  OPTION4_zoom_rank%dto%d_k%d.png (zoom panel 2)\n", offset + 1, offset + topN2, k_side))
cat(sprintf("  OPTION5_pruned_top%d_only.png (discordant subset panel)\n", top_only))




###### DEC 24 END ######



library(stringdist)
library(ape)

# Let's say you have these from your `safe_tree_comparison()`:
hclust_names <- result$hclust_names
raxml_names <- result$raxml_names

# 1. Build a best-name mapping using minimal string distance
name_mapping <- sapply(hclust_names, function(name_h) {
  distances <- stringdist::stringdist(name_h, raxml_names, method = "jw")  # Jaro-Winkler is good for fuzzy matches
  raxml_names[which.min(distances)]
}, USE.NAMES = TRUE)

# Optional: Show mapping table
mapping_df <- data.frame(
  hclust = names(name_mapping),
  matched_raxml = unname(name_mapping)
)
print(mapping_df)

# 2. Fix the RAxML tree by renaming the tips
fix_raxml_tree_labels <- function(tree, name_mapping) {
  tree$tip.label <- sapply(tree$tip.label, function(old_name) {
    mapped <- names(name_mapping)[which(name_mapping == old_name)]
    if (length(mapped) > 0) mapped[1] else old_name
  })
  return(tree)
}

# Example:
tree_raxml <- read.tree("raxml_192_output.raxml_192_fixed.bestTree")
tree_raxml_fixed <- fix_raxml_tree_labels(tree_raxml, name_mapping)
plot(tree_raxml_fixed)

write.tree(tree_raxml_fixed, file = "RAxML_fixed_labels_192.tree")
