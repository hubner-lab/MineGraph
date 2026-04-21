#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)

# -----------------------------
# Input
# -----------------------------
infile  <- "gene_PAV_matrix_var_gt0.tsv"
outfile <- "gene_PAV_heatmap.png"

# -----------------------------
# Read data
# -----------------------------
mat <- fread(infile)

meta_cols <- c("species", "Family")
gene_cols <- setdiff(names(mat), meta_cols)

# Ensure numeric (safety)
mat[, (gene_cols) := lapply(.SD, as.numeric), .SDcols = gene_cols]

# -----------------------------
# Optional: order species by Family
# -----------------------------
mat <- mat %>%
  arrange(Family, species)

# -----------------------------
# Long format for ggplot
# -----------------------------
long <- mat %>%
  pivot_longer(
    cols = all_of(gene_cols),
    names_to = "gene",
    values_to = "PAV"
  )

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(long, aes(x = gene, y = species, fill = factor(PAV))) +
  geom_tile(color = NA) +
  scale_fill_manual(
    values = c("0" = "red", "1" = "blue"),
    name   = "PAV"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y  = element_text(size = 6),
    panel.grid   = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 10)
  ) +
  labs(
    x = "Gene",
    y = "Species",
    title = "Presence–Absence Variation (Variable Genes Only)"
  )

# -----------------------------
# Save
# -----------------------------
ggsave(
  outfile,
  plot   = p,
  width  = max(8, length(gene_cols) * 0.4),
  height = max(14, nrow(mat) * 0.05),
  dpi    = 300
)
