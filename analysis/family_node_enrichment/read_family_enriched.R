#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Extract per-family node sequences from a PGGB GFA into FASTA
# with FULL support for oriented node IDs like: 123+ / 123- / 123
#
# New features requested:
#  1) --min_len CLI argument (default 20)
#  2) Plot node length distributions per family (PNG)
#  3) Keep short nodes (< min_len) but write them to a separate FASTA
#
# Inputs expected (from your enrichment outputs):
#   results_dir/<Family>/family_nodes_{enriched|unique|core_within_family}.tsv
# TSV must contain a column named: node
# Nodes can be: 70867+, 70867-, 70867, node.70867, node.70867+
#
# Outputs:
#   out_dir/<Family>/family_nodes_<which>_minlen<min_len>.fa          (kept >= min_len)
#   out_dir/<Family>/family_nodes_<which>_short_lt<min_len>.fa       (short < min_len)  [optional]
#   out_dir/node_lengths_<which>.tsv                                 (all nodes lengths)
#   out_dir/node_length_distributions_<which>.png                    (plot)
#   out_dir/ALL_FAMILIES_missing_nodes_<which>.tsv                   (missing)
#
# Usage:
#   Rscript read_Family_enriched.r \
#     --gfa 192_graph_gfa.gfa \
#     --results_dir results_families \
#     --which enriched \
#     --out_dir family_fastas_enriched \
#     --min_len 20 \
#     --plot TRUE \
#     --write_short TRUE
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(Biostrings)   # reverseComplement(DNAString())
  library(ggplot2)
  library(fs)
})

# --------------------------- CLI ----------------------------

opt_list <- list(
  make_option(c("--gfa"), type="character", help="Input PGGB graph in GFA1 (.gfa or .gfa.gz)"),
  make_option(c("--results_dir"), type="character", default="results_families",
              help="Root folder containing <Family>/family_nodes_*.tsv [default %default]"),
  make_option(c("--which"), type="character", default="enriched",
              help="Which node set per family: enriched | unique | core_within_family [default %default]"),
  make_option(c("--out_dir"), type="character", default="family_node_fastas",
              help="Output directory [default %default]"),
  make_option(c("--min_len"), type="integer", default=20,
              help="Minimum node length (bp) to keep in main FASTA [default %default]"),
  make_option(c("--write_short"), type="logical", default=TRUE,
              help="If TRUE, write short nodes (<min_len) to a separate FASTA per family [default %default]"),
  make_option(c("--one_fasta_per_family"), type="logical", default=TRUE,
              help="Write one FASTA per family (recommended) [default %default]"),
  make_option(c("--also_write_per_node"), type="logical", default=FALSE,
              help="Also write one FASTA per oriented node in each family folder [default %default]"),
  make_option(c("--wrap"), type="integer", default=80,
              help="FASTA line wrap width [default %default]"),
  make_option(c("--chunk_lines"), type="integer", default=500000,
              help="How many GFA lines per chunk [default %default]"),
  make_option(c("--plot"), type="logical", default=TRUE,
              help="If TRUE, make a per-family length distribution plot [default %default]"),
  make_option(c("--plot_max_families"), type="integer", default=60,
              help="Max families to include in a faceted plot (avoid huge figures) [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

allowed <- c("enriched", "unique", "core_within_family")
if (is.null(opt$gfa)) stop("Missing --gfa")
if (!file.exists(opt$gfa)) stop("GFA file not found: ", opt$gfa)
if (!fs::dir_exists(opt$results_dir)) stop("results_dir not found: ", opt$results_dir)
if (!(opt$which %in% allowed)) stop("--which must be one of: ", paste(allowed, collapse=", "))
if (is.na(opt$min_len) || opt$min_len < 1) stop("--min_len must be >= 1")

fs::dir_create(opt$out_dir)

# -------------------- Helpers -------------------------------

open_gfa_con <- function(gfa_path) {
  if (grepl("\\.gz$", gfa_path, ignore.case = TRUE)) gzfile(gfa_path, open = "rt") else file(gfa_path, open = "rt")
}

# Parses: "70867+", "70867-", "70867", "node.70867", "node.70867-"
# Returns list(segment="70867", strand="+"/"-")
parse_oriented_node <- function(x, default_strand = "+") {
  x <- trimws(as.character(x))
  x <- sub("^node\\.", "", x)
  
  m <- regexec("^([^\\s\\t]+?)([+-])?$", x)
  r <- regmatches(x, m)[[1]]
  if (length(r) == 0) stop("Could not parse node id: '", x, "'")
  
  seg <- r[2]
  strand <- r[3]
  if (is.na(strand) || strand == "") strand <- default_strand
  
  list(segment = seg, strand = strand)
}

wrap_seq <- function(seq, width = 80L) {
  if (is.na(seq) || seq == "*" || seq == "") return(character(0))
  starts <- seq(1, nchar(seq), by = width)
  vapply(starts, function(i) substr(seq, i, min(i + width - 1L, nchar(seq))), character(1))
}

write_fasta <- function(headers, seqs, out_fa, width = 80L) {
  stopifnot(length(headers) == length(seqs))
  con <- file(out_fa, open = "wt")
  on.exit(close(con), add = TRUE)
  
  for (i in seq_along(headers)) {
    writeLines(paste0(">", headers[[i]]), con)
    s <- seqs[[i]]
    if (!is.na(s) && s != "*" && s != "") {
      writeLines(wrap_seq(s, width), con)
    }
  }
}

apply_strand <- function(seq, strand) {
  if (is.na(seq) || seq == "*" || seq == "") return(seq)
  if (strand == "+") return(seq)
  if (strand == "-") return(as.character(reverseComplement(DNAString(seq))))
  stop("Unexpected strand: ", strand)
}

pick_tsv <- function(fam_dir, which_key) {
  fname <- switch(which_key,
                  enriched = "family_nodes_enriched.tsv",
                  unique   = "family_nodes_unique.tsv",
                  core_within_family = "family_nodes_core_within_family.tsv")
  fpath <- fs::path(fam_dir, fname)
  if (!file.exists(fpath)) return(NA_character_)
  fpath
}

# -------------------- 1) Collect per-family nodes ------------

cat("Reading per-family node lists from:", opt$results_dir, "\n")

family_dirs <- fs::dir_ls(opt$results_dir, type = "directory", recurse = FALSE)
if (length(family_dirs) == 0) stop("No family directories found in: ", opt$results_dir)

family_nodes <- list()
needed_segments_env <- new.env(parent = emptyenv(), hash = TRUE)

for (fam_dir in family_dirs) {
  fam <- fs::path_file(fam_dir)
  tsv <- pick_tsv(fam_dir, opt$which)
  if (is.na(tsv)) next
  
  dt <- fread(tsv)
  if (!("node" %in% names(dt))) {
    warning("Skipping ", tsv, " (no 'node' column)")
    next
  }
  
  nodes <- unique(dt$node)
  parsed <- lapply(nodes, parse_oriented_node)
  
  segs <- unique(vapply(parsed, `[[`, character(1), "segment"))
  family_nodes[[fam]] <- list(nodes = nodes, parsed = parsed, segs = segs, tsv = tsv)
  
  for (s in segs) assign(s, TRUE, envir = needed_segments_env)
}

if (length(family_nodes) == 0) {
  stop("No usable family node TSVs found for --which=", opt$which, " under: ", opt$results_dir)
}

needed_segments <- ls(needed_segments_env)
cat("Families loaded:", length(family_nodes), "\n")
cat("Total unique segments needed:", length(needed_segments), "\n")

# -------------------- 2) Stream-parse the GFA -----------------
# Build segName -> sequence only for needed segments.

seq_map <- new.env(parent = emptyenv(), hash = TRUE)

cat("Parsing GFA (streaming):", opt$gfa, "\n")
con <- open_gfa_con(opt$gfa)
on.exit(close(con), add = TRUE)

found <- 0L
chunk_n <- opt$chunk_lines

repeat {
  lines <- readLines(con, n = chunk_n)
  if (length(lines) == 0) break
  
  s_lines <- lines[startsWith(lines, "S\t")]
  if (length(s_lines) == 0) next
  
  parts <- strsplit(s_lines, "\t", fixed = TRUE)
  
  for (p in parts) {
    if (length(p) < 3) next
    seg <- p[[2]]
    
    if (exists(seg, envir = needed_segments_env, inherits = FALSE) &&
        !exists(seg, envir = seq_map, inherits = FALSE)) {
      assign(seg, p[[3]], envir = seq_map)
      found <- found + 1L
    }
  }
  
  if (found >= length(needed_segments)) break
}

cat("Found sequences for", found, "/", length(needed_segments), "needed segments\n")

# -------------------- 3) Write FASTA + collect lengths --------

cat("Writing FASTA to:", opt$out_dir, "\n")

missing_report <- list()
len_rows <- list()
len_i <- 0L

for (fam in names(family_nodes)) {
  nodes <- family_nodes[[fam]]$nodes
  
  fam_out_dir <- fs::path(opt$out_dir, fam)
  fs::dir_create(fam_out_dir)
  
  kept_headers <- character(0)
  kept_seqs    <- character(0)
  
  short_headers <- character(0)
  short_seqs    <- character(0)
  
  missing <- character(0)
  filtered_short <- 0L
  kept_n <- 0L
  
  for (node in nodes) {
    pr <- parse_oriented_node(node)
    seg <- pr$segment
    strand <- pr$strand
    
    if (!exists(seg, envir = seq_map, inherits = FALSE)) {
      missing <- c(missing, node)
      
      len_i <- len_i + 1L
      len_rows[[len_i]] <- data.table(
        Family = fam, node = node, segment = seg, strand = strand,
        length = NA_integer_, bucket = "missing"
      )
      next
    }
    
    base_seq <- get(seg, envir = seq_map, inherits = FALSE)
    seq <- apply_strand(base_seq, strand)
    L <- nchar(seq)
    
    bucket <- if (L >= opt$min_len) "kept" else "short"
    
    len_i <- len_i + 1L
    len_rows[[len_i]] <- data.table(
      Family = fam, node = node, segment = seg, strand = strand,
      length = as.integer(L), bucket = bucket
    )
    
    hdr <- paste0(
      node,
      "|seg=", seg,
      "|strand=", strand,
      "|len=", L,
      "|family=", fam,
      "|set=", opt$which
    )
    
    if (L >= opt$min_len) {
      kept_headers <- c(kept_headers, hdr)
      kept_seqs    <- c(kept_seqs, seq)
      kept_n <- kept_n + 1L
      
      if (isTRUE(opt$also_write_per_node)) {
        out_one <- fs::path(fam_out_dir, paste0(node, ".fa"))
        write_fasta(hdr, seq, out_one, width = opt$wrap)
      }
    } else {
      filtered_short <- filtered_short + 1L
      if (isTRUE(opt$write_short)) {
        short_headers <- c(short_headers, hdr)
        short_seqs    <- c(short_seqs, seq)
      }
    }
  }
  
  if (isTRUE(opt$one_fasta_per_family)) {
    out_kept <- fs::path(fam_out_dir, paste0("family_nodes_", opt$which, "_minlen", opt$min_len, ".fa"))
    write_fasta(kept_headers, kept_seqs, out_kept, width = opt$wrap)
    
    if (isTRUE(opt$write_short)) {
      out_short <- fs::path(fam_out_dir, paste0("family_nodes_", opt$which, "_short_lt", opt$min_len, ".fa"))
      write_fasta(short_headers, short_seqs, out_short, width = opt$wrap)
    }
  }
  
  if (length(missing) > 0) {
    miss_tsv <- fs::path(fam_out_dir, paste0("missing_nodes_", opt$which, ".tsv"))
    fwrite(data.table(node = missing), miss_tsv, sep = "\t")
    missing_report[[fam]] <- length(missing)
  } else {
    missing_report[[fam]] <- 0L
  }
  
  cat(
    "Family:", fam,
    "kept(>=", opt$min_len, "bp):", kept_n,
    "short(<", opt$min_len, "bp):", filtered_short,
    "missing:", length(missing), "\n"
  )
}

# -------------------- 4) Write summary tables -----------------

len_dt <- rbindlist(len_rows, use.names = TRUE, fill = TRUE)
len_tsv <- fs::path(opt$out_dir, paste0("node_lengths_", opt$which, ".tsv"))
fwrite(len_dt, len_tsv, sep = "\t")

miss_all <- data.table(
  Family = names(missing_report),
  missing_nodes = as.integer(unlist(missing_report))
)
miss_tsv <- fs::path(opt$out_dir, paste0("ALL_FAMILIES_missing_nodes_", opt$which, ".tsv"))
fwrite(miss_all, miss_tsv, sep = "\t")

# -------------------- 5) Plot distributions per family --------

if (isTRUE(opt$plot)) {
  
  # ---- Keep only valid measured nodes (kept/short) ----
  plot_dt_all <- len_dt[
    bucket %in% c("kept", "short") &
      !is.na(length) &
      length > 0
  ]
  
  # ---- Count the 1-bp spike (we will exclude it from the histogram) ----
  ones_dt <- plot_dt_all[length == 1, .N, by = Family][order(-N)]
  
  # ---- Exclude length == 1 from the plotted distribution ----
  plot_dt <- plot_dt_all[length > 1]
  
  if (nrow(plot_dt) == 0) {
    warning("No nodes with length > 1 available for plotting after filtering.")
  } else {
    
    # If too many families, keep the top families by number of nodes (kept+short)
    fam_counts <- plot_dt[, .N, by = Family][order(-N)]
    keep_fams <- head(fam_counts$Family, opt$plot_max_families)
    plot_dt2 <- plot_dt[Family %in% keep_fams]
    
    # Ensure a stable order in facets
    plot_dt2[, Family := factor(Family, levels = keep_fams)]
    
    # ---- Log-spaced bins (critical fix: avoids huge first bin on linear bins) ----
    log_breaks <- c(2, 5, 10, 20, 50, 100, 200, 500, 1e3, 2e3, 5e3, 1e4, 2e4)
    
    p <- ggplot(plot_dt2, aes(x = length, fill = bucket)) +
      geom_histogram(
        breaks = log_breaks,
        position = "identity",
        alpha = 0.7,
        color = "black",
        linewidth = 0.1
      ) +
      scale_x_log10(
        breaks = log_breaks,
        labels = scales::comma
      ) +
      scale_fill_manual(
        name = "Node length",
        values = c(kept = "#F8766D", short = "#00BFC4"),
        labels = c(
          kept  = "\u2265 50 bp",
          short = "< 50 bp"
        )
      ) +
      facet_wrap(~Family, scales = "free_y") +
      labs(
        title = paste0(
          "Node length distributions per family\n",
          "set=", opt$which,
          ", kept min_len=", opt$min_len, " bp; length=1 bp excluded from histogram"
        ),
        x = "Node length (bp)",
        y = "Count"
      ) +
      theme_bw() +
      theme(
        legend.position = "top",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 13),
        
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.text.x  = element_text(size = 13, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 13),
        
        strip.text   = element_text(size = 14, face = "bold"),
        plot.title   = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    
    
    plot_file <- fs::path(opt$out_dir, paste0("node_length_distributions_", opt$which, ".png"))
    ggsave(plot_file, p, width = 18, height = 10, dpi = 300)
    cat("Wrote plot:", plot_file, "\n")
    
    if (nrow(fam_counts) > opt$plot_max_families) {
      cat(
        "Note: plot includes only top", opt$plot_max_families,
        "families by node count (to keep figure manageable).\n"
      )
    }
  }
  
  # ---- Report (do not hide) how many 1-bp nodes were excluded ----
  if (nrow(ones_dt) > 0) {
    cat("\nNodes with length = 1 bp (excluded from histogram, still in TSV):\n")
    print(ones_dt)
  } else {
    cat("\nNo nodes with length = 1 bp detected.\n")
  }
}

cat("\nDONE.\n")
cat("Lengths TSV:  ", len_tsv, "\n")
cat("Missing TSV:  ", miss_tsv, "\n")
cat("Output dir:   ", opt$out_dir, "\n")
