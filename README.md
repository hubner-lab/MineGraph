<p align="center">
  <img src="examples/logo.png" alt="MineGraph logo" width="420">
</p>

# MineGraph — Plastid & Mitochondrial Graph-Pangenome Toolkit

---

## Overview

MineGraph is a Docker-based pipeline for constructing, analysing, extracting, converting, and visualising graph-pangenomes of plastid and mitochondrial genomes.

**What it does:**

- Auto-tunes PGGB parameters (percent identity, segment length) from your data
- Constructs a variation graph with PGGB / ODGI
- Computes graph statistics (Cloud / Shell / Core node classification, π, density)
- Generates interactive and publication-quality plots
- Extracts subgraphs and inspects path composition
- Converts GFA graphs to FASTA or VG format
- Builds phylogenetic trees (graph-based Jaccard/Ward or RAxML MSA-based)
- Launches interactive visualisation via SequenceTubeMap

---

## System Requirements

| Requirement | Notes |
|-------------|-------|
| Docker      | Required to pull and run the `rakanhaib/opggb` image |
| Python ≥ 3.9 | With `pandas` installed on the **host** (`pip install pandas`) |
| ≥ 16 GB RAM | Recommended for typical large datasets |

---

## Architecture

```
Host machine
└── MineGraph.py          ← entry point (pure host-side Python)
    ├── construct  ──────── docker run rakanhaib/opggb → src/construct.py
    │                           └── calls sibling scripts directly (no docker-in-docker)
    ├── extract    ──────── docker run rakanhaib/opggb → src/extract.py
    ├── convert    ──────── docker run rakanhaib/opggb → vg / gfa2fasta
    └── tubemap    ──────── docker run rakanhaib/sequencetubemap
```

---

## Quick Start

```bash
# Full pipeline (most common usage)
python MineGraph.py construct -- \
  --data_dir  ./mito_Data \
  --output_dir ./out \
  --metadata  ./samples.csv \
  --threads   32 \
  --phyl-tree \
  --plots

# Just build the graph, skip stats
python MineGraph.py construct -- \
  --data_dir ./mito_Data --output_dir ./out --metadata ./samples.csv \
  --mode construct-graph

# Extract a subgraph
python MineGraph.py extract -- subgraph \
  -i out/pggb_output/graph.smooth.final.gfa \
  -w wanted_paths.txt \
  -o subgraph.gfa

# Convert to FASTA
python MineGraph.py convert fasta -i graph.gfa -o graph.fasta

# Visualise
python MineGraph.py tubemap -i graph.gfa
# then open http://localhost:3210
```

---

## Commands

### 1 · `construct` — Build a Pangenome Graph

```
python MineGraph.py construct -- [OPTIONS]
```

All options after `--` are forwarded to `src/construct.py` inside the container.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--data_dir` | path | **required** | Directory containing input FASTA files |
| `--output_dir` | path | **required** | Output root directory |
| `--metadata` | file | **required** | `.csv` or `.xlsx` listing FASTA filenames (first column) |
| `--threads` | int | `16` | CPU threads |
| `--mode` | choice | `all` | `all` · `extract-tr` · `construct-graph` |
| `--phyl-tree` | flag | off | Generate a phylogenetic tree |
| `--tree_type` | choice | `graph` | `graph` (Jaccard/Ward on ODGI matrix) or `msa` (RAxML) |
| `--tree_pars` | int | `10` | RAxML parsimonious starting trees |
| `--tree_bs` | int | `10` | RAxML bootstrap replicates |
| `--plots` | flag | off | Generate all statistical plots |
| `--only-stats` | flag | off | XLSX statistics only; skip plots and tree |
| `--quantile` | float | `25` | Consensus node presence threshold (%) |
| `--top_n` | int | `50` | Top-N nodes shown in interactive plot |
| `--window_size` | int | `1000` | Sliding window size for similarity analysis |
| `--sc_th` | int | `5` | Cloud/Core boundary (%). Nodes in ≤ threshold% of haplotypes → Cloud; ≥ (100−threshold)% → Core |
| `--sq_view` | flag | off | Launch SequenceTubeMap after construction |

**Workflow modes:**

| Mode | Steps run |
|------|-----------|
| `all` (default) | Prepare → RepeatMask → PGGB → Statistics |
| `extract-tr` | Prepare → RepeatMask only |
| `construct-graph` | Prepare → RepeatMask → PGGB only |

**Metadata file format** (CSV example):

```
sample_A.fasta
sample_B.fasta
sample_C.fa
```

One filename per row, first column. Additional columns are ignored.

---

### 2 · `extract` — Subgraph Extraction & Inspection

```
python MineGraph.py extract -- <subcommand> [OPTIONS]
```

#### A) Extract a subgraph

```bash
python MineGraph.py extract -- subgraph \
  -i graph.gfa  -w wanted.txt  -o subgraph.gfa
```

`wanted.txt` — one path name per line. The first listed path is used as the extraction anchor.

#### B) Graph statistics

```bash
python MineGraph.py extract -- stats -i graph.gfa -o stats.txt
```

#### C) List paths

```bash
python MineGraph.py extract -- paths -i graph.gfa -o paths.txt [--prefix Elymus_]
```

---

### 3 · `convert` — Format Conversion

#### GFA → FASTA

```bash
python MineGraph.py convert fasta -i graph.gfa -o graph.fasta
```

All segment sequences are concatenated into a single FASTA record labelled `<stem>#1`.

#### GFA → VG

```bash
python MineGraph.py convert vg -i graph.gfa -o graph.vg
```

Uses `vg convert` internally.

Add `--mkdir` to any convert command to auto-create missing output directories.

---

### 4 · `SequenceTubeMap` — Interactive Visualisation

```bash
python MineGraph.py tubemap -i graph.gfa
```

Then open **http://localhost:3210** in your browser. The container is automatically stopped when you press `Ctrl+C`.

<p align="center">
  <img src="examples/stm.gif" alt="SequenceTubeMap tutorial" width="900">
</p>

---

## Global Flags

| Flag | Description |
|------|-------------|
| `--dry-run` | Print the docker command without running it |
| `--mkdir` | Auto-create missing host-side output directories |
| `-m host:container` | Add extra Docker volume mounts (repeatable) |
| `--image IMAGE` | Override the opggb Docker image (default: `rakanhaib/opggb:latest`) |
| `--tubemap-image IMAGE` | Override the SequenceTubeMap image (default: `rakanhaib/sequencetubemap:latest`) |

---

## Output Structure

```
<output_dir>/
├── params.yaml                      # Auto-tuned PGGB parameters
├── panSN_output.fasta.gz            # Merged, renamed, compressed FASTA
├── downsampled_panSN_output.fasta   # Subsampled FASTA used for mash
├── pggb_output/
│   ├── *.smooth.final.gfa           # Final variation graph
│   ├── *.smooth.final.vcf           # Variant calls
│   ├── *.smooth.maf                 # Multiple alignment format
│   ├── *.alignments.wfmash.paf      # Pairwise alignments
│   ├── *.log                        # Per-tool logs
│   └── multiqc_report.html          # MultiQC summary
└── MineGraph_output/
    ├── statistics/
    │   ├── graph_stats_graph_stats.csv      # Node/edge/diversity summary
    │   └── graph_Node_Plot_frequency*       # Histogram CSV + PNG
    ├── plots/
    │   ├── graph_top_N_interactive.html     # Interactive PyVis graph
    │   ├── node_histogram_by_paths.png      # Cloud/Shell/Core bar chart
    │   ├── similarity.*                     # Window-based similarity heatmap
    │   └── paths.*                          # Path statistics heatmap
    ├── phylogenetics_msa/
    │   ├── Consensus_sequence.fasta
    │   ├── odgi_matrix.tsv
    │   ├── graph_phylo_tree.tree            # Newick (graph-based)
    │   ├── graph_phylo_tree.png
    │   └── MSA_result.fasta                 # (msa mode only)
    └── gfa_convert/
        ├── gfa_to_vg.vg
        └── gfa_to_fasta.fasta
```

---

## Node Classification

MineGraph classifies each graph node by the fraction of haplotypes that contain it, controlled by `--sc_th` (default 5 %):

| Category | Criterion | Meaning |
|----------|-----------|---------|
| **Cloud** | ≤ `sc_th` % of haplotypes | Private / rare regions |
| **Shell** | between Cloud and Core | Variably shared regions |
| **Core** | ≥ `(100 − sc_th)` % of haplotypes | Conserved backbone |

---

## Phylogenetic Tree Modes

| `--tree_type` | Method | When to use |
|---------------|--------|-------------|
| `graph` (default) | Jaccard distance on ODGI node-presence matrix → Ward clustering → Newick | Fast; no MSA required; works well for divergent datasets |
| `msa` | MAF → concatenated MSA → RAxML-NG (`GTR+G`) | Slower; use when sequence-level model fit matters |

---

## Citations

> Garrison E. et al. Building pangenome graphs. *Nature Methods* (2024). https://doi.org/10.1038/s41592-024-02430-3

> Beyer W. et al. SequenceTubeMap: visualization for graph-based genomes. *Bioinformatics* (2019). https://doi.org/10.1093/bioinformatics/btz597

> Vorbrugg S. et al. Gretl – Variation GRaph Evaluation TooLkit. *bioRxiv* (2024). https://doi.org/10.1101/2024.03.04.580974

---

MineGraph provides a fast, reproducible framework for plastid and mitochondrial graph-pangenome analysis.
