
<p align="center">
  <img src="examples/logo.png" alt="MineGraph logo" width="260">
</p>

# MineGraph: Plastid and mitochondrial graph-pangenome analysis toolkit

## Overview

MineGraph is a computational pipeline designed to construct, analyze, and visualize graph-pangenomes of small genomes such as plastids and mitochondrial genomes. MineGraph includes a full implementation of graph-pangenome construction and analysis and incorporates an automatic optimization process for key parameters, graph construction, analysis, and visualization, thus offering a comprehensive solution for comparative genomics and variant calling. MineGraph was optimized for speed and can handle hundreds of genomes efficiently with a single command.

MineGraph is operated within Docker containerized environments, eliminating complex installation steps and dependency management.

---

## Key Terminology

- **Graph**: A data structure where nodes represent sequence segments and edges define adjacency relationships.
- **Variation Graph**: A non-overlapping sequence graph where nodes represent sequence segments and edges show connections, useful for identifying genetic variation.
- **Pangenome Variation Graph**: A generalized multiple sequence alignment revealing conserved and divergent genomic regions.
- **PGGB (PanGenome Graph Builder)**: A toolchain used to construct pangenome graphs from sequence data.
- **PanSN-spec**: A standardized naming convention for labeling sequences in a pangenome.
- **GFA (Graphical Fragment Assembly)**: A file format used to represent sequence graphs.

---

## Key Functionalities and Advantages

### 1. Automated Parameter Optimization

MineGraph extracts and optimizes critical parameters required for accurate graph construction, including **mapping identity minimum (-p)** and **segment length (-s)**.

- **Maximum Divergence Calculation**: Using MASH triangle distance, MineGraph estimates the optimal percentage identity threshold, representing divergence among input sequences.
- **Handling Repetitive Regions**: Repetitive elements are common in plant organellar genomes. MineGraph uses RepeatMasker to identify the longest tandem repeat (TR) and sets this as the PGGB segment length, improving graph resolution in repetitive regions.

---

### 2. Sequence Preparation and Compression

MineGraph organizes, renames, and compresses sequences following the **PanSN-spec** convention. The resulting indexed inputs ensure consistent sequence naming across FASTA, GFA, VG, and downstream analyses.

---

### 3. Comprehensive Graph Analysis and Outputs

After graph construction, MineGraph performs extensive downstream analyses and generates multiple outputs, including:

- **Graph statistics** (node/edge counts, average node degree)
- **Polymorphism summary** (VCF of SNPs, MNPs, and indels)
- **Graph visualization** (circular plots with node size proportional to sequence length)
- **Node count histograms**
- **Metrics dataframes** (node IDs, counts, sizes)
- **Consensus sequence**
- **Multiple Sequence Alignment (MSA)**
- **Phylogenetic tree**
- **Graph format conversion** (GFA → FASTA, GFA → VG)
- **Interactive visualization via SequenceTubeMap**

---

### 4. Docker Integration

MineGraph runs fully within Docker containers, ensuring reproducibility and portability across systems without requiring complex installations.

---

## System Requirements

- **Docker**
- **Python (with pandas installed)**

---

## Getting Started

Prepare a CSV or XLSX metadata file containing **one column listing FASTA filenames**:

```
Avena_sativa.fasta
Triticum.fasta
Zea_Mays.fasta
```

Place the FASTA files in a directory (e.g., `./my_data/`) and provide this directory as input when running MineGraph.

---

## Usage

### Full MineGraph pipeline

```bash
python MineGraph.py construct -- \
  --data_dir <fasta_files_dir> \
  --output_dir <out_dir> \
  --metadata <csv_or_xlsx>
```

---

### Graph format conversion

#### GFA → FASTA

Produces a linearized FASTA sequence by concatenating all GFA segments.

```bash
python MineGraph.py convert fasta \
  -i results/pggb_output/example.final.gfa \
  -o results/gfa_convert/example.fasta
```

#### GFA → VG

Converts a GFA graph into VG format (stdout redirected internally).

```bash
python MineGraph.py convert vg \
  -i results/pggb_output/example.final.gfa \
  -o results/gfa_convert/example.vg
```

Use `--mkdir` to auto-create output directories if needed.

---

## Interactive Graph Visualization (SequenceTubeMap)

```bash
python MineGraph.py tubemap -i results/pggb_output/example.final.gfa
```

Open:

```
http://localhost:3210
```

The SequenceTubeMap container is safely stopped and removed on exit or `Ctrl+C`.

---

## Output Structure

### MineGraph_output

```
MineGraph_output/
├── statistics/          # Graph-level and node-level statistics (XLSX, TSV)
├── plots/               # Circular plots, histograms, summary figures
├── phylogenetics_msa/   # MSA files and phylogenetic trees
├── gfa_convert/         # GFA → FASTA / VG conversions
```

### pggb_output

The `pggb_output` directory contains raw and processed outputs generated by PGGB and ODGI.
File names are parameterized and timestamped; below we list **only file types / suffixes**.

```
pggb_output/
├── *.final.gfa
├── *.final.vcf
├── *.final.vcf.stats
├── *.final.og
├── *.final.og.lay
├── *.final.og.lay.draw.png
├── *.final.og.lay.tsv
├── *.final.og.stats.yaml
├── *.final.og.viz_*.png
├── *.maf
├── *.paf
├── *.params.yml
├── *.log
├── *.fix.affixes.tsv.gz
├── *.seqwish.og.stats.yaml
├── multiqc_report.html
├── multiqc_data/
├── multiqc_config.yaml
```

---

MineGraph provides a powerful, efficient, and reproducible framework for plastid and mitochondrial graph-pangenome analysis, delivering optimized pangenome graphs and comprehensive statistics for comparative genomics.
