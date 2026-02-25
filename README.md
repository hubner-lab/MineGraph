<p align="center">
  <img src="examples/logo.png" alt="MineGraph logo" width="260">
</p>

# MineGraph: Plastid and mitochondrial graph-pangenome analysis toolkit

---

## Overview

MineGraph is a Docker-based computational pipeline designed to construct, analyze, extract, convert, and visualize graph-pangenomes of plastid and mitochondrial genomes.

It integrates:

- Automatic PGGB parameter optimization  
- Graph construction (PGGB / ODGI)  
- Subgraph extraction  
- Graph format conversion (GFA ↔ FASTA / VG)  
- Graph statistics and path inspection  
- Interactive visualization via SequenceTubeMap  

MineGraph is optimized for speed and reproducibility and can handle hundreds of genomes with a single command.

---

## System Requirements

- Docker
- Python 3 (with pandas installed)

---

# 🔧 Usage Overview

MineGraph supports four main commands:

```
construct   → build full pangenome graph
extract     → subgraph extraction / stats / path inspection
convert     → GFA → FASTA or VG
tubemap     → interactive visualization
```

---

# 1️⃣ Construct a Pangenome Graph

```bash
python MineGraph.py construct --   --data_dir <fasta_files_dir>   --output_dir <out_dir>   --metadata <csv_or_xlsx>
```

This runs the full PGGB pipeline inside the `opggb` Docker container.

---

# 2️⃣ Extract Subgraphs & Inspect Graphs

The `extract` command wraps `/src/extract.py` inside the container.

## A) Subgraph Extraction

```bash
python MineGraph.py extract -- subgraph   -i graph.gfa   -w wanted.txt   -o subgraph.gfa
```

- `-i` input GFA
- `-w` wanted paths / anchors file
- `-o` output subgraph

Use `--mkdir` to auto-create missing output directories.

---

## B) Graph Statistics

```bash
python MineGraph.py extract -- stats   -i graph.gfa
```

Produces node/edge/path statistics for the graph.

---

## C) Inspect Graph Paths

```bash
python MineGraph.py extract -- paths   -i graph.gfa
```

Lists paths contained in the GFA graph.

You can forward additional extract.py arguments after `extract --`.

---

# 3️⃣ Convert Graph Formats

## GFA → FASTA

Linearizes all GFA `S` segments into one concatenated FASTA sequence.

```bash
python MineGraph.py convert fasta   -i graph.gfa   -o graph.fasta
```

---

## GFA → VG

Converts GFA to VG format using `vg convert` (stdout safely redirected internally).

```bash
python MineGraph.py convert vg   -i graph.gfa   -o graph.vg
```

Add `--mkdir` if output directory does not exist.

---

# 4️⃣ Interactive Visualization (SequenceTubeMap)

Launch SequenceTubeMap locally:

```bash
python MineGraph.py tubemap -i graph.gfa
```

Then open:

```
http://localhost:3210
```

The container is automatically stopped and removed when exiting or pressing `Ctrl+C`.

---

# 🎬 SequenceTubeMap Tutorial

<p align="center">
  <img src="examples/sequencetubemap_tutorial.gif" alt="SequenceTubeMap tutorial" width="900">
</p>

This tutorial demonstrates:

- Loading the graph
- Zooming and navigating
- Inspecting nodes and paths
- Exploring structural variation

Ensure the GIF file exists inside the `examples/` directory.

---

# 📂 Output Structure

## MineGraph_output

```
MineGraph_output/
├── statistics/
├── plots/
├── phylogenetics_msa/
├── gfa_convert/
```

## pggb_output

```
pggb_output/
├── *.final.gfa
├── *.final.vcf
├── *.final.og
├── *.maf
├── *.paf
├── *.params.yml
├── *.log
├── multiqc_report.html
```

---

# 🔎 Helpful Flags

```
--dry-run     → print docker command without running
--mkdir       → auto-create missing output directories
-m host:container  → add extra Docker mounts
```

---

MineGraph provides a powerful, efficient, and reproducible framework for plastid and mitochondrial graph-pangenome analysis.
