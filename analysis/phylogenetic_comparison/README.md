# Phylogenetic Comparison Analysis Pipeline
## Poales Pangenome Graph — Reviewer Validation Suite

**Purpose:** This directory contains all scripts that address the reviewer demand for "deep phylogenetic analysis" to validate the graph-derived tree. The original pipeline produced only a Ward-linkage PAV tree with no statistical support. These analyses supersede and extend that.

---

## Background: What the Reviewer Requested

> Reviewer statement: "The phylogenetic analyses presented rely on Ward clustering with no bootstrap support. Deep phylogenetic analysis is required to validate the graph tree's power relative to SNP-based phylogenies."

**Original pipeline limitation:** The 192-taxon graph-PAV tree (`192_graph_pav_tree_named.nwk`) was built by Ward-linkage hierarchical clustering on the PAV presence/absence matrix. Ward clustering produces a tree topology but provides **no statistical support metrics** — no bootstrap values, no concordance factors. A Genome Biology reviewer correctly flagged this as insufficient.

**Solution implemented:** Five complementary analyses (A1–A5) were added, plus extended comparisons (A6, Ward-vs-BS, 709-taxon extension).

---

## Scripts and Their Purpose

| Script | Analysis | Output files |
|--------|----------|-------------|
| `bootstrap_summary.R` | **A1** Bootstrap support on new BS tree | `bootstrap_support_summary.csv`, `FigSX_bootstrap_histogram.png` |
| `tree_distance_matrix.R` | **A2** Pairwise nRF, CID, Cophenetic r, Quartet distance (6×6 matrix) | `pairwise_tree_distances.csv` (+ individual CSVs) |
| `clade_recovery.R` | **A3** Monophyly of 20 established Poales clades per tree | `clade_recovery_table.csv`, `FigSX_clade_recovery_bars.png` |
| `run_scf.sh` | **A4+A5** IQ-TREE sCF on PAV binary alignment (MK+G4 model) | `pav_scf_minegraph_bs.cf.*`, `pav_scf_pav_ward.cf.*` |
| `scf_visualization.R` | **A4+A5** Visualize sCF distributions and BS×sCF scatter | `scf_summary_statistics.csv`, `scf_branch_table.csv`, `FigSX_scf_vs_bootstrap.png`, `FigSX_scf_histogram.png` |
| `pav_to_fasta.py` | **A5 prep** Convert odgi_matrix.tsv to binary FASTA (A=present, C=absent) | `pav_binary_variable.fasta` |
| `process_709_tree.py` | **709 tree** Strip `_1_0` panSN suffix, compute BS statistics | `709_graph_minegraph_bs_tree_named.nwk`, `709_bs_support_summary.csv` |
| `ward_bs_ir_analysis.R` | **Ward validation + IR** Topology agreement Ward vs BS; BS/sCF at family MRCAs vs IR size | `ward_vs_bs_comparison.csv`, `family_bs_ir_table.csv`, `FigSX_ward_bs_ir_analysis.png` |
| `compare_709_vs_published.R` | **Extended comparison** 709 tree + 192 trees vs GPWG2/Givnish/BK2015 | `709_vs_published_distances.csv`, `192_vs_published_extended.csv`, `FigSX_709_vs_published.png` |

---

## Key Results (as of 2026-04-14)

### A1 — Bootstrap Support (192-taxon MineGraph BS tree)
- Mean BS = **95.2%**, Median = 100%, SD = 10.3%
- 82.1% of nodes ≥ 90%; 64.2% = 100%
- Tree: `phylogenetic_trees/192_graph_minegraph_bs_tree_named.nwk`
- Built from: ODGI PAV matrix (`odgi_matrix.tsv`, 192 taxa × all nodes)
- Method: MineGraph `generate_graph_based_tree` — Jaccard distance (`scipy.spatial.distance.pdist`), UPGMA (`scipy.cluster.hierarchy.linkage(method="average")`), Felsenstein PAV column-resampling bootstrap, 1,000 replicates. NOT IQ-TREE; NOT consensus MSA.

### A2 — Pairwise Tree Distances (6×6 matrix)
- MineGraph_BS vs PAV_Ward: nRF = 0.28, r = 0.956 (very similar topology)
- MineGraph_BS vs RAxML: nRF = 0.44, r = −0.038
- MineGraph_BS vs GPWG2: nRF = 0.62, r = 0.667
- MineGraph_BS vs Givnish: nRF = 0.60, r = 0.078
- Files: `pairwise_tree_distances.csv`, `tree_distance_matrix_nRF.csv`, etc.

### A3 — Clade Recovery (20 clades × 6 trees)
- MineGraph_BS: Y=38.9%, Y+P+=100%
- RAxML: Y=44.4%, Y+P+=100%
- GPWG2: Y=27.8%, Y+P+=83.3%
- All major families present in graph trees
- File: `clade_recovery_table.csv`

### A4+A5 — Site Concordance Factor (PAV binary alignment, MK+G4)
- MineGraph_BS reference: **88.8% branches above 1/3 random threshold**, mean sCF = 49.0%
- PAV_Ward reference: 86.7% above 1/3 threshold, mean sCF = 48.3%
- r(BS, sCF) = 0.333 — BS and sCF are partially independent
- PAV alignment: 192 taxa × 285,728 variable binary sites (A=present, C=absent)
- IQ-TREE command: `--scf 1000 -m MK+G4 --seed 42`
- Files: `pav_scf_minegraph_bs.cf.stat`, `pav_scf_pav_ward.cf.stat`

### Ward vs BS Comparison
- nRF = **0.280** (moderate difference), cophenetic r = **0.956** (very high)
- 71.7% of BS-tree bipartitions are present in Ward tree
- **Interpretation:** Ward clustering captured the correct topology. The BS tree confirms and statistically validates the Ward groupings. The Ward tree was not "wrong" — it lacked support quantification.
- File: `ward_vs_bs_comparison.csv`

### Family BS vs IR Size
- ALL 7 families with recoverable MRCAs: BS = **100%** (saturated)
- sCF varies more informationally: Eriocaulaceae 98.5% → Typhaceae 56.6%
- sCF correlates inversely with IR size (Cyperaceae: largest IR=37.4kb, sCF=61%; Eriocaulaceae: IR=26.4kb, sCF=98.5%)
- Cyperaceae is paraphyletic in PAV/Ward tree (large IR may confound PAV signal)
- File: `family_bs_ir_table.csv`

### 709-Taxon Tree
- Mean BS = **94.9%**, Median = 100%, SD = 13.4%
- 87.0% of nodes ≥ 90%; 71.9% = 100%
- vs BK2015: nRF=0.553, cophenetic r=**0.813**, 45% shared bipartitions
- vs Givnish2018: nRF=0.633, r=0.694
- vs GPWG2: nRF=0.711, r=0.639
- File: `709_graph_minegraph_bs_tree_named.nwk`

---

## Method Notes

### Why Felsenstein Bootstrap (not Bayesian posterior)?
Felsenstein (1985) bootstrap resamples alignment columns (here: graph nodes for BS tree, binary PAV sites for sCF). It is widely accepted, computationally tractable for large datasets, and directly comparable across trees. Bayesian methods were not used because we do not have a prior on graph structural variation.

### Why MK+G4 model for PAV sCF?
The Mk (Lewis 2001) model is the standard for binary/morphological characters. It assumes equal rates for 0→1 and 1→0 transitions. +G4 adds discrete Gamma rate variation (4 categories) to accommodate heterogeneous node evolutionary rates. This is the only appropriate model for a binary presence/absence matrix; sequence models (GTR, HKY) are not applicable.

### Why not sCF on the consensus MSA?
The MineGraph-BS tree is built from the ODGI PAV matrix (Jaccard UPGMA), not from the consensus MSA. The consensus sequence (~16.9 kb) is a MineGraph output artefact — it represents only ~11.5% of the chloroplast and is NOT used for tree construction in the current pipeline. For sCF, we chose the **PAV binary alignment** (285,728 sites) because:
1. It is independent of the consensus sequence (tests PAV signal, not sequence signal)
2. It has far more sites (285k vs 16.9k), giving better sCF precision
3. It tests whether structural variation (presence/absence of graph nodes) corroborates the tree topology independently of sequence similarity

### IQ-TREE sCF vs BS: what do they measure?
- **Bootstrap (BS):** resampling stability of tree topology across random subsets of sites. High BS = the same bipartition appears in many resampled trees.
- **sCF:** fraction of variable sites that unambiguously support a given branch over its two alternatives. sCF > 33.3% = above random chance. These measure different aspects of support: BS tests sampling stability; sCF tests concordance.

---

## File Locations

All output files are in `phylogenetic_trees/` or `manuscript/figures/`.

**Tree files:**
- `192_graph_minegraph_bs_tree_named.nwk` — primary 192-taxon BS tree (MineGraph Felsenstein 1000 BS)
- `192_graph_pav_tree_named.nwk` — original Ward PAV tree (no BS)
- `192_graph_RAxML_tree_named.nwk` — RAxML gene-based tree
- `709_graph_minegraph_bs_tree_named.nwk` — 709-taxon full dataset BS tree

**Published reference trees:**
- `GPWG2_2012_poaceae_renamed.nwk` — 90/213 tips
- `Givnish2018_monocot_renamed.nwk` — 63/208 tips (monocot plastome)
- `BK2015_MCC_POALES_renamed.nwk` — 88/223 tips (Poales time-calibrated)

---

## Execution Order (if re-running from scratch)

```bash
# 1. Bootstrap summary (Analysis A1)
Rscript scripts/phylo_comparison/bootstrap_summary.R

# 2. PAV binary FASTA (Analysis A5 prep)
python3 scripts/phylo_comparison/pav_to_fasta.py

# 3. IQ-TREE sCF runs (A4+A5) — takes ~20 min
bash scripts/phylo_comparison/run_scf.sh

# 4. sCF visualization (A4+A5)
Rscript scripts/phylo_comparison/scf_visualization.R

# 5. Tree distance matrix (A2)
Rscript scripts/phylo_comparison/tree_distance_matrix.R

# 6. Clade recovery (A3)
Rscript scripts/phylo_comparison/clade_recovery.R

# 7. Process 709 tree
python3 scripts/phylo_comparison/process_709_tree.py

# 8. Ward vs BS + IR analysis
Rscript scripts/phylo_comparison/ward_bs_ir_analysis.R

# 9. 709 vs published comparison
Rscript scripts/phylo_comparison/compare_709_vs_published.R
```

---

## Connection to IR Analysis

See `rebuttal/IR_GRAPH_ANALYSIS.md` for the full IR analysis. Key connection points:
- 63% of core graph nodes are estimated to be IR-derived
- Cyperaceae has the largest IR (37.4 kb), Xyridaceae the smallest (17.4 kb)
- Family sCF values (PAV) are lower for large-IR families: Cyperaceae sCF=61%, vs Eriocaulaceae sCF=98.5%
- Bootstrap support is uniformly 100% for all major families (BS saturated) — sCF is the more discriminating metric
- The Ward PAV tree shows Cyperaceae as paraphyletic — consistent with large IR creating confounding PAV signal when IR boundary differences are misinterpreted as genuine presence/absence variation

---

*Last updated: 2026-04-14 — analyses A1–A5, Ward-vs-BS, 709-taxon extension, IR-bootstrap analysis*
