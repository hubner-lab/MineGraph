# Analysis Scripts

Downstream analysis scripts used to generate every figure and quantitative result reported in the accompanying manuscript. Each folder is named after the analysis it performs and can be run independently of the construction pipeline once the MineGraph outputs (`odgi_matrix.tsv`, `graph_phylo_tree.tree`, `graph_stats.xlsx`, etc.) are available.

## Folder index

| Folder | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| [`phylogenetic_comparison/`](phylogenetic_comparison/) | Compare the graph phylogeny against published trees (GPWG2, BK2015, Givnish, Saarela) | graph `.nwk`, published `.tre` | normalised Robinson-Foulds distances, clade-recovery bars, summary PDF |
| [`displacement_analysis/`](displacement_analysis/) | Per-taxon displacement metric and pairwise displacement matrices | two rooted trees in Newick | displacement TSVs, universal-taxa tables, displacement × bootstrap integration |
| [`site_concordance/`](site_concordance/) | Site concordance factors (SCF) on the MAF alignment | MAF, species tree | SCF histogram, SCF vs bootstrap scatter |
| [`pav_linkage_disequilibrium/`](pav_linkage_disequilibrium/) | Pairwise LD (`r²`) between PAV nodes and t-SNE sensitivity under LD pruning | `odgi_matrix.tsv` | LD decay, r² distribution, pruned-node retention, t-SNE grid |
| [`compartment_provenance/`](compartment_provenance/) | BLAST the MineGraph consensus against a chloroplast reference partitioned into LSC/IR/SSC | `Consensus_sequence.fasta`, reference compartments | per-family LSC/IR/SSC proportions, composite panel |
| [`gene_pav/`](gene_pav/) | Gene-level PAV heatmaps, family-stratified plots, annotation-free validation | injected graph, gene BED, `odgi_matrix.tsv` | gene-presence heatmap, top-variable gene heatmap, validation scatter |
| [`topology_analysis/`](topology_analysis/) | Graph-native topology metrics (path convergence, bubble density, symmetry, graph edit distance) | `.gfa` | topology TSVs, edit-distance matrix |
| [`go_enrichment/`](go_enrichment/) | topGO enrichment of consensus vs non-consensus node-mapped genes | gene2go, study/population IDs | enrichment TSVs, bubble plots |
| [`family_node_enrichment/`](family_node_enrichment/) | Per-family enrichment of family-unique graph nodes | `odgi_matrix.tsv`, metadata | family-specific node tables and plots |
| [`ward_clustering/`](ward_clustering/) | Ward-linkage clustering coloured by bootstrap and IR-content | `odgi_matrix.tsv` | Ward dendrogram with BS × IR overlay |
| [`pav_tsne/`](pav_tsne/) | t-SNE embedding of the PAV matrix coloured by Poaceae tribe | `odgi_matrix.tsv`, metadata | t-SNE PNG |
| [`figure_composition/`](figure_composition/) | Composite figure assembly (multi-panel layouts) | individual panels | Figure 2 composite PNG |

## Running the scripts

All scripts expect the MineGraph output layout produced by `MineGraph.py construct` (see the top-level [README](../README.md#output-structure)). Paths at the top of each script are written as `../graphs/192_graph/...` relative to the repository root — adjust to your output location.

- **Python (≥3.9):** `pandas`, `numpy`, `matplotlib`, `seaborn`, `scikit-learn`, `scipy`, `networkx`, `Bio` (for annotation-free gene PAV)
- **R (≥4.1):** `ape`, `phangorn`, `phytools`, `ggplot2`, `ggtree`, `treeio`, `topGO`, `ComplexHeatmap`, `Rtsne`, `vegan`

Install Python deps with:

```bash
pip install pandas numpy matplotlib seaborn scikit-learn scipy networkx biopython
```

Install R deps inside an R session:

```r
install.packages(c("ape","phangorn","phytools","ggplot2","Rtsne","vegan"))
BiocManager::install(c("ggtree","treeio","topGO","ComplexHeatmap"))
```

Individual folders that have non-trivial multi-script workflows (e.g. `phylogenetic_comparison/`) carry their own `README.md`.
