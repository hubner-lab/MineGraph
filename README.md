
# MineGraph: Plastid and mitochondrial graph-pangenome analysis toolkit

**Overview**  
MineGraph is a computational pipeline designed to construct, analyze, and visualize graph-pangenomes of small genomes such as plastids and mitochondrial genomes. MineGraph include a full implementation of graph-pangenome construction and analysis and include automatic optimization process for key parfamters, graph construction, analysis, and visualization, thus offering a comprehensive solution for comparative genomics and variat calling. MineGraph was also optimized for speed and can handle hundreds of genomes efficiently with a single command. 
MineGraph is operated within a Docker contained environment, thus eliminating complex installation steps. 

### Key Terminology

- **Graph**: A data structure where nodes (vertices) represent entities, and edges define relationships between them.
- **Variation Graph**: A non-overlapping sequence graph where nodes represent sequence segments and edges show connections, useful in identifying genetic variations.
- **Pangenome Variation Graph**: A form of generic multiple sequence alignment. It reveals similarities where genomes traverse the same parts of the graph and differences where they do not.
- **PGGB (PanGenome Graph Builder)**: A tool used to build pangenome graphs from sequence data.
- **PanSN-spec**: A naming convention for labeling sequences in a pangenome, retaining this information across formats.
- **GFA (Graphical Fragment Assembly) Format**: A format used to represent sequence graphs.
  

### Key Functionalities and Advantages

1. **Automated Parameter Optimization**  
   MineGraph extracts and optimizes critical parameters, such as **mapping identity minimum (-p)** and **segment length (-s)**, which are essential for accurate graph construction. This automation ensures optimal graph structure, streamlining the process with minimal configuration.

   - **Maximum Divergence Calculation**: Using MASH Triangle, MineGraph calculates the ideal percentage identity, representing the divergence among input sequences.
   - **Handling Repetitive Regions**: Transposable elements and repetitive sequences are common in plant genomes, complicating pangenome analysis. MineGraph uses RepeatMasker to identify the longest tandem repeats (TRs) and sets this as the segment length for PGGB, improving graph accuracy in repetitive regions.

2. **Sequence Preparation and Compression**  
   MineGraph organizes, renames, and compresses sequences following the PanSN-spec convention. The compressed and indexed files serve as streamlined inputs for PGGB, enabling consistent and traceable sequence naming across formats.

3. **Comprehensive Graph Analysis and Outputs**  
   Once the pangenome graph is generated, MineGraph performs a robust statistical analysis and generates a variety of outputs, including:
   - **Graph Statistics**: Details like node and edge counts, and average node degree.
   - **Polymorphism Summary**: A VCF file summarizing Single Nucleotide Polymorphisms (SNPs), Multiple Nucleotide Polymorphisms (MNPs), and indels.
   - **Graph Visualization**: A circular plot where node size reflects sequence size, and edge width represents connection weight.
   - **Node Count Histogram**: A histogram representing node count distribution across the graph.
   - **Metrics Dataframes**: Detailed data including node IDs, counts, and sizes.
   - **Consensus Sequence**: The derived consensus sequence of the graph.
   - **Multiple Sequence Alignment (MSA)**: Alignment file reflecting the consensus across input sequences.
   - **Phylogenetic Tree**: A tree showing evolutionary relationships based on the pangenome graph.
   - **GFA Conversion**: Conversion of GFA format into FASTA and VG for compatibility with additional tools.
   - **SequenceTubeMap Integration:** Enables compact graph visualization, showing polymorphisms and homology via SequenceTubeMap built-in GUI.

4. **Docker Integration**  
   MineGraph is designed to run in a Docker container, simplifying setup by eliminating software dependency issues. It’s compatible with most systems, making it accessible to users without requiring extensive technical expertise.

### System Requirements

- **Docker installed**
- **Python installed with pandas package**

### Getting Started

To start using MineGraph, create a CSV or XLSX file listing the names of your FASTA files with a one column. Example:

```plaintext
Avena_sativa.fasta
Triticum.fasta
Zea_Mays.fasta
...
```

Place the FASTA files in a folder (e.g., `./my_data/`) in the current directory, then provide the folder as an argument when running MineGraph.

### Usage
**example 1 :** running the following command will run MineGraph pipeline and output all the definded workflows e.g pggb graph, phylogenetic tree, all the statistcs and plots. 

```bash
python MineGraph.py --data_dir <fasta_files_dir> --output_dir <out/dir>--metadata <csv_file>
```
**example 2 :** running the following command will run MineGraph pipeline and only outputs the TRs analysis results (optimized segment lenght parameter)
a params.yaml will be saved to the output dir with the input data prepared, hence named according PanSN-spec convention.

```bash
python MineGraph.py --data_dir input_fasta/ --output_dir results/ --metadata metadata.xlsx --mode extract-t
```

**example 3 :** running the following command will run MineGraph pipeline and only outputs PGGB graph after passing the optimized parameters + TRs output
a params.yaml will be saved to the output dir with the input data prepared and hence named according PanSN-spec convention.

```bash
python MineGraph.py --data_dir input_fasta/ --output_dir results/ --metadata metadata.xlsx --mode construct-graph
```

**example 4 :** running the following command will run MineGraph pipeline up to constructing the graph but plotting only the statistcal Datasheets (xlsx) of the graph.

```bash
python MineGraph.py --data_dir input_fasta/ --output_dir results/ --metadata metadata.xlsx --only-stats
```


if --view parameter is passed you should be able to access SequenceTubeMap api using http://localhost:3210/
the graph should be present on a custom track named "gfa_to_vg.vg.xg"
![](https://github.com/hubner-lab/MineGraph/blob/main/examples/stm.gif)
see : https://github.com/vgteam/sequenceTubeMap/blob/master/public/help/help.md

## Arguments

| Argument       | Required | Default   | Description                                                                                          |
|----------------|----------|-----------|------------------------------------------------------------------------------------------------------|
| `--data_dir`   | Yes      |           | Directory containing the input FASTA files to be processed.                                          |
| `--output_dir` | Yes      |           | Directory where the output files (e.g., MSA, trees, visualizations) will be saved.                   |
| `--metadata`   | Yes      |           | Metadata file (in CSV or XLSX format) listing the FASTA files to be processed.                       |
| `--threads`    | No       | 16        | Number of threads to use for parallel processing.                                                    |
| `--tree_pars`  | No       | 10        | Number of parsimonious trees to generate.                                                            |
| `--tree_bs`    | No       | 10        | Number of bootstrap trees to generate for assessing phylogenetic tree confidence.                    |
| `--quantile`   | No       | 50        | Consensus nodes percentage of presence. For example, `100` means nodes must appear in all paths.     |
| `--top_n`      | No       | 1000      | The top `N` node sizes to visualize in the output, sorted by size.                                   |
| `--window_size`| No       | 1000      | The sliding window size for similarity plot.                                                         |
| `--view`       | No       | False     | Runs SequenceTubeMap pipeline for visualizing the output graph.                                      |
| `--mode`       | No       | `all`     | Specifies which steps of the workflow to run. Options: `all` (full workflow), `extract-tr` (only extract tandem repeats), `construct-graph` (build the graph but skip statistics). |
| `--phyl-tree`  | No       | True      | Enables multiple sequence alignment (MSA) and phylogenetic tree generation.                         |
| `--plots`      | No       | True      | Enables plotting of statistical graphs during analysis.                                              |
| `--only-stats` | No       | False     | Runs only the statistical analysis step, generating XLSX files while skipping MSA, tree, and plots.  |

---



- **Examples:**
  ```bash
  python MineGraph.py --data_dir ./input_fasta/ --output_dir ./results/ --metadata metadata.csv --threads 32 --tree_pars 20 --tree_bs 20 --quantile 75 --top_n 500 --view
  python MineGraph.py --data_dir input_fasta/ --output_dir results/ --metadata metadata.xlsx --mode construct-graph
  python MineGraph.py --data_dir input_fasta/ --output_dir results/ --metadata metadata.xlsx --mode extract-tr
  python MineGraph.py --data_dir input_fasta/ --output_dir results/ --metadata metadata.xlsx --only-stats
  ```

The final output will be organized in the following folders:
- `/path/to/your/output/MineGraph_output`: Contains the pipeline’s final output files.
- `/path/to/your/output/pggb_output`: Includes the GFA file, VCF file, and an ODGI drawing of the final graph.

---

MineGraph provides a powerful, efficient approach to analyzing plant plastid and mitochondria genomes, equipping researchers with optimized pangenome graphs and comprehensive statistics for their genomic studies.
