import argparse
import os
import subprocess

from utils import (
    get_nodes_edges_paths, get_graph_stats, get_nodes_counts_histogram, maf_to_msa,
    draw_distribution, plot_node_histogram_by_paths, extract_graph_stats,
    extract_consensus_path, draw_graph_as_circle_interactive, compute_min_node_count,
    write_consensus_sequence, draw_consensus_path, plot_graph_statistics,
    render_phylotree, write_odgi_matrix
)





def analyze_graph_and_alignments(threads="8", tree_pars=10, tree_bs=10, input_dir="/data/pggb_output",
                                 output_dir="/data/MineGraph_output", quantile=0.25, top_n=50,
                                 window_size=1000, phyl_tree=True, plots=True, only_stats=False):
    os.makedirs(output_dir, exist_ok=True)
    stats_dir = os.path.join(output_dir, "MineGraph_output", "statistics")
    graph_dir = os.path.join(output_dir, "MineGraph_output", "plots")
    msa_dir = os.path.join(output_dir, "MineGraph_output", "phylogenetics_msa")
    os.makedirs(stats_dir, exist_ok=True)
    os.makedirs(graph_dir, exist_ok=True)
    os.makedirs(msa_dir, exist_ok=True)

    # Locate required files
    gfa_file, paf_file, vcf_file, maf_file = None, None, None, None
    for file in os.listdir(input_dir):
        if file.endswith("smooth.final.gfa"):
            gfa_file = os.path.join(input_dir, file)
        elif file.endswith("alignments.wfmash.paf"):
            paf_file = os.path.join(input_dir, file)
        elif file.endswith(".vcf.stats"):
            vcf_file = os.path.join(input_dir, file)
        elif file.endswith(".smooth.maf"):
            maf_file = os.path.join(input_dir, file)

    if not gfa_file or not paf_file or not maf_file:
        print("Error: Required GFA, PAF, or MAF files not found.")
        return

    # Extract graph statistics
    nodes, edges, paths = get_nodes_edges_paths(gfa_file)
    nodes_num, edges_num, avg_degree, num_paths, edges_count = get_graph_stats(nodes, edges, paths)
    node_count_hist = get_nodes_counts_histogram(paths)

    # Assign count_dist directly without conditionals
    count_dist = draw_distribution(node_count_hist, os.path.join(stats_dir, "graph_Node_Plot_frequency"))

    if not only_stats:
        print("Generating graphical representations...")
        draw_graph_as_circle_interactive(nodes, edges_count, node_count_hist,
                                         os.path.join(graph_dir, "graph_top_N_interactive"), top_n)

    # Save statistics
    extract_graph_stats(node_count_hist, edges_count, nodes_num, edges_num, avg_degree, num_paths,
                        os.path.join(stats_dir, "graph_stats"), vcf_file, count_dist)

    # Skip MSA and Phylogenetic tree if only_stats is enabled
    if only_stats:
        print("Skipping MSA and Phylogenetic tree generation (only stats mode enabled).")
        return

    if plots:
        plot_graph_statistics(gfa_file, graph_dir, window_size)
        plot_node_histogram_by_paths(count_dist, num_paths, os.path.join(graph_dir, "node_histogram_by_paths.png"))

    # Step 6: Write Consensus path sequence to file
    cons_path, cons_seq = extract_consensus_path(nodes, edges, node_count_hist,
                                                 edges_count,
                                                 min_node_count=compute_min_node_count(node_count_hist,
                                                                                       num_paths,quantile=quantile))
    write_consensus_sequence(msa_dir, cons_seq)
    draw_consensus_path(nodes, edges_count, cons_path, node_count_hist,
                            os.path.join(msa_dir, "graph_consensus_interactive"), top_n=top_n)
    write_odgi_matrix(gfa_file, msa_dir)
    # Generate MSA and Phylogenetic tree if enabled
    if phyl_tree:
        print("Generating MSA and Phylogenetic Tree...")
        output_fasta = os.path.join(msa_dir, "MSA_result.fasta")
        maf_to_msa(maf_file, output_fasta)

        # Generate Phylogenetic Tree using RAxML
        raxml_output = os.path.join(output_dir, "graph.raxml.support")
        if not os.path.exists(raxml_output):
            print(f"Generating phylogenetic tree with RAxML ({threads} threads, search1 mode)...")
            raxml_command = [
                "raxml-ng", "--all", "--msa", output_fasta,
                "--model", "GTR+G",
                "--tree", f"pars{{{tree_pars}}}",
                "--bs-trees", f"{tree_bs}",
                "--threads", f"{threads}",
                "--prefix", os.path.join(output_dir, "graph")
            ]
            subprocess.run(raxml_command, check=True)

        # Render the tree
        if os.path.exists(raxml_output):
            render_phylotree(raxml_output, os.path.join(msa_dir, "graph_phylo_tree"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Graph and Alignment Analysis")
    parser.add_argument("--threads", default="8", help="Number of threads to use")
    parser.add_argument("--tree_pars", type=int, default=10, help="Number of parsimonious trees")
    parser.add_argument("--tree_bs", type=int, default=10, help="Number of bootstrap trees")
    parser.add_argument("--input_dir", default="/data/pggb_output", help="Input directory")
    parser.add_argument("--output_dir", default="/data/MineGraph_output", help="Output directory")
    parser.add_argument("--quantile", type=float, default=0.25, help="Consensus nodes top percentage (default: 25%)")
    parser.add_argument("--top_n", type=int, default="50", help="Top N nodes sizes to visualize")
    parser.add_argument("--window_size", type=int, default=100, help="Sliding window size (default: 100)")
    parser.add_argument("--phyl-tree", action="store_false", help="Generate MSA and Phylogenetic Tree (default: False)")
    parser.add_argument("--plots", action="store_false", help="Generate statistical plots (default: False)")
    parser.add_argument("--only-stats", action="store_false", default=False,
                        help="Generate only XLSX statistics (default: True)")

    args = parser.parse_args()

    analyze_graph_and_alignments(
        threads=args.threads,
        tree_pars=args.tree_pars,
        tree_bs=args.tree_bs,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        quantile=args.quantile,
        top_n=args.top_n,
        window_size=args.window_size,
        phyl_tree=args.phyl_tree,
        plots=args.plots,
        only_stats=args.only_stats
    )
