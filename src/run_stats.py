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


# ---------------------------------------------------------------------------
# Sub-step sentinel files (relative to output_dir inside the container).
# Used by MineGraph.py (host side) to decide which --skip-* flags to pass.
# Kept here as a single source of truth so host and container agree.
# ---------------------------------------------------------------------------
SENTINEL = {
    "core_stats": "MineGraph_output/statistics/graph_stats.xlsx",
    "interactive": "MineGraph_output/plots/graph_top_N_interactive.html",
    "plots":       "MineGraph_output/plots/node_histogram_by_paths.png",
    "odgi":        "MineGraph_output/phylogenetics_msa/odgi_matrix.tsv",
    "tree_graph":  "MineGraph_output/phylogenetics_msa/graph_phylo_tree.nwk",
    "tree_msa":    "MineGraph_output/phylogenetics_msa/MSA_result.fasta",
    "tree_render": "MineGraph_output/phylogenetics_msa/graph_phylo_tree.png",
}


def analyze_graph_and_alignments(threads="16", tree_pars=10, tree_bs=1000,
                                 input_dir="/data/pggb_output",
                                 output_dir="/data/MineGraph_output", quantile=0.25, top_n=50,
                                 window_size=1000, phyl_tree=True, plots=True, only_stats=False,
                                 tree_type='graph', softcore_threshold=5, convert_gfa=False,
                                 bootstrap_method="felsenstein",
                                 ufboot_min_rep=1000, ufboot_max_rep=10000,
                                 ufboot_convergence=0.99, ufboot_batch=100,
                                 # ── intra-step resume flags ──────────────────────────────────
                                 skip_core_stats=False, skip_plots=False,
                                 skip_odgi=False, skip_tree=False):
    os.makedirs(output_dir, exist_ok=True)
    stats_dir = os.path.join(output_dir, "MineGraph_output", "statistics")
    graph_dir = os.path.join(output_dir, "MineGraph_output", "plots")
    msa_dir   = os.path.join(output_dir, "MineGraph_output", "phylogenetics_msa")
    vg_dir    = os.path.join(output_dir, "MineGraph_output", "gfa_convert")
    os.makedirs(stats_dir, exist_ok=True)
    os.makedirs(graph_dir, exist_ok=True)
    os.makedirs(msa_dir,   exist_ok=True)
    os.makedirs(vg_dir,    exist_ok=True)

    # Locate required input files
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

    # ── Always parse the GFA — needed by every downstream sub-step ────────────
    nodes, edges, paths = get_nodes_edges_paths(gfa_file)
    nodes_num, edges_num, avg_degree, num_paths, edges_count = get_graph_stats(nodes, edges, paths)
    node_count_hist = get_nodes_counts_histogram(paths)
    count_dist = draw_distribution(
        node_count_hist,
        os.path.join(stats_dir, "graph_Node_Plot_frequency"),
    )

    # ── Sub-step 5a: core statistics (distribution plot + XLSX) ───────────────
    if skip_core_stats:
        print("[RESUME] Sub-step 5a (core stats) — already done, skipping.")
    else:
        if not only_stats:
            print("Generating graphical representations...")
            draw_graph_as_circle_interactive(
                nodes, edges_count, node_count_hist,
                os.path.join(graph_dir, "graph_top_N_interactive"), top_n,
            )
        extract_graph_stats(
            node_count_hist, edges_count, nodes_num, edges_num, avg_degree, num_paths,
            os.path.join(stats_dir, "graph_stats"), vcf_file, count_dist, softcore_threshold,
        )

    # Skip MSA and tree entirely when only_stats is set
    if only_stats:
        print("Skipping MSA and Phylogenetic tree generation (only stats mode enabled).")
        return

    # ── Sub-step 5b: statistical plots ────────────────────────────────────────
    if skip_plots:
        print("[RESUME] Sub-step 5b (statistical plots) — already done, skipping.")
    elif plots:
        plot_graph_statistics(gfa_file, graph_dir, window_size)
        plot_node_histogram_by_paths(
            count_dist, num_paths,
            os.path.join(graph_dir, "node_histogram_by_paths.png"),
            int(softcore_threshold),
        )

    # ── Sub-step 5c: consensus sequence + ODGI matrix ─────────────────────────
    if skip_odgi:
        print("[RESUME] Sub-step 5c (consensus + ODGI matrix) — already done, skipping.")
    else:
        cons_path, cons_seq, cycle_stat = extract_consensus_path(
            nodes, edges, node_count_hist, edges_count, paths,
            min_node_count=compute_min_node_count(node_count_hist, num_paths, quantile=quantile),
        )
        write_consensus_sequence(msa_dir, cons_seq)
        if convert_gfa:
            from utils import convert_gfa as _convert_gfa
            _convert_gfa(gfa_file, vg_dir)
        write_odgi_matrix(gfa_file, msa_dir)

    # ── Sub-step 5d: phylogenetic tree ────────────────────────────────────────
    if not phyl_tree:
        return

    odgi_tsv = os.path.join(msa_dir, "odgi_matrix.tsv")

    if skip_tree:
        print("[RESUME] Sub-step 5d (phylogenetic tree) — already done, skipping.")
        return

    if tree_type == "graph":
        if tree_pars != 10:
            print(
                f"[WARN] --tree_pars ({tree_pars}) is only used with --tree_type msa "
                f"(RAxML parsimony starts). It has no effect on the graph-based UPGMA tree."
            )
        print(
            f"Generating GRAPH-BASED phylogenetic tree from ODGI matrix "
            f"({tree_bs} bootstrap replicates)..."
        )
        if not os.path.exists(odgi_tsv):
            print("[ERROR] ODGI matrix not found and --skip-odgi was not set. "
                  "Re-run without --skip-odgi to regenerate it.")
            return
        from utils import generate_graph_based_tree, render_phylotree
        res = generate_graph_based_tree(
            odgi_tsv,
            output_prefix=os.path.join(msa_dir, "graph_phylo_tree"),
            n_bootstrap=tree_bs,
            n_threads=int(threads),
            bootstrap_method=bootstrap_method,
            ufboot_min_rep=ufboot_min_rep,
            ufboot_max_rep=ufboot_max_rep,
            ufboot_convergence=ufboot_convergence,
            ufboot_batch=ufboot_batch,
        )
        render_phylotree(res["tree_file"], os.path.join(msa_dir, "graph_phylo_tree"))

    else:
        # MSA / RAxML path
        print("Generating MSA and Phylogenetic Tree...")
        output_fasta = os.path.join(msa_dir, "MSA_result.fasta")
        if not skip_odgi:   # MSA depends on MAF conversion done in 5c
            maf_to_msa(maf_file, output_fasta)
        elif not os.path.exists(output_fasta):
            print("[ERROR] MSA FASTA not found and --skip-odgi was set. "
                  "Re-run without --skip-odgi.")
            return

        raxml_prefix = os.path.join(output_dir, "graph")
        raxml_output = f"{raxml_prefix}.raxml.support"
        print(
            f"Generating phylogenetic tree with RAxML "
            f"({threads} threads, pars={tree_pars}, bs={tree_bs})..."
        )
        raxml_command = [
            "raxml-ng", "--all", "--msa", output_fasta,
            "--model", "GTR+G",
            "--tree", f"pars{{{tree_pars}}}",
            "--bs-trees", str(tree_bs),
            "--threads", str(threads),
            "--prefix", raxml_prefix,
            "--redo",
        ]
        subprocess.run(raxml_command, check=True)

        if os.path.exists(raxml_output):
            render_phylotree(raxml_output, os.path.join(msa_dir, "graph_phylo_tree"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Graph and Alignment Analysis")
    parser.add_argument("--threads",     default="16",  help="Number of threads to use")
    parser.add_argument("--tree_pars",   type=int, default=10, help="Number of parsimonious trees")
    parser.add_argument("--tree_bs",     type=int, default=1000,
                        help="Bootstrap replicates: for --tree_type graph, Felsenstein bootstrap "
                             "replicates on the PAV matrix (default: 1000); "
                             "for --tree_type msa, RAxML --bs-trees.")
    parser.add_argument("--input_dir",   default="/data/pggb_output",    help="Input directory")
    parser.add_argument("--output_dir",  default="/data/MineGraph_output", help="Output directory")
    parser.add_argument("--quantile",    type=float, default=0.25,
                        help="Consensus nodes top percentage (default: 25%%)")
    parser.add_argument("--top_n",       type=int, default=50,
                        help="Top N nodes sizes to visualize")
    parser.add_argument("--window_size", type=int, default=1000,
                        help="Sliding window size (default: 1000)")
    parser.add_argument("--phyl-tree",   action="store_true", dest="phyl_tree",
                        help="Generate Phylogenetic Tree (default: False)")
    parser.add_argument("--plots",       action="store_true",
                        help="Generate statistical plots (default: False)")
    parser.add_argument("--only-stats",  action="store_true", dest="only_stats",
                        help="Generate only XLSX statistics, skip MSA and tree (default: False)")
    parser.add_argument("--tree_type",   choices=["graph", "msa"], default="graph",
                        help="'graph' (UPGMA from ODGI matrix, default) or 'msa' (RAxML on MSA).")
    parser.add_argument("--softcore_threshold", type=int, default=5,
                        help="Percentage threshold for unique/core node classification (default: 5).")
    parser.add_argument("--convert-gfa", action="store_true", dest="convert_gfa", default=False,
                        help="Convert GFA to VG and FASTA (default: off). Requires vg.")
    parser.add_argument("--bootstrap-method", dest="bootstrap_method",
                        choices=["felsenstein", "ufboot"], default="felsenstein",
                        help="Bootstrap algorithm: 'felsenstein' (default) or 'ufboot' (UFBoot-PAV).")
    parser.add_argument("--ufboot-min-rep",     type=int,   dest="ufboot_min_rep",     default=1000)
    parser.add_argument("--ufboot-max-rep",     type=int,   dest="ufboot_max_rep",     default=10000)
    parser.add_argument("--ufboot-convergence", type=float, dest="ufboot_convergence", default=0.99)
    parser.add_argument("--ufboot-batch",       type=int,   dest="ufboot_batch",       default=100)

    # ── Intra-step resume flags (set by MineGraph.py based on existing outputs) ──
    parser.add_argument("--skip-core-stats", action="store_true", dest="skip_core_stats",
                        default=False,
                        help="[Resume] Skip core stats XLSX + interactive viz (already done).")
    parser.add_argument("--skip-plots",      action="store_true", dest="skip_plots",
                        default=False,
                        help="[Resume] Skip statistical plots (already done).")
    parser.add_argument("--skip-odgi",       action="store_true", dest="skip_odgi",
                        default=False,
                        help="[Resume] Skip consensus sequence + ODGI matrix (already done).")
    parser.add_argument("--skip-tree",       action="store_true", dest="skip_tree",
                        default=False,
                        help="[Resume] Skip phylogenetic tree generation (already done).")

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
        only_stats=args.only_stats,
        tree_type=args.tree_type,
        softcore_threshold=args.softcore_threshold,
        convert_gfa=args.convert_gfa,
        bootstrap_method=args.bootstrap_method,
        ufboot_min_rep=args.ufboot_min_rep,
        ufboot_max_rep=args.ufboot_max_rep,
        ufboot_convergence=args.ufboot_convergence,
        ufboot_batch=args.ufboot_batch,
        skip_core_stats=args.skip_core_stats,
        skip_plots=args.skip_plots,
        skip_odgi=args.skip_odgi,
        skip_tree=args.skip_tree,
    )
