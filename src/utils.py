import csv
import os
import re
import subprocess
from collections import defaultdict
import igraph as ig
import numpy as np
import pandas as pd
from ete3 import Tree
from ete3.treeview import TreeStyle
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
from cyvcf2 import VCF
from matplotlib import rcParams

# Global font and text settings for publication
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Liberation Sans', 'Arial'],
    'font.size': 18,
    'axes.titlesize': 22,
    'axes.labelsize': 20,
    'legend.fontsize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
})


####################################
"""

Utility script to handle bunch of macro function and analysis

"""
####################################


# Extract nodes, edges, and paths from a GFA file
def get_nodes_edges_paths(graph_file):
    """
    a function that given the graph file path it opens it and processes the data to extract the nodes, edges, paths
    nodes as a dict() where key is a node id as number (the graph is topologically sorted) ,
                            value is the corresponding DNA segment as string
    edges as dict() where key is a tuple of nodes ids with the strand of each node e.g. (1+, 2-)  means
                    an edge from the forward reading node 1 DNA fragment to the backward reading node 2 DNA fragment
                    if node 1 is ACT , node 2 is GAC so the edge still connects both node 1 and 2 but when walking on
                    the graph we get ACTCAG if it was (1-, 2-) so we get TCACAG and so on,
                        values are
    :param graph_file: graph file path
    :return:
    """
    print("Extracting nodes, edges, and paths from GFA file...")
    with open(graph_file, 'r') as graph_in:
        lines = graph_in.readlines()
        nodes = {}
        edges = {}
        paths = {}
        for line in lines:
            if line.startswith('L'):
                data = line[2:].split('\t')[:-1]
                n0, st0, n1, st1 = data[0], data[1], data[2], data[3]
                data = [n0 + st0, n1 + st1]
                if data[0] not in edges:
                    edges[data[0]] = [data[1]]
                else:
                    edges[data[0]].append(data[1])
            if line.startswith('S'):
                data = line[2:].split('\t')
                nodes[data[0]] = data[1]
            if line.startswith('P'):
                data = line[2:].split('\t')
                paths[data[0]] = [data[1].split(',')]
                paths[data[0]].append([str(sum(len(nodes[i[:-1]][:-1]) for i in paths[data[0]][0]))])
    print("Extraction complete.")
    return nodes, edges, paths

# Generate sequence based on wanted nodes and graph nodes reference
def generate_seq(wanted_nodes_list, nodes):
    seq = ''
    for link in wanted_nodes_list:
        if link.endswith('+'):
            seq += nodes[link[:-1]][:-1]
        else:
            seq += nodes[link[:-1]][::-1][:-1]
    return ''.join(seq.split())


def truncate_alignment_sequences(alignments):
    # Find the length of the shortest sequence
    min_length = min(len(seq) for seq in alignments.values())

    # Truncate all sequences to the minimum length
    for seq_name in alignments:
        alignments[seq_name] = alignments[seq_name][:min_length]

    return alignments


# Calculate histogram of node support
def get_nodes_counts_histogram(paths):
    """
    paths: dict[path_name] -> (list_of_steps, ...)
           where list_of_steps looks like ['123+', '124-', '123+', ...]
    Returns:
      my_hist[node_id] = [[fwd_step_count, rev_step_count], path_presence_count]
    """
    print("Generating data for nodes counts ...")
    my_hist = {}

    for path_name, (node_steps, *_) in paths.items():
        # 1) STEP COUNTS (counts every visit; keeps your fwd/rev tallies)
        for step in node_steps:
            node_id = step[:-1]
            direction = step[-1]
            if node_id not in my_hist:
                my_hist[node_id] = [[0, 0], 0]
            if direction == '+':
                my_hist[node_id][0][0] += 1
            else:
                my_hist[node_id][0][1] += 1

        # 2) PATH PRESENCE (counts each path at most once per node; ignores strand)
        nodes_in_this_path = {s[:-1] for s in node_steps}  # dedup by node_id
        for node_id in nodes_in_this_path:
            if node_id not in my_hist:
                my_hist[node_id] = [[0, 0], 0]
            my_hist[node_id][1] += 1

    print("Done.")
    return my_hist


# Compute graph statistics such as node and edge count
def get_graph_stats(nodes, edges, paths):
    print("Computing graph statistics...")
    nodes_num = len(nodes)
    edges_num = sum([len(edges[key]) for key in edges.keys()])
    avg_degree = sum(len(i) for i in edges.values()) / nodes_num
    num_paths = len(paths)
    edges_count = {}
    for path in paths:
        for i in range(len(paths[path][0]) - 2):
            edge = (paths[path][0][i][:-1], paths[path][0][i + 1][:-1])
            edges_count[edge] = edges_count.get(edge, 0) + 1
    print("Graph statistics computation complete.")
    return nodes_num, edges_num, avg_degree, num_paths, edges_count

# Create an igraph object from nodes and edges
def creat_igraph(nodes, edges_count):
    print("Creating igraph object from nodes and edges...")
    G = ig.Graph(directed=True)
    G.add_vertices(list(nodes.keys()))
    G.add_edges([(edge[0], edge[1]) for edge in edges_count.keys()])
    print("igraph creation complete.")
    return G

# Render a phylogenetic tree using ete3
def render_phylotree(tree_file, output_name):
    os.environ["QT_QPA_PLATFORM"] = "offscreen"  # Use offscreen mode for rendering
    print("Rendering the Phylogenetic Tree")
    with open(tree_file) as my_tree:
        t = Tree(my_tree.read())
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_support = True
        t.render('{}.png'.format(output_name), tree_style=ts, dpi=600)
    print("Rendering completed.")

# Process PAF file, extracting alignment sequences with specific reference
def generate_MSA_using_graph(num_paths, nodes_num, paths, nodes, output_fasta):
    print("Processing PAF file with reference...")
    # Create a boolean matrix for path presence
    x = np.full((num_paths, nodes_num), False, dtype=bool)
    cnt = 0
    ord_names_list = []

    # Precompute which paths contain each node (optimized loop)
    path_node_map = {i: set(paths[i][0]) for i in paths.keys()}

    for i in paths.keys():
        ord_names_list.append(i)
        for n in paths[i][0]:
            x[cnt][int(n[:-1]) - 1] = True
        cnt += 1

    # Generate the FASTA sequence output
    seeq = []
    for cnt, i in enumerate(ord_names_list):
        seq_lines = [f">{i}\n"]
        current_path_nodes = path_node_map[i]  # Cached set of nodes in the current path
        seq = []

        for j in range(nodes_num):
            node_id = str(j + 1)
            if not x[cnt][j]:
                seq.append("-" * len(nodes[node_id][:-1]))  # Node is absent
            else:
                # Node is present, determine strand
                if f"{node_id}+" in current_path_nodes:
                    seq.append(nodes[node_id][:-1])  # Forward strand
                else:
                    seq.append(nodes[node_id][:-1][::-1])  # Reverse complement strand

        seq_lines.append("".join(seq))
        seq_lines.append("\n")
        seeq.extend(seq_lines)
    with open(output_fasta, 'w') as f:
        f.writelines(seeq)

# Visualization and output functions
def save_node_distribution(my_hist, graph_name):
    """
    Saves the node count distribution information to a CSV file.

    Args:
        my_hist (dict): Dictionary where keys are node counts and values are their frequencies.
        graph_name (str): Name of the graph to label and save the output file.
    """
    output_csv = f"{graph_name}_Node_Distribution.csv"
    with open(output_csv, mode="w", newline="") as csvfile:
        fieldnames = ["Node Count", "Frequency"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for count, frequency in sorted(my_hist.items()):
            writer.writerow({"Node Count": count, "Frequency": frequency})

    print(f"Node distribution data saved as {output_csv}")


def draw_distribution(node_count_hist, graph_name):
    """
    Plots the distribution of node degrees in a histogram with a logarithmic y-axis.

    Args:
        node_count_hist (dict): Dictionary where keys are nodes and values are their counts.
        graph_name (str): Name to label the graph and save the output file.
    """
    print(f"Drawing Nodes frequency (node count in all paths)...")

    # Create histogram data
    my_hist = {}
    graph_nodes_counts = [node_count_hist[i][1] for i in node_count_hist.keys()]
    for count in graph_nodes_counts:
        if count in my_hist:
            my_hist[count] += 1
        else:
            my_hist[count] = 1

    # Save node distribution data to CSV
    save_node_distribution(my_hist, graph_name)

    # Generate bins up to the max count value
    max_node_count = max(my_hist.keys())
    bins = range(0, max_node_count + 1)

    # Plot histogram with logarithmic y-axis
    plt.figure(figsize=(12, 7))
    plt.hist(list(my_hist.keys()), bins=bins, weights=list(my_hist.values()), edgecolor='black', alpha=0.7)
    plt.yscale('log')  # Apply logarithmic scale to y-axis
    plt.title("Node Frequency Histogram (Log Scale)")
    plt.xlabel("Node Count (How many times a node appears in all paths)")
    plt.ylabel("Frequency (Log Scale)")
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Add a description text box to the plot
    description = (
        "This histogram shows the frequency of node count values in the graph. "
        "Each bin represents a unique node count (x-axis), which is the number of paths a node is part of. "
        "The y-axis (logarithmic scale) indicates how many nodes share the same count value."
    )
    plt.gcf().text(0.5, -0.15, description, wrap=True, ha='center', alpha=0.8)

    # Save the plot
    output_path = f"{graph_name}_Node_Distribution.png"
    plt.savefig(output_path, dpi=600, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"Node frequency histogram saved as {output_path}")
    return my_hist


def save_all_nodes_info(largest_nodes_indices, nodes, node_sizes, node_count_hist, graph_name):
    """
    Saves information about all nodes in a CSV file.

    Args:
        largest_nodes_indices (list): List of all node indices, sorted by size.
        nodes (dict): Dictionary of nodes and their associated sequences.
        node_sizes (list): List of node sizes.
        node_count_hist (dict): Dictionary containing the count for each node.
        graph_name (str): Name of the graph to label and save the output file.
    """
    output_csv = f"{graph_name}_all_nodes_info.csv"
    all_nodes_info = []

    for idx in largest_nodes_indices:
        node_id = list(nodes.keys())[idx]
        size = node_sizes[idx]
        count = node_count_hist.get(node_id, [0, 0, 0])[1]  # Total count (forward + reverse visits)
        forward_visits, reverse_visits = node_count_hist.get(node_id, [[0, 0], 0])[0]

        all_nodes_info.append({
            "Node ID": node_id,
            "Size": size,
            "Count": count,
            "Forward Visits": forward_visits,
            "Reverse Visits": reverse_visits
        })

    # Write data to CSV
    with open(output_csv, mode="w", newline="") as csvfile:
        fieldnames = ["Node ID", "Size", "Count", "Forward Visits", "Reverse Visits"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_nodes_info)

    print(f"All nodes information saved as {output_csv}")


def draw_graph_as_circle(nodes, edges_count, node_count_hist, graph_name):
    print("Drawing graph as a circular layout...")

    # Create igraph graph with nodes and edges
    edges_list = list(edges_count.keys())
    g = ig.Graph(directed=True)
    g.add_vertices(list(nodes.keys()))
    g.add_edges(edges_list)

    # Normalize node sizes to the range [0, 1] and then apply scaling factor
    node_sizes = [len(nodes[node]) - 1 for node in nodes]
    max_size = max(node_sizes) if node_sizes else 1
    min_size = min(node_sizes) if node_sizes else 1

    # Scale nodes to range [0, 1]
    normalized_node_sizes = [(size - min_size) / (max_size - min_size) if max_size > min_size else 0 for size in
                             node_sizes]
    scaling_factor = 180  # Choose an appropriate scaling factor for node visibility
    min_node_display_size = 5  # Minimum size for visibility

    # Apply scaling factor
    scaled_node_sizes = [size * scaling_factor + min_node_display_size for size in normalized_node_sizes]
    for i, size in enumerate(scaled_node_sizes):
        g.vs[i]["size"] = size

    # Select the top 50 largest nodes for labeling
    largest_nodes_indices = sorted(range(len(node_sizes)), key=lambda i: node_sizes[i], reverse=True)
    top_50 = largest_nodes_indices[:50]
    vertex_labels = [""] * len(node_sizes)  # Default to no label

    # Set labels only for the top 50 largest nodes
    for i in top_50:
        vertex_labels[i] = g.vs[i]["name"]  # Assuming node name is the label

    g.vs["label"] = vertex_labels  # Apply labels to vertices

    # Save top 50 nodes information to CSV
    #save_all_nodes_info(largest_nodes_indices, nodes, node_sizes, node_count_hist, graph_name)

    # Normalize edge weights within the range [0, 1]
    edge_weights = list(edges_count.values())
    max_weight = max(edge_weights) if edge_weights else 1
    min_weight = min(edge_weights) if edge_weights else 1

    # Scale edges to range [0, 1]
    normalized_edge_weights = [(weight - min_weight) / (max_weight - min_weight) if max_weight > min_weight else 0 for
                               weight in edge_weights]
    edge_scaling_factor = 0.5  # Reduced scaling for edge width
    min_edge_display_weight = 0.2  # Minimum edge width for visibility

    # Apply scaling factor
    scaled_edge_weights = [weight * edge_scaling_factor + min_edge_display_weight for weight in normalized_edge_weights]
    for i, weight in enumerate(scaled_edge_weights):
        g.es[i]["weight"] = weight

    # Compute circular layout
    layout = g.layout("circle")

    # Visual style settings for the plot
    visual_style = {
        "vertex_size": g.vs["size"],
        "vertex_label": g.vs["label"],
        "vertex_color": "skyblue",
        "edge_width": g.es["weight"],
        "layout": layout,
        "bbox": (2400, 1800),
        "margin": 100,
        "autocurve": True,
        "edge_curved": False,
        "edge_arrow_size": 0.1,
        "arrow_width": 0.1,
    }

    # Save the graph as an image
    output_path = f"{graph_name}_igraph_plot_stretched.png"
    expl_path = f"{graph_name}_igraph_plot_stretched_README.txt"
    with open(expl_path, "w") as f:
        f.write("Graph Visualization\n The graph nodes were scattered around a circle\n"
                "Node sizes indicate sequence-node size\n"
                "Edge widths indicate connection strength\nTop 50 nodes labeled")
    ig.plot(g, output_path, **visual_style)
    print(f"Graph saved as {output_path}")


def draw_graph_as_circle_interactive(nodes, edges_count, node_count_hist, graph_name, top_n=1000):
    """
    Generates an interactive HTML visualization of the graph where node info is displayed on hover.

    Args:
        nodes (dict): Dictionary of nodes and their associated sequences.
        edges_count (dict): Dictionary of edges and their weights.
        node_count_hist (dict): Dictionary of node counts (forward, reverse, total).
        graph_name (str): Name of the graph for saving the HTML file.
        top_n (int): Number of top nodes to visualize based on size.
    """
    print("Generating interactive HTML visualization of the graph...")

    # Calculate node sizes
    node_sizes = [len(nodes[node]) - 1 for node in nodes]

    # Sort nodes by size and select the top N
    largest_nodes_indices = sorted(range(len(node_sizes)), key=lambda i: node_sizes[i], reverse=True)[:top_n]
    selected_nodes = {list(nodes.keys())[i]: nodes[list(nodes.keys())[i]] for i in largest_nodes_indices}

    # Filter edges to include only those connecting the selected nodes
    selected_edges_count = {
        edge: weight
        for edge, weight in edges_count.items()
        if edge[0] in selected_nodes and edge[1] in selected_nodes
    }

    # Initialize the pyvis Network
    net = Network(height="800px", width="100%", directed=True, notebook=False)

    # Add selected nodes with scaled sizes but display original sizes in tooltips
    for idx, (node_id, node_seq) in enumerate(selected_nodes.items()):
        scaled_size = (node_sizes[largest_nodes_indices[idx]] - min(node_sizes)) / (max(node_sizes) - min(node_sizes))
        scaled_size = scaled_size * 180 + 5  # Scale size for clarity
        original_size = node_sizes[largest_nodes_indices[idx]]  # Original size
        count = node_count_hist.get(node_id, [0, 0, 0])[1]
        forward_visits, reverse_visits = node_count_hist.get(node_id, [[0, 0], 0])[0]

        tooltip_text = (
            f"Node ID: {node_id}\n"
            f"Original Size: {original_size}\n"
            f"Scaled Size (for visualization): {scaled_size:.2f}\n"
            f"Count: {count}\n"
            f"Forward Visits: {forward_visits}\n"
            f"Reverse Visits: {reverse_visits}"
        )

        net.add_node(
            node_id,
            label=node_id,
            size=scaled_size,
            title=tooltip_text,  # Display original size in tooltip
        )

    # Add selected edges
    for edge, weight in selected_edges_count.items():
        net.add_edge(
            edge[0],
            edge[1],
            value=weight,
            title=f"Edge Weight: {weight}"
        )

    # Customize layout
    net.repulsion(node_distance=300, central_gravity=0.1, spring_length=100, spring_strength=0.01)

    # Generate the HTML file
    output_html = os.path.join(f"{graph_name}_interactive.html")
    net.write_html(output_html, notebook=False)

    # Add explanation text to the generated HTML file
    explanation_text = """
    <div style="padding: 10px; font-family: Arial, sans-serif; background-color: #f9f9f9; border: 1px solid #ddd; margin-bottom: 20px;">
        <h2>Graph Visualization Explanation</h2>
        <p>
            This interactive graph represents the top N nodes layout of nodes and edges, offering insights into the underlying relationships
            and structures. Each node corresponds to a sequence or entity in the dataset, with its size proportional to the sequence length
            or another relevant metric.
        </p>
        <p>
            The edges indicate connections or interactions between nodes, with their widths scaled to represent the strength or weight of these
            connections. Hovering over a node reveals detailed information, including:
        </p>
        <ul>
            <li><b>Node ID:</b> A unique identifier for the node.</li>
            <li><b>Size:</b> The scaled size of the node, reflecting its relative importance or size.</li>
            <li><b>Count:</b> The total number of occurrences or interactions for this node.</li>
            <li><b>Forward and Reverse Visits:</b> The directional frequency of interactions involving the node.</li>
        </ul>
        <p>
            Use the interactive features to explore the graph: zoom in, pan around, and hover over elements to uncover details. 
            This visualization provides a powerful way to understand the complex interplay of elements within the dataset.
        </p>
    </div>
    """

    # Insert explanation text into the HTML file
    with open(output_html, "r") as file:
        html_content = file.read()

    updated_html = html_content.replace(
        "<body>",
        f"<body>{explanation_text}"  # Add explanation text right after the body tag
    )

    with open(output_html, "w") as file:
        file.write(updated_html)

    print(f"Interactive graph visualization saved to {output_html}")


def plot_node_histogram_by_paths(
    node_histogram,
    num_paths,
    output_file="node_histogram_by_paths.png",
    softcore_threshold: int = 5
):
    """
    Plot the full node presence distribution colored by category:
      - Unique (≤ threshold%)
      - Shared (middle range)
      - Core (≥ 100−threshold%)

    Automatically switches to log scale on Y-axis if node counts
    span several orders of magnitude.
    """
    import math
    from matplotlib.patches import Patch


    # --- Determine cutoffs ---
    low_cut = math.ceil((softcore_threshold / 100) * num_paths)
    high_cut = math.floor(((100 - softcore_threshold) / 100) * num_paths)

    # --- Prepare data ---
    x = list(sorted(node_histogram.keys()))
    y = np.array([node_histogram[c] for c in x], dtype=float)

    # --- Assign colors per bar ---
    colors = []
    for c in x:
        if c <= low_cut:
            colors.append("#1f77b4")  # Unique (blue)
        elif c >= high_cut:
            colors.append("#2ca02c")  # Core (green)
        else:
            colors.append("#ff7f0e")  # Shared (orange)

    # --- Check if we need log scale ---
    if np.max(y) / max(np.min(y[y > 0]), 1) > 1e3:
        use_log = True
        print("Wide value range detected — using logarithmic Y scale.")
    else:
        use_log = False

    # --- Plot the histogram ---
    plt.figure(figsize=(10, 4))
    plt.bar(x, y, color=colors, edgecolor="none", width=1.0)

    plt.xlabel("Number of Paths (Haplotypes)")
    plt.ylabel("Number of Nodes")
    #plt.title(f"Node Distribution by Path Presence (Threshold = {softcore_threshold}%)", fontsize=12)

    if use_log:
        plt.yscale("log")

    # --- Add legend ---
    legend_elements = [
        Patch(facecolor="#1f77b4", edgecolor="black", label=f"Unique (≤ {softcore_threshold}%)"),
        Patch(facecolor="#ff7f0e", edgecolor="black", label=f"Shared ({softcore_threshold}%–{100 - softcore_threshold}%)"),
        Patch(facecolor="#2ca02c", edgecolor="black", label=f"Core (≥ {100 - softcore_threshold}%)")
    ]
    plt.legend(handles=legend_elements, loc="upper right", frameon=False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Colored node distribution histogram saved to {output_file}")





def convert_gfa(gfa_file, output_dir="/data/MineGraph_output/"):
    """Convert GFA file to VG and FASTA format, saving outputs to specified directory."""

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # VG output file path
    vg_output_file = os.path.join(output_dir, "gfa_to_vg.vg")
    vg_command = ["vg", "convert", "-g", gfa_file]

    # Redirect VG output to specified file
    print(f"Converting GFA to VG: {' '.join(vg_command)} > {vg_output_file}")
    with open(vg_output_file, 'w') as vg_out:
        subprocess.run(vg_command, stdout=vg_out, check=True)

    # FASTA output file path
    fasta_output_file = os.path.join(output_dir, "gfa_to_fasta.fasta")
    gfa2fasta_command = ["python", "gfa2fasta.py", gfa_file, fasta_output_file]

    print(f"Converting GFA to FASTA: {' '.join(gfa2fasta_command)}")
    subprocess.run(gfa2fasta_command, check=True)

    print(f"Conversion completed. VG saved to {vg_output_file}, FASTA saved to {fasta_output_file}")



def parse_paf_line(line):
    fields = line.strip().split()
    query_seq_name = fields[0]
    query_start = int(fields[2])
    query_end = int(fields[3])
    strand = fields[4]
    target_seq_name = fields[5]
    target_start = int(fields[7])
    target_end = int(fields[8])
    cigar = fields[-1].split(':')[-1]

    return {
        "query_seq_name": query_seq_name,
        "query_start": query_start,
        "query_end": query_end,
        "strand": strand,
        "target_seq_name": target_seq_name,
        "target_start": target_start,
        "target_end": target_end,
        "cigar": cigar
    }


def restore_alignment(cigar, ref_seq, query_seq):
    aligned_ref = []
    aligned_query = []
    ref_index = 0
    query_index = 0

    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

    for length, operation in cigar_tuples:
        length = int(length)

        if operation in {'M', '=', 'X'}:
            aligned_ref.append(ref_seq[ref_index:ref_index + length])
            aligned_query.append(query_seq[query_index:query_index + length])
            ref_index += length
            query_index += length

        elif operation == 'I':
            aligned_ref.append('-' * length)
            aligned_query.append(query_seq[query_index:query_index + length])
            query_index += length

        elif operation == 'D':
            aligned_ref.append(ref_seq[ref_index:ref_index + length])
            aligned_query.append('-' * length)
            ref_index += length

    return ''.join(aligned_ref), ''.join(aligned_query)



def parse_vcf_stats(file_path, output_file="/data/MineGraph_output/polymorphism_summary.txt"):
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True)

    stats = {
        "samples": None,
        "total_records": None,
        "no_alts": None,
        "snp_count": None,
        "mnp_count": None,
        "indel_count": None,
        "other_count": None,
        "multiallelic_sites": None,
        "multiallelic_snp_sites": None,
        "transition_transversion": {"transitions": None, "transversions": None, "ratio": None},
        "substitution_types": {},
    }

    with open(file_path, "r") as f:
        for line in f:
            # Parse summary stats
            if line.startswith("SN"):
                parts = line.strip().split("\t")
                key = parts[2]
                value = int(parts[3])

                if key == "number of samples:":
                    stats["samples"] = value
                elif key == "number of records:":
                    stats["total_records"] = value
                elif key == "number of no-ALTs:":
                    stats["no_alts"] = value
                elif key == "number of SNPs:":
                    stats["snp_count"] = value
                elif key == "number of MNPs:":
                    stats["mnp_count"] = value
                elif key == "number of indels:":
                    stats["indel_count"] = value
                elif key == "number of others:":
                    stats["other_count"] = value
                elif key == "number of multiallelic sites:":
                    stats["multiallelic_sites"] = value
                elif key == "number of multiallelic SNP sites:":
                    stats["multiallelic_snp_sites"] = value

            # Parse transitions/transversions ratio
            elif line.startswith("TSTV"):
                parts = line.strip().split("\t")
                stats["transition_transversion"]["transitions"] = int(parts[2])
                stats["transition_transversion"]["transversions"] = int(parts[3])
                stats["transition_transversion"]["ratio"] = float(parts[4])

            # Parse substitution types
            elif line.startswith("ST"):
                parts = line.strip().split("\t")
                substitution_type = parts[2]
                count = int(parts[3])
                stats["substitution_types"][substitution_type] = count

    # Write the polymorphism summary to a file
    with open(output_file, "w") as out_f:
        out_f.write("Polymorphism Summary:\n")
        out_f.write(f"Total Samples: {stats['samples']}\n")
        out_f.write(f"Total Records: {stats['total_records']}\n")
        out_f.write(f"SNPs: {stats['snp_count']}\n")
        out_f.write(f"MNPs: {stats['mnp_count']}\n")
        out_f.write(f"Indels: {stats['indel_count']}\n")
        out_f.write(f"Other Types: {stats['other_count']}\n")
        out_f.write(f"Multiallelic Sites: {stats['multiallelic_sites']}\n")
        out_f.write(f"Multiallelic SNP Sites: {stats['multiallelic_snp_sites']}\n")
        out_f.write(f"Transition/Transversion Ratio: {stats['transition_transversion']['ratio']}\n\n")

        out_f.write("Substitution Types:\n")
        for sub_type, count in stats["substitution_types"].items():
            out_f.write(f"{sub_type}: {count}\n")


def parse_paf(paf_file):
    alignments = defaultdict(list)
    with open(paf_file, 'r') as f:
        for line in f:
            alignment = parse_paf_line(line)
            alignments[alignment["query_seq_name"]].append(alignment)
    return alignments


def process_path(path_name, path_data, nodes, alignments, max_length):
    nucleotide_counts = np.zeros((max_length, 4), dtype=int)  # A, C, G, T counts per position
    node_ids = path_data  # List of nodes for this path
    path_sequence = generate_seq(node_ids, nodes)

    if path_name in alignments:
        for alignment in alignments[path_name]:
            query_start = alignment["query_start"]
            query_end = alignment["query_end"]
            target_start = alignment["target_start"]

            for i in range(query_start, min(query_end, len(path_sequence))):
                position = target_start + i - query_start
                if position < max_length:
                    base = path_sequence[i]
                    if base == 'A':
                        nucleotide_counts[position][0] += 1
                    elif base == 'C':
                        nucleotide_counts[position][1] += 1
                    elif base == 'G':
                        nucleotide_counts[position][2] += 1
                    elif base == 'T':
                        nucleotide_counts[position][3] += 1
    return nucleotide_counts


def compute_min_node_count(node_count_hist,paths_num, quantile=0.25):
    """
    Compute the optimal MIN_NODE_COUNT based on node count distribution.

    Parameters:
        node_count_hist (dict): Dictionary where keys are node IDs and values are lists of the form [[forward_count, reverse_count], total_count].
        quantile (float): Quantile value to set the threshold Consensus nodes For example, 100 means nodes must appear in all paths.

    Returns:
        int: Optimal minimum node count.
        :param paths_num: graph paths total number
    """
    # # Step 1: Compute total counts for each node
    # total_counts = []
    # for count in node_count_hist.values():
    #     if isinstance(count, list) and len(count) == 2 and isinstance(count[0], list) and len(count[0]) == 2:
    #         # Extract total count
    #         total_counts.append(count[1])
    #     else:
    #         raise ValueError(f"Invalid format in node_count_hist value: {count}")
    #
    # # Step 2: Analyze the distribution of node counts
    # if not total_counts:
    #     raise ValueError("Node count histogram is empty or incorrectly formatted.")
    #
    # median = np.median(total_counts)
    # std_dev = np.std(total_counts)
    # quantile_threshold = math.ceil(quantile * paths_num)  # Top quantile nodes
    #
    # # Step 3: Set the threshold
    # # Use a combination of median, standard deviation, and quantile threshold
    # min_node_count = max(int(median + std_dev), int(quantile_threshold))
    #
    # return min_node_count
    return int(quantile * paths_num)




import networkx as nx
from collections import defaultdict
from Bio.Seq import Seq


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


from collections import defaultdict
import re

def extract_consensus_path(
        nodes,                    # {node_id(str): sequence(str)}
        edges,                    # unused (kept for compatibility)
        node_count_hist,          # {node_id: [[fwd, rev], path_presence]}
        edges_count,              # unused (kept for compatibility)
        paths,                    # {haplo_name: [oriented_nodes_list], ...}
        min_node_count=2,         # filter by #paths that contain node (not steps)
        majority_threshold=0.5    # unused here (kept for compatibility)
):
    """
    Simplified consensus builder with NO cycle handling and NO coordinates.
    Nodes are ordered purely by their IDs (topological order assumed == ID order).

    Returns:
        (consensus_path, consensus_sequence, meta)
          - consensus_path: list of oriented node IDs (e.g. '123+')
          - consensus_sequence: concatenated sequence (respecting orientation)
          - meta: small dict with notes about the run
    """

    # --- helpers ---
    def reverse_complement(seq: str) -> str:
        comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return seq.translate(comp)[::-1]

    def id_sort_key(nid: str):
        """
        Sort by trailing integer if present (e.g., 'u54' -> 54, '123' -> 123),
        else fallback to full string for stable lexicographic order.
        """
        s = str(nid)
        m = re.search(r'(\d+)$', s)
        return (0, int(m.group(1))) if m else (1, s)

    # --- Step 1: select valid nodes by path-presence ---
    valid_nodes = set()
    for nid, info in node_count_hist.items():
        # expect [[fwd, rev], presence]
        presence = 0
        if isinstance(info, (list, tuple)) and len(info) >= 2:
            presence = info[1]
        if presence >= min_node_count:
            valid_nodes.add(str(nid))

    # --- Step 2: lock dominant orientation per node ---
    dominant_orient = {}
    for nid, info in node_count_hist.items():
        nid = str(nid)
        if nid not in valid_nodes:
            continue
        fwd = rev = 0
        if isinstance(info, (list, tuple)) and len(info) >= 1 and isinstance(info[0], (list, tuple)) and len(info[0]) == 2:
            fwd, rev = info[0]
        dominant_orient[nid] = '+' if fwd >= rev else '-'

    # --- Step 3: order nodes purely by ID (topological-by-ID assumed) ---
    ordered_nids = sorted(valid_nodes, key=id_sort_key)

    # --- Step 4: build consensus path + sequence in that order ---
    consensus_path = []
    consensus_sequence = []

    for nid in ordered_nids:
        seq = nodes.get(nid, "")[:-1]
        if not seq:
            # skip if no sequence available
            continue
        strand = dominant_orient.get(nid, '+')
        final_seq = seq if strand == '+' else reverse_complement(seq)
        consensus_path.append(f"{nid}{strand}")
        consensus_sequence.append(final_seq)

    consensus_sequence = "".join(consensus_sequence)

    meta = {
        "cycles_processed": False,
        "ordering": "by_id",
        "n_nodes_in_path": len(consensus_path),
        "min_node_count": min_node_count,
    }

    return consensus_path, consensus_sequence, meta




def preprocess_edges(edges):
    """Preprocess the edges to map node IDs to their outgoing edges."""
    outgoing_edges_map = {}
    for edge, targets in edges.items():
        node_id = edge.rstrip("+-")  # Remove strand annotation
        if node_id not in outgoing_edges_map:
            outgoing_edges_map[node_id] = []
        outgoing_edges_map[node_id].extend(targets)
    return outgoing_edges_map

def reverse_complement(sequence):
    """Returns the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return sequence.translate(complement)[::-1]


def write_consensus_sequence(output_dir, cons_seq):
    """
    Writes the consensus sequence to a FASTA file in the specified output directory.

    Parameters:
        output_dir (str): Path to the directory where the output file will be saved.
        cons_seq (str): The consensus sequence to write.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define the output file path
    output_file = os.path.join(output_dir, "Consensus_sequence.fasta")

    # Write the consensus sequence to the file
    with open(output_file, "w") as f:
        f.write(">Consensus_sequence\n")
        f.write(cons_seq)

    print(f"Consensus sequence saved to {output_file}")


def draw_consensus_path(nodes, edges_count, consensus_path, node_count_hist, graph_name, top_n=1000):
    """
    Generates an interactive HTML visualization of the top_n largest nodes (by sequence length)
    and their edges in the consensus path, arranged dynamically using a circular layout.

    Args:
        nodes (dict): Dictionary of node IDs and their sequences.
        edges_count (dict): Dictionary of edge coverage counts (edge tuple -> count).
        consensus_path (list): List of node IDs forming the consensus path, sorted in order.
        node_count_hist (dict): Dictionary with node visit counts, structured as [[forward_count, reverse_count], total_count].
        graph_name (str): Name of the graph for saving the HTML file.
        top_n (int): Maximum number of nodes to display, selected by sequence length.
    """
    print("Generating consensus path visualization...")

    # Step 1: Sort nodes by sequence length in descending order
    sorted_nodes = sorted(consensus_path, key=lambda node: len(nodes.get(node, "")), reverse=True)[:top_n]

    # Step 2: Initialize NetworkX graph
    G = nx.Graph()
    for node in sorted_nodes:
        G.add_node(node)
    for node1, node2 in edges_count:
        if node1 in sorted_nodes and node2 in sorted_nodes:
            G.add_edge(node1, node2, weight=edges_count[(node1, node2)])

    # Step 3: Circular layout with bigger circle spacing
    scale_factor = 10000 + (top_n * 50)
    pos = nx.circular_layout(G, scale=scale_factor)

    # Step 4: Initialize pyvis Network
    net = Network(height="1000px", width="100%", directed=True, notebook=False, bgcolor="white", font_color="black")

    # Step 5: Dynamic node sizes and colors using a color dictionary
    min_size, max_size = 50, 350  # Larger sizes
    sequence_lengths = [len(nodes[node[:-1]]) for node in sorted_nodes]
    seq_min, seq_max = min(sequence_lengths), max(sequence_lengths)

    # Define a simple color palette for distinct ranges
    size_color_map = {
        0: "rgba(255, 0, 0, 0.8)",     # Red
        1: "rgba(255, 165, 0, 0.8)",   # Orange
        2: "rgba(255, 255, 0, 0.8)",   # Yellow
        3: "rgba(0, 128, 0, 0.8)",     # Green
        4: "rgba(0, 0, 255, 0.8)",     # Blue
        5: "rgba(75, 0, 130, 0.8)",    # Indigo
        6: "rgba(148, 0, 211, 0.8)"    # Violet
    }

    def get_color(seq_len):
        """Map sequence lengths to color ranges."""
        range_step = (seq_max - seq_min) / len(size_color_map)
        index = min(int((seq_len - seq_min) / range_step), len(size_color_map) - 1)
        return size_color_map[index]

    for node in sorted_nodes:
        seq_len = len(nodes[node[:-1]])
        node_size = min_size + (seq_len - seq_min) / (seq_max - seq_min) * (max_size - min_size)
        color = get_color(seq_len)
        count = node_count_hist.get(node, [[0, 0], 0])[1]
        tooltip_text = (
            f"Node ID: {node}\n"
            f"Sequence Length: {seq_len}\n"
            f"Node Size (scaled): {node_size:.2f}\n"
            f"Node Count: {count}"
        )

        net.add_node(
            node,
            label=node,
            title=tooltip_text,
            size=node_size,
            x=pos[node][0],
            y=-pos[node][1],
            color=color,
            physics=False
        )

    # Step 6: Add edges with dynamically scaled thickness
    edge_min_width, edge_max_width = 20, 100  # Increased edge width range
    weight_min = min(edges_count.values(), default=1)
    weight_max = max(edges_count.values(), default=1)

    for edge in G.edges:
        weight = edges_count.get(edge, 1)
        scaled_width = edge_min_width + (weight - weight_min) / (weight_max - weight_min) * (edge_max_width - edge_min_width)
        net.add_edge(
            edge[0],
            edge[1],
            width=scaled_width,
            title=f"Edge Weight: {weight}",
            color="rgba(100, 200, 150, 0.8)"
        )

    # Step 7: Disable physics for fixed layout
    net.toggle_physics(False)

    # Step 8: Generate the HTML file
    output_html = f"{graph_name}_top{top_n}.html"
    net.write_html(output_html, notebook=False)

    # Step 9: Add explanation text
    explanation_text = f"""
    <div style="padding: 10px; font-family: Arial, sans-serif; background-color: #f9f9f9; border: 1px solid #ddd; margin-bottom: 20px;">
        <h2>Top {top_n} Nodes by Sequence Length (Circular Layout)</h2>
        <p>
            This visualization represents the <b>top {top_n}</b> largest nodes (by sequence length) and their edges in the consensus path.
            Nodes are scaled up 4x for visibility and colored based on their sequence length ranges.
            Edge thickness increases with edge weight for better representation.
        </p>
        <p>
            Hover over a node to see detailed information, including:
            <ul>
                <li>Node ID</li>
                <li>Sequence Length</li>
                <li>Scaled Size</li>
                <li>Node Count</li>
            </ul>
        </p>
        <p>
            Use interactive features to zoom, pan, and explore the graph.
        </p>
    </div>
    """
    with open(output_html, "r") as file:
        html_content = file.read()

    updated_html = html_content.replace("<body>", f"<body>{explanation_text}")
    with open(output_html, "w") as file:
        file.write(updated_html)

    print(f"Consensus path visualization saved to {output_html}")


def process_paf_file_with_reference(paf_file, reference_name, sequences_dict, output_fasta):
    print("Processing PAF file with reference...")
    alignments = defaultdict(lambda: '')
    with open(paf_file, 'r') as paf:
        for line in paf:
            paf_data = parse_paf_line(line)
            if paf_data['query_seq_name'] == reference_name:
                ref_seq = sequences_dict[paf_data['query_seq_name']][paf_data['query_start']:paf_data['query_end']]
                query_seq = sequences_dict[paf_data['target_seq_name']][paf_data['target_start']:paf_data['target_end']]
                if paf_data['strand'] == '-':
                    query_seq = query_seq[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))
                aligned_ref, aligned_query = restore_alignment(paf_data['cigar'], ref_seq, query_seq)
                alignments[paf_data['query_seq_name']] += aligned_ref
                alignments[paf_data['target_seq_name']] += aligned_query
            elif paf_data['target_seq_name'] == reference_name:
                ref_seq = sequences_dict[paf_data['target_seq_name']][paf_data['target_start']:paf_data['target_end']]
                query_seq = sequences_dict[paf_data['query_seq_name']][paf_data['query_start']:paf_data['query_end']]
                if paf_data['strand'] == '-':
                    query_seq = query_seq[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))
                aligned_ref, aligned_query = restore_alignment(paf_data['cigar'], ref_seq, query_seq)
                alignments[paf_data['target_seq_name']] += aligned_ref
                alignments[paf_data['query_seq_name']] += aligned_query
    alignments = truncate_alignment_sequences(alignments)
    with open(output_fasta, 'w') as f:
        for seq_name, seq in alignments.items():
            f.write(f">{seq_name}\n{seq}\n")


def calculate_graph_density(num_nodes, num_edges):
    """
    Calculate the density of a graph.

    Args:
        num_nodes (int): Number of nodes in the graph.
        num_edges (int): Number of edges in the graph.

    Returns:
        float: Graph density.
    """
    if num_nodes < 2:
        return 0  # A graph with fewer than 2 nodes has no edges

    max_edges = (num_nodes * (num_nodes - 1)) / 2  # Maximum possible edges
    density = num_edges / max_edges  # Density calculation
    return round(density * 1000, 12)


def process_maf_to_best_msa(input_maf, output_msa, target_length=100000):
    """
    Process a MAF file and generate the best MSA with fixed-length lines.

    Args:
        input_maf (str): Path to the input MAF file.
        output_msa (str): Path to save the MSA.
        target_length (int): Target length for the final MSA region.

    Returns:
        None
    """
    print("Processing MAF file...")
    blocks = []
    current_block = {}
    block_start = None

    # Read and parse MAF file
    with open(input_maf, "r") as maf_file:
        for line in maf_file:
            line = line.strip()
            if line.startswith("a"):
                # Save the current block
                if current_block:
                    blocks.append((block_start, current_block))
                current_block = {}
                block_start = None
            elif line.startswith("s "):
                parts = line.split()
                species = parts[1].split("#")[0]
                start = int(parts[2])
                size = int(parts[3])
                sequence = parts[6]

                end = start + size
                if block_start is None or start < block_start:
                    block_start = start

                current_block[species] = {"start": start, "end": end, "sequence": sequence}

        # Save the final block
        if current_block:
            blocks.append((block_start, current_block))

    # Sort blocks by starting coordinate
    blocks.sort(key=lambda x: x[0])

    print("Selecting best segments...")
    # Select the best blocks to cover the target length
    selected_blocks = []
    total_length = 0

    for block_start, block_sequences in blocks:
        # Determine block range
        block_length = max(seq["end"] for seq in block_sequences.values()) - block_start

        if total_length + block_length > target_length:
            # Truncate the block to fit the remaining length
            remaining_length = target_length - total_length
            for species in block_sequences:
                seq_data = block_sequences[species]
                seq_data["sequence"] = seq_data["sequence"][:remaining_length]
                seq_data["end"] = seq_data["start"] + remaining_length

        selected_blocks.append(block_sequences)
        total_length += block_length

        if total_length >= target_length:
            break

    print("Constructing fixed-length MSA...")
    # Generate MSA from selected blocks
    sequences = defaultdict(str)
    for block_sequences in selected_blocks:
        block_length = max(len(seq["sequence"]) for seq in block_sequences.values())
        for species, seq_data in block_sequences.items():
            sequences[species] += seq_data["sequence"].ljust(block_length, "-")
        # Add gaps for species not in the current block
        missing_species = set(sequences.keys()) - set(block_sequences.keys())
        for species in missing_species:
            sequences[species] += "-" * block_length

    # Ensure all sequences are of exact length
    for species in sequences:
        sequences[species] = sequences[species][:target_length].ljust(target_length, "-")

    print("Writing MSA...")
    # Save the MSA to file
    with open(output_msa, "w") as output_file:
        for species, sequence in sequences.items():
            output_file.write(f">{species}\n{sequence}\n")

    print(f"Fixed-length MSA saved to {output_msa}.")

def parse_maf(maf_file):
    """
    Parse a MAF file into blocks.

    Args:
        maf_file (str): Path to the input MAF file.

    Returns:
        list: List of MAF blocks.
    """
    blocks = []
    current_block = []
    with open(maf_file, "r") as file:
        for line in file:
            if line.startswith("a"):  # Start of a new block
                if current_block:
                    blocks.append(current_block)
                current_block = [line.strip()]
            elif line.startswith("s"):  # Sequence line
                current_block.append(line.strip())
        if current_block:
            blocks.append(current_block)
    return blocks


def sort_maf_blocks(blocks):
    """
    Sort MAF blocks by the starting position of the first sequence in each block.

    Args:
        blocks (list): List of MAF blocks.

    Returns:
        list: Sorted list of MAF blocks.
    """
    def get_start_pos(block):
        for line in block[1:]:
            if line.startswith("s"):
                return int(line.split()[2])  # Extract start position
        return float("inf")  # If no sequence line, place at end

    return sorted(blocks, key=get_start_pos)


def generate_msa(sorted_blocks, output_file):
    """
    Generate a multiple sequence alignment (MSA) from sorted MAF blocks.

    Args:
        sorted_blocks (list): List of sorted MAF blocks.
        output_file (str): Path to save the output MSA file.
    """
    # Collect all unique species names
    species = set()
    for block in sorted_blocks:
        for seq in block[1:]:
            species.add(seq.split()[1].split("#")[0])  # Extract species name

    species = sorted(species)  # Consistent ordering of species
    msa_data = defaultdict(list)

    # Create MSA from blocks
    for block in sorted_blocks:
        block_sequences = {seq.split()[1].split("#")[0]: seq.split()[-1] for seq in block[1:]}
        max_len = max(len(seq) for seq in block_sequences.values())

        for sp in species:
            if sp in block_sequences:
                msa_data[sp].append(block_sequences[sp])
            else:
                msa_data[sp].append('-' * max_len)  # Fill with gaps

    # Concatenate sequences for each species
    concatenated_msa = {sp: ''.join(msa_data[sp]) for sp in species}

    # Write to output file in FASTA format
    with open(output_file, 'w') as out:
        for sp, seq in concatenated_msa.items():
            out.write(f">{sp}\n{seq}\n")

    print(f"MSA saved to {output_file}")


def maf_to_msa(maf_file, output_file):
    """
    Wrapper function to convert a MAF file into an MSA.

    Args:
        maf_file (str): Path to the input MAF file.
        output_file (str): Path to save the output MSA file.
    """

    blocks = parse_maf(maf_file)
    sorted_blocks = sort_maf_blocks(blocks)
    generate_msa(sorted_blocks, output_file)


def plot_graph_statistics(gfa_file, output_dir, w_size):
    """
    Generate and plot graphs based on the provided GFA file using Gretl and specified Python scripts.

    Args:
        gfa_file (str): Path to the input GFA file.
        output_dir (str): Directory to save the output files.

    Returns:
        None
        :param threads: threads to use
    """
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Define output file paths
    similarity_output = os.path.join(output_dir, "similarity")
    path_output = os.path.join(output_dir, "paths")

    # 2. On a genomic (linear) scale, which regions are similar to each other? (window + window.py)
    try:
        print("Generating similarity plot...")
        subprocess.run(
            f"gretl window -s {w_size} -w {w_size} -g {gfa_file} -o - | python3 window.py -i - -o {similarity_output}",
            shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        print(f"Similarity plot saved to {similarity_output}")
    except subprocess.CalledProcessError as e:
        print(f"Error generating similarity plot: {e}")
        exit(1)

    try:
        print("Generating path statistics plot...")
        subprocess.run(
            f"gretl stats -p -g {gfa_file} -o - | python3 path.py -i - -o {path_output}",
            shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        print(f"Path statistics plot saved to {path_output}")
    except subprocess.CalledProcessError as e:
        print(f"Error generating path statistics plot: {e}")
        exit(1)

    print("All plots generated successfully.")



def extract_graph_stats(
    node_count_hist,
    edges_count,
    nodes_num,
    edges_num,
    avg_degree,
    num_haplotypes,
    name,
    vcf_file,
    c_dist,
    softcore_threshold: int = 5
):
    """
    Extract and summarize graph statistics using 3 node categories:
    - Unique: nodes present in ≤ softcore_threshold % of haplotypes
    - Shared: nodes present in between unique/core cutoffs
    - Core:   nodes present in ≥ (100 - softcore_threshold) % of haplotypes

    Also exports results as <name>_graph_stats.csv.
    """
    import math
    import pandas as pd
    import os

    low_cut = math.ceil((softcore_threshold / 100) * num_haplotypes)
    high_cut = math.floor(((100 - softcore_threshold) / 100) * num_haplotypes)

    print(f"Classifying nodes using threshold {softcore_threshold}%:")
    print(f"  Unique: ≤ {low_cut} haplotypes")
    print(f"  Shared: > {low_cut} and < {high_cut}")
    print(f"  Core: ≥ {high_cut} haplotypes")

    # --- Classification ---
    unique_nodes = sum(v for k, v in c_dist.items() if k <= low_cut)
    shared_nodes = sum(v for k, v in c_dist.items() if low_cut < k < high_cut)
    core_nodes   = sum(v for k, v in c_dist.items() if k >= high_cut)

    graph_stats = {
        "Statistic": [
            "Total number of nodes",
            f"Unique nodes (≤ {softcore_threshold}%)",
            f"Shared nodes ({softcore_threshold}%–{100 - softcore_threshold}%)",
            f"Core nodes (≥ {100 - softcore_threshold}%)",
            "Total number of edges",
            "Average node degree",
            "Graph density",
            "Number of haplotypes (paths)",
            "Average nucleotide diversity (π)"
        ],
        "Value": [
            nodes_num,
            unique_nodes,
            shared_nodes,
            core_nodes,
            edges_num,
            f"{avg_degree:.2f}",
            f"{calculate_graph_density(nodes_num, edges_num):.12f}",
            num_haplotypes,
            calculate_genomewide_pi('.'.join(vcf_file.split('.')[:-1]))
        ],
        "Explanation": [
            "Total number of sequence nodes in the variation graph.",
            f"Nodes present in ≤ {softcore_threshold}% of haplotypes (private/rare regions).",
            f"Nodes present in > {softcore_threshold}% and < {100 - softcore_threshold}% of haplotypes (variable regions).",
            f"Nodes present in ≥ {100 - softcore_threshold}% of haplotypes (highly conserved core).",
            "Total number of graph edges connecting nodes.",
            "Average number of connections per node.",
            "Proportion of observed vs. possible connections.",
            "Number of haplotype paths represented in the graph.",
            "Genome-wide average nucleotide diversity across paths."
        ]
    }

    # --- Export to CSV ---
    df = pd.DataFrame(graph_stats)

    # Use directory of `name`, not CWD
    base = os.path.splitext(os.path.basename(name))[0]  # "graph_stats"
    output_file = os.path.join(os.path.dirname(name), f"{base}_graph_stats.csv")

    df.to_csv(output_file, index=False)
    print(f"Graph statistics successfully exported to: {output_file}")


def calculate_genomewide_pi(vcf_file):
    vcf = VCF(vcf_file)

    pi_values = []

    for variant in vcf:
        if not variant.is_snp or len(variant.ALT) != 1: #skip non-biallelic SNPs
            continue

        genotypes = variant.genotypes
        alleles = [allele for g in genotypes for allele in g[:2] if allele != -1]

        if not alleles:
            continue

        ref_count = alleles.count(0)
        alt_count = alleles.count(1)
        total = ref_count + alt_count

        if total == 0:
            continue

        p = alt_count / total
        pi = 2 * p * (1 - p)
        pi_values.append(pi)

    genome_avg_pi = sum(pi_values) / len(pi_values) if pi_values else 0
    return genome_avg_pi

def write_odgi_matrix(gfa_file, msa_out_dir):
    """
    a fucntion to write a  |V|x|Paths| binary present/absent segments matrix
    :param gfa_file: gfa file directory
    :param msa_out_dir: output tsv file directory
    :return:
    """
    odgi_cmd = ['odgi', 'paths', '-i', gfa_file, '-H']
    try:
        odgi_res = subprocess.run(odgi_cmd, capture_output=True, text=True)
        if odgi_res.returncode != 0:
            raise OSError('odgi failed')
        with open(msa_out_dir + '/odgi_matrix.tsv', 'w') as fout:
            fout.write(odgi_res.stdout)
    except Exception:
        print("error writing odgi_matrix.tsv")
        raise Exception




from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage

def _linkage_to_newick(Z, labels):
    """
    Convert a SciPy linkage matrix to Newick format.
    """
    from collections import defaultdict
    n = len(labels)
    children = defaultdict(list)
    for i, (c1, c2, dist, _) in enumerate(Z, start=n):
        children[i] = [int(c1), int(c2)]

    def build_newick(node):
        if node < n:
            name = labels[node].replace("(", "_").replace(")", "_").replace(":", "_").replace(",", "_")
            return f"{name}"
        left, right = children[node]
        return f"({build_newick(left)},{build_newick(right)})"

    root = len(Z) + n - 1
    return build_newick(root) + ";"

def _read_odgi_matrix(odgi_matrix_tsv):
    """
    Reads an ODGI matrix TSV (from write_odgi_matrix) and returns paths and data matrix.
    """
    df = pd.read_csv(odgi_matrix_tsv, sep="\t")
    path_col = "path.name" if "path.name" in df.columns else df.columns[0]
    paths = df[path_col].astype(str).tolist()
    X = df.drop(columns=[path_col])
    # Keep only variable columns
    X = X.loc[:, X.var(axis=0, numeric_only=True) > 0]
    if X.empty:
        raise ValueError("ODGI matrix has no variable columns.")
    return paths, X.values

def generate_graph_based_tree(
    odgi_matrix_tsv,
    output_prefix="graph_phylo_tree",
    n_pcs=20,
    standardize=True,
    random_state=42
):
    """
    Generates a phylogenetic tree from an ODGI path-node matrix.

    Steps:
      1. Load ODGI matrix (paths × nodes)
      2. Filter constant columns
      3. PCA (dim reduction)
      4. Ward hierarchical clustering
      5. Save Newick tree

    Args:
        odgi_matrix_tsv (str): Path to odgi_matrix.tsv
        output_prefix (str): Output prefix for the Newick file
        n_pcs (int): Number of principal components to use
        standardize (bool): Whether to standardize columns before PCA
        random_state (int): Random seed for reproducibility
    """
    print(f"[graph-tree] Reading ODGI matrix: {odgi_matrix_tsv}")
    paths, X = _read_odgi_matrix(odgi_matrix_tsv)

    # Standardize and run PCA
    if standardize:
        X = StandardScaler(with_mean=True, with_std=True).fit_transform(X)

    k = max(2, min(n_pcs, X.shape[1], X.shape[0] - 1))
    print(f"[graph-tree] Running PCA (n_components={k}) on {X.shape[0]} paths × {X.shape[1]} nodes...")
    Xp = PCA(n_components=k, random_state=random_state).fit_transform(X)

    # Ward.D2 hierarchical clustering
    print("[graph-tree] Building hierarchical clustering (Ward)...")
    Z = linkage(Xp, method="ward", metric="euclidean")

    # Convert to Newick format
    print("[graph-tree] Writing Newick tree...")
    newick = _linkage_to_newick(Z, labels=paths)
    tree_file = f"{output_prefix}.tree"
    with open(tree_file, "w") as f:
        f.write(newick)
    print(f"[graph-tree] ✅ Saved: {tree_file}")

    return {"tree_file": tree_file}
