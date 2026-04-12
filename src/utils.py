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
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
import math
from matplotlib.patches import Patch
from collections import defaultdict
from Bio.Seq import Seq
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator




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



import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

def plot_node_histogram_by_paths(
    node_histogram,
    num_paths,
    output_file="node_histogram_by_paths.png",
    softcore_threshold: int = 5,
    max_bars: int = 120      # controls automatic bin merging
):
    """
    Plot node presence distribution colored by category:
      - Cloud (≤ threshold%)
      - Shell (middle range)
      - Core (≥ 100−threshold%)

    Improvements:
      • X-axis always integer (haplotype counts)
      • Dynamic bin merging for large num_paths
      • Sparse, readable x-axis ticks (no clutter)
      • Automatic log-scale for wide Y ranges
    """

    # --------------------------------------------------
    # 1. Determine cutoffs (in path-count units)
    # --------------------------------------------------
    low_cut  = math.ceil((softcore_threshold / 100) * num_paths)
    high_cut = math.floor(((100 - softcore_threshold) / 100) * num_paths)

    # --------------------------------------------------
    # 2. Prepare raw data
    # --------------------------------------------------
    x_raw = sorted(int(k) for k in node_histogram.keys())
    y_raw = np.array([node_histogram[c] for c in x_raw], dtype=float)

    n_raw = len(x_raw)

    # --------------------------------------------------
    # 3. Decide bin width dynamically
    # --------------------------------------------------
    if n_raw > max_bars:
        bin_width = math.ceil(n_raw / max_bars)
    else:
        bin_width = 1

    # --------------------------------------------------
    # 4. Merge bins (adjacent observed values)
    # --------------------------------------------------
    x_plot, y_plot = [], []

    i = 0
    while i < n_raw:
        j = min(i + bin_width, n_raw)
        xs = x_raw[i:j]
        ys = y_raw[i:j]

        x_plot.append((xs[0] + xs[-1]) / 2.0)
        y_plot.append(np.sum(ys))

        i = j

    x_plot = np.array(x_plot)
    y_plot = np.array(y_plot)

    # --------------------------------------------------
    # 5. Assign colors by bin center
    # --------------------------------------------------
    colors = []
    for c in x_plot:
        if c <= low_cut:
            colors.append("#1f77b4")   # Cloud
        elif c >= high_cut:
            colors.append("#2ca02c")   # Core
        else:
            colors.append("#ff7f0e")   # Shell

    # --------------------------------------------------
    # 6. Decide log scale
    # --------------------------------------------------
    positive = y_plot[y_plot > 0]
    use_log = (
        len(positive) > 0
        and np.max(positive) / max(np.min(positive), 1) > 1e3
    )

    # --------------------------------------------------
    # 7. Plot
    # --------------------------------------------------
    plt.figure(figsize=(12, 4))

    bar_width = 0.9 * (
        np.median(np.diff(x_plot)) if len(x_plot) > 1 else 1.0
    )

    plt.bar(
        x_plot,
        y_plot,
        width=bar_width,
        color=colors,
        edgecolor="none"
    )

    ax = plt.gca()

    plt.xlabel("Number of Paths (Haplotypes)")
    plt.ylabel("Number of Nodes")

    if use_log:
        plt.yscale("log")

    # --------------------------------------------------
    # 8. CLEAN X-AXIS: integer ticks, sparse
    # --------------------------------------------------
    ax.xaxis.set_major_locator(
        MaxNLocator(integer=True, nbins=10)
    )
    ax.set_xlim(0.5, num_paths + 0.5)

    # --------------------------------------------------
    # 9. Legend
    # --------------------------------------------------
    legend_elements = [
        Patch(facecolor="#1f77b4", edgecolor="black",
              label=f"Cloud (≤ {softcore_threshold}%)"),
        Patch(facecolor="#ff7f0e", edgecolor="black",
              label=f"Shell ({softcore_threshold}%–{100 - softcore_threshold}%)"),
        Patch(facecolor="#2ca02c", edgecolor="black",
              label=f"Core (≥ {100 - softcore_threshold}%)"),
    ]

    plt.legend(
        handles=legend_elements,
        loc="upper right",
        frameon=False
    )

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(
        f"Saved {output_file} | "
        f"bins={len(x_plot)}, bin_width={bin_width}, logY={use_log}"
    )





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







def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())




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
    print(f"  Cloud: ≤ {low_cut} haplotypes")
    print(f"  Shell: > {low_cut} and < {high_cut}")
    print(f"  Core: ≥ {high_cut} haplotypes")

    # --- Classification ---
    unique_nodes = sum(v for k, v in c_dist.items() if k <= low_cut)
    shared_nodes = sum(v for k, v in c_dist.items() if low_cut < k < high_cut)
    core_nodes   = sum(v for k, v in c_dist.items() if k >= high_cut)

    graph_stats = {
        "Statistic": [
            "Total number of nodes",
            f"Cloud nodes (≤ {softcore_threshold}%)",
            f"Shell nodes ({softcore_threshold}%–{100 - softcore_threshold}%)",
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



def _linkage_to_newick(Z, labels):
    """
    Convert a SciPy linkage matrix to Newick format with branch lengths.
    Branch length for each child = parent merge height - child merge height
    (leaf height = 0).
    """
    from collections import defaultdict
    n = len(labels)
    children = defaultdict(list)
    heights = {}

    # Leaf heights are 0
    for i in range(n):
        heights[i] = 0.0

    for i, (c1, c2, dist, _) in enumerate(Z, start=n):
        children[i] = [int(c1), int(c2)]
        heights[i] = float(dist)

    def sanitize(name):
        for ch in ("(", ")", ":", ",", "#", ";"):
            name = name.replace(ch, "_")
        return name

    def build_newick(node, parent_height):
        branch_len = parent_height - heights[node]
        if node < n:
            name = sanitize(labels[node])
            return f"{name}:{branch_len:.6f}"
        left, right = children[node]
        node_height = heights[node]
        return (
            f"({build_newick(left, node_height)},"
            f"{build_newick(right, node_height)}):{branch_len:.6f}"
        )

    root = len(Z) + n - 1
    root_height = heights[root]
    left, right = children[root]
    return f"({build_newick(left, root_height)},{build_newick(right, root_height)});"

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


def read_odgi_as_presence_sets(tsv_path):
    paths = []
    node_sets = []  # list of Python sets of node indices

    # if gzipped, handle that here (gzip.open)
    with open(tsv_path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        # header[0] is probably something like "path", rest are nodes
        n_nodes = len(header) - 1

        print(f"[graph-tree] Header read: {n_nodes} nodes")

        for line_idx, line in enumerate(f, start=1):
            parts = line.rstrip("\n").split("\t")
            path_name = parts[0]
            # Skip path.length (col 1) and path.step.count (col 2); nodes start at col 3
            values = parts[3:]

            # Build set of node indices where value != "0"
            present = {j for j, v in enumerate(values) if v != "0" and v != ""}

            paths.append(path_name)
            node_sets.append(present)

            if line_idx % 5 == 0:
                print(f"[graph-tree] Processed {line_idx} paths")

    return paths, node_sets, n_nodes

def compute_jaccard_distance_matrix(node_sets):
    n = len(node_sets)
    D = np.zeros((n, n), dtype=float)

    sizes = [len(s) for s in node_sets]

    for i in range(n):
        for j in range(i + 1, n):
            inter = len(node_sets[i] & node_sets[j])
            union = sizes[i] + sizes[j] - inter
            if union == 0:
                d = 0.0
            else:
                d = 1.0 - inter / union
            D[i, j] = D[j, i] = d

    return D


def read_odgi_as_binary_matrix(tsv_path):
    """
    Read an ODGI presence/absence matrix TSV and return a dense binary numpy array.

    The TSV format produced by ``odgi paths -H`` has the columns:
        path.name | path.length | path.step.count | node_0 | node_1 | ...

    This function skips the three metadata columns and returns a
    (n_paths × n_nodes) boolean matrix where True = node present in path.

    Returns
    -------
    paths : list[str]
    binary_matrix : np.ndarray, shape (n_paths, n_nodes), dtype bool
    """
    paths = []
    rows = []
    with open(tsv_path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # Columns: path.name(0) | path.length(1) | path.step.count(2) | nodes(3…)
        n_nodes = max(len(header) - 3, 0)

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            paths.append(parts[0])
            values = parts[3:]  # skip metadata columns
            row = np.array(
                [v != "0" and v != "" for v in values[:n_nodes]],
                dtype=bool,
            )
            # Pad with False if a row is shorter than expected
            if len(row) < n_nodes:
                row = np.concatenate([row, np.zeros(n_nodes - len(row), dtype=bool)])
            rows.append(row)

    if rows:
        binary_matrix = np.vstack(rows)
    else:
        binary_matrix = np.zeros((0, n_nodes), dtype=bool)
    return paths, binary_matrix


def _filter_informative_nodes(binary_matrix):
    """
    Remove phylogenetically uninformative columns from a PAV binary matrix.

    Columns removed:
    - All-zero  : node absent in every path      (no information)
    - Singleton : node present in exactly 1 path  (autapomorphic; cannot support any split)
    - All-one   : node present in every path      (invariant core node; no variation)

    A column is *informative* only when ``1 < col_sum < n_paths``, i.e. it
    distinguishes at least two paths from at least one other.

    Returns
    -------
    filtered_matrix : np.ndarray, shape (n_paths, n_informative)
    kept_mask       : np.ndarray[bool], shape (n_nodes,)  — True for kept columns
    stats           : dict with counts for logging
    """
    n_paths = binary_matrix.shape[0]
    col_sums = binary_matrix.sum(axis=0)

    kept_mask = (col_sums > 1) & (col_sums < n_paths)

    n_total      = binary_matrix.shape[1]
    n_all_zero   = int((col_sums == 0).sum())
    n_singleton  = int((col_sums == 1).sum())
    n_universal  = int((col_sums == n_paths).sum())
    n_informative = int(kept_mask.sum())

    stats = {
        "total_nodes":      n_total,
        "removed_all_zero": n_all_zero,
        "removed_singleton": n_singleton,
        "removed_universal": n_universal,
        "informative_nodes": n_informative,
    }
    return binary_matrix[:, kept_mask], kept_mask, stats


def _get_bipartitions(Z, n_leaves):
    """
    Extract canonical bipartitions from a SciPy linkage matrix.

    Each split is stored as the SMALLER half only (canonical form), as a plain
    frozenset of integer leaf indices.  Storing one side instead of both halves
    halves hashing cost and removes the need for a frozenset-of-frozensets.

    Returns
    -------
    set[frozenset[int]]
    """
    all_leaves = frozenset(range(n_leaves))
    subtree    = {i: frozenset([i]) for i in range(n_leaves)}
    half       = n_leaves / 2
    bipartitions = set()

    for idx, (c1, c2, _, _) in enumerate(Z, start=n_leaves):
        merged       = subtree[int(c1)] | subtree[int(c2)]
        subtree[idx] = merged
        canonical    = merged if len(merged) <= half else (all_leaves - merged)
        bipartitions.add(canonical)

    return bipartitions


def _bootstrap_worker_mmap(rep, mat_path, idx_path, n_paths, half, all_leaves):
    """
    Joblib/loky worker for one Felsenstein replicate.

    Loads the PAV matrix and column-index array from memory-mapped numpy files
    written to /tmp (tmpfs = RAM inside Docker) — zero-copy read via OS page
    cache, no large-array pickling, no fork() lock-inheritance issues.
    """
    _mat = np.load(mat_path, mmap_mode='r')
    _col = np.load(idx_path, mmap_mode='r')[rep]
    return _one_bootstrap_replicate(_mat, _col, n_paths, half, all_leaves)


def _one_bootstrap_replicate(mat, col_idx, n_paths, half, all_leaves):
    """
    Run a single bootstrap replicate.

    Defined at module level so it is picklable (required by ProcessPoolExecutor
    if ever switched from threads to processes).  All arguments are read-only;
    ``mat`` is shared safely across threads because scipy C code releases the GIL.
    """
    from scipy.spatial.distance import pdist as _pdist
    from scipy.cluster.hierarchy import linkage as _linkage

    boot_mat  = mat[:, col_idx]
    condensed = _pdist(boot_mat, metric="jaccard")
    np.nan_to_num(condensed, copy=False, nan=1.0)
    Z_boot    = _linkage(condensed, method="average")

    subtree = {i: frozenset([i]) for i in range(n_paths)}
    biparts = set()
    for idx, (c1, c2, _, _) in enumerate(Z_boot, start=n_paths):
        merged       = subtree[int(c1)] | subtree[int(c2)]
        subtree[idx] = merged
        biparts.add(merged if len(merged) <= half else (all_leaves - merged))
    return biparts


def _bootstrap_upgma(filtered_matrix, n_bootstrap, seed=42, n_threads=1, tmp_dir=None):
    """
    Fast Felsenstein bootstrap on the PAV binary matrix for UPGMA trees.

    Speed-ups over the naive Python-loop implementation:

    1. **Vectorised Jaccard** — ``scipy.spatial.distance.pdist`` with
       ``metric='jaccard'`` is C-backed and returns the condensed distance
       vector directly, bypassing the O(n² × k) Python set-intersection loop
       and the subsequent ``squareform`` round-trip.  Typical speed-up: ×50–200.

    2. **Parallel replicates** — when ``n_threads > 1``, replicates are distributed
       across ``joblib.Parallel`` with the ``loky`` backend.  ``loky`` uses a
       forkserver-style clean worker pool that avoids the BLAS/threading lock
       inheritance that makes ``fork`` silently deadlock inside Docker.  The PAV
       matrix and column-index array are written once to ``/tmp`` (tmpfs = RAM
       inside the container) as ``.npy`` files; each worker loads them with
       ``mmap_mode='r'``, so the OS page cache shares the pages — effectively
       zero-copy across all workers.  Only a small integer ``rep`` is pickled
       per task.  Speed-up is approximately O(n_bootstrap / n_threads).

    3. **Bulk index generation** — all n_bootstrap column-resample arrays are
       generated in a single ``numpy.random.Generator.integers`` call before
       the loop.

    4. **Canonical bipartitions** — each split is represented as the smaller-half
       frozenset, halving hashing cost.

    Parameters
    ----------
    filtered_matrix : np.ndarray[bool], shape (n_paths, n_informative_nodes)
    n_bootstrap     : int
    seed            : int — RNG seed for reproducibility
    n_threads       : int — worker processes for parallel replicates (default: 1)

    Returns
    -------
    list[set[frozenset[int]]]  — one canonical bipartition set per replicate
    """
    import os, tempfile
    from joblib import Parallel, delayed

    rng              = np.random.default_rng(seed)
    n_paths, n_nodes = filtered_matrix.shape
    half             = n_paths / 2
    all_leaves       = frozenset(range(n_paths))

    # Generate all bootstrap column-index arrays upfront (single numpy call)
    all_col_idx = rng.integers(0, n_nodes, size=(n_bootstrap, n_nodes))

    # Convert bool → float64 once; scipy pdist('jaccard') requires numeric input.
    mat = filtered_matrix.astype(np.float64)

    n_workers = max(1, int(n_threads))

    if n_workers == 1:
        # Serial path — no process-pool overhead
        boot_bipartitions = []
        for rep in range(n_bootstrap):
            if (rep + 1) % 100 == 0:
                print(f"[graph-tree] Bootstrap replicate {rep + 1}/{n_bootstrap}...")
            boot_bipartitions.append(
                _one_bootstrap_replicate(mat, all_col_idx[rep], n_paths, half, all_leaves)
            )
    else:
        print(f"[graph-tree] Bootstrap: {n_bootstrap} replicates across {n_workers} workers "
              f"(joblib/loky)...")
        # Write arrays to /tmp (tmpfs = RAM inside Docker) as .npy files.
        # Workers load them with mmap_mode='r' — OS page cache shares the pages
        # across all worker processes, giving zero-copy read without fork().
        mat_tmp = tempfile.mktemp(suffix='_boot_mat.npy', dir=tmp_dir)
        idx_tmp = tempfile.mktemp(suffix='_boot_idx.npy', dir=tmp_dir)
        try:
            np.save(mat_tmp, mat)
            np.save(idx_tmp, all_col_idx)
            boot_bipartitions = Parallel(
                n_jobs=n_workers, backend='loky', verbose=5,
            )(
                delayed(_bootstrap_worker_mmap)(
                    rep, mat_tmp, idx_tmp, n_paths, half, all_leaves
                )
                for rep in range(n_bootstrap)
            )
        finally:
            for _f in (mat_tmp, idx_tmp):
                try:
                    os.unlink(_f)
                except OSError:
                    pass

    return boot_bipartitions


def _compute_bootstrap_support(Z, n_leaves, bootstrap_bipartitions):
    """
    Map bootstrap frequencies onto the internal nodes of the original UPGMA tree.

    Uses ``collections.Counter`` (C-backed hash table) to tally all bootstrap
    bipartitions in a single pass, then looks up each original-tree split in
    O(1).  This avoids the quadratic inner loop of the naive approach.

    Returns
    -------
    dict[int, int]  — node index → support percentage 0–100
    """
    from collections import Counter

    n_bootstrap = len(bootstrap_bipartitions)
    if n_bootstrap == 0:
        return {}

    # ── One-pass tally of all bootstrap bipartitions ─────────────────────────
    boot_counter = Counter()
    for bp_set in bootstrap_bipartitions:
        boot_counter.update(bp_set)   # batch-adds every key in the set

    # ── Look up support for each split in the original tree ──────────────────
    all_leaves   = frozenset(range(n_leaves))
    half         = n_leaves / 2
    subtree      = {i: frozenset([i]) for i in range(n_leaves)}
    support_values = {}

    for idx, (c1, c2, _, _) in enumerate(Z, start=n_leaves):
        merged       = subtree[int(c1)] | subtree[int(c2)]
        subtree[idx] = merged
        canonical    = merged if len(merged) <= half else (all_leaves - merged)
        support_values[idx] = int(100 * boot_counter.get(canonical, 0) / n_bootstrap)

    return support_values


# =============================================================================
# UFBoot-PAV — ultrafast adaptive bootstrap for UPGMA on binary PAV matrices
# =============================================================================

def _jaccard_condensed_rell(X_float, Xt_float, m):
    """
    Compute the condensed Jaccard distance vector directly from a column
    multiplicity vector, without creating an explicit resampled matrix.

    **Mathematical equivalence with Felsenstein bootstrap**

    Standard column-resampling bootstrap builds boot_mat = X[:, col_idx]
    (with repeated columns) and calls pdist(boot_mat, 'jaccard').
    This function computes exactly the same distances via weighted inner products:

        inter[i,j]  = Σ_c  m[c] · X[i,c] · X[j,c]  =  (X · diag(m)) @ X.T
        count_i     = Σ_c  m[c] · X[i,c]             =  X @ m
        union[i,j]  = count_i + count_j − inter[i,j]
        d[i,j]      = 1 − inter[i,j] / max(union[i,j], ε)

    **Why it is faster than pdist on the resampled matrix**

    • No explicit resampled matrix allocation (avoids n_paths × n_nodes copy).
    • The inner-product step ``(X * m) @ X.T`` is dispatched as BLAS DGEMM,
      which uses SIMD instructions and — when numpy is linked against MKL or
      OpenBLAS — multiple threads automatically, without any extra thread
      management on our side.
    • ``X @ m`` is BLAS DGEMV (highly vectorised).

    Parameters
    ----------
    X_float  : np.ndarray[float64], shape (n_paths, n_nodes)
    Xt_float : np.ndarray[float64], shape (n_nodes, n_paths) — contiguous X.T
    m        : np.ndarray[float64], shape (n_nodes,)  — column multiplicity vector

    Returns
    -------
    np.ndarray[float64], condensed distance vector, length n_paths*(n_paths-1)//2
    """
    inter = (X_float * m) @ Xt_float        # BLAS DGEMM — (n_paths × n_nodes) × (n_nodes × n_paths)
    count = X_float @ m                      # BLAS DGEMV — (n_paths,)
    union = count[:, None] + count[None, :] - inter
    np.maximum(union, 1e-10, out=union)      # in-place: avoid divide-by-zero
    D = 1.0 - inter / union
    # Diagonal is analytically 0 (inter[i,i]==count[i], union[i,i]==count[i])
    # but floating-point may introduce tiny residuals — zero it explicitly
    np.fill_diagonal(D, 0.0)
    n = D.shape[0]
    i_idx, j_idx = np.triu_indices(n, k=1)
    return D[i_idx, j_idx]


def _ufboot_worker_mmap(m, x_path, xt_path, n_paths, half, all_leaves):
    """
    Joblib/loky worker for one UFBoot-PAV replicate.

    Loads X and Xt from memory-mapped .npy files (written to /tmp = tmpfs).
    Only the multiplicity vector ``m`` (one float64 row, ~400 KB) is pickled
    per task — the large matrices are shared via OS page cache.
    """
    from scipy.cluster.hierarchy import linkage as _lnk
    _X  = np.load(x_path,  mmap_mode='r')
    _Xt = np.load(xt_path, mmap_mode='r')
    condensed = _jaccard_condensed_rell(_X, _Xt, m)
    Z_boot    = _lnk(condensed, method="average")
    subtree   = {i: frozenset([i]) for i in range(n_paths)}
    biparts   = set()
    for idx, (c1, c2, _, _) in enumerate(Z_boot, start=n_paths):
        merged       = subtree[int(c1)] | subtree[int(c2)]
        subtree[idx] = merged
        biparts.add(merged if len(merged) <= half else (all_leaves - merged))
    return biparts


def _ufboot_pav(
    filtered_matrix,
    n_threads=1,
    seed=42,
    min_replicates=1000,
    max_replicates=10000,
    batch_size=100,
    convergence_threshold=0.99,
    n_consecutive_converged=2,
    tmp_dir=None,
):
    """
    UFBoot-PAV: Ultrafast bootstrap with convergence-based stopping for UPGMA
    trees built on Jaccard distances from binary PAV matrices.

    Algorithmic design
    ------------------
    Adapts two core innovations from **UFBoot2** (Hoang et al. 2018, *MBE* 35:518)
    to the UPGMA / Jaccard-PAV setting:

    1. **RELL-weighted Jaccard** (replaces per-replicate pdist on resampled matrix)

       Each bootstrap replicate draws a multiplicity vector
       ``m ~ Multinomial(n_nodes, uniform)``.  The resampled Jaccard distances
       are computed analytically via BLAS matrix multiplication
       (see ``_jaccard_condensed_rell``), which is *mathematically identical*
       to building an explicit resampled matrix and calling pdist but avoids the
       column-copy and uses hardware-level SIMD + multi-threading.

    2. **Convergence-based stopping** (UFBoot2 criterion, Hoang et al. 2018)

       After every ``batch_size`` replicates the Pearson correlation between the
       current split-support vector and the previous one is computed.  When
       ``r ≥ convergence_threshold`` (default 0.99, same as UFBoot2) for
       ``n_consecutive_converged`` (default 2) consecutive batches AND at least
       ``min_replicates`` have been run, the algorithm stops.  This prevents
       running an arbitrary fixed number of replicates: fast-converging datasets
       finish early, slow-converging ones keep running up to ``max_replicates``.

    Relationship to standard Felsenstein bootstrap
    -----------------------------------------------
    The distances computed per replicate are **identical** to the standard
    Felsenstein column-resampling bootstrap.  UFBoot-PAV is a computational
    re-implementation, not a statistical approximation.  Support values carry
    the same frequentist interpretation.

    Parameters
    ----------
    filtered_matrix         : np.ndarray[bool], (n_paths × n_informative_nodes)
    n_threads               : int   — parallel workers (BLAS inside each replicate
                                      also benefits from numpy's thread pool)
    seed                    : int   — RNG seed
    min_replicates          : int   — minimum reps before convergence testing (1000)
    max_replicates          : int   — hard upper cap (10000)
    batch_size              : int   — reps per convergence-check round (100)
    convergence_threshold   : float — Pearson r threshold (0.99, UFBoot2 default)
    n_consecutive_converged : int   — consecutive rounds required (2, UFBoot2 default)

    Returns
    -------
    dict:
        ``Z``              — SciPy linkage matrix of the original UPGMA tree
        ``support_values`` — dict {internal_node_idx: support_pct 0-100}
        ``n_replicates``   — int, total replicates run
        ``converged``      — bool
    """
    import os, tempfile
    from joblib import Parallel, delayed
    from scipy.spatial.distance import pdist as _pdist
    from scipy.cluster.hierarchy import linkage as _linkage
    from collections import Counter

    rng      = np.random.default_rng(seed)
    n_paths, n_nodes = filtered_matrix.shape
    half     = n_paths / 2
    all_leaves = frozenset(range(n_paths))
    n_workers  = max(1, int(n_threads))

    # ── Pre-computation ───────────────────────────────────────────────────────
    X_float  = filtered_matrix.astype(np.float64)
    Xt_float = np.ascontiguousarray(X_float.T)

    # Original UPGMA tree
    condensed_orig = _pdist(X_float, metric="jaccard")
    np.nan_to_num(condensed_orig, copy=False, nan=1.0)
    Z = _linkage(condensed_orig, method="average")

    # Build ordered list of original bipartitions (stable ordering for Pearson r)
    subtree_tmp = {i: frozenset([i]) for i in range(n_paths)}
    orig_biparts = []
    for idx, (c1, c2, _, _) in enumerate(Z, start=n_paths):
        merged = subtree_tmp[int(c1)] | subtree_tmp[int(c2)]
        subtree_tmp[idx] = merged
        orig_biparts.append(merged if len(merged) <= half else (all_leaves - merged))

    # ── Write X and Xt to /tmp once (tmpfs = RAM); workers mmap them ──────────
    x_tmp  = tempfile.mktemp(suffix='_rell_x.npy',  dir=tmp_dir)
    xt_tmp = tempfile.mktemp(suffix='_rell_xt.npy', dir=tmp_dir)
    np.save(x_tmp,  X_float)
    np.save(xt_tmp, Xt_float)

    # ── Bootstrap loop ────────────────────────────────────────────────────────
    boot_counter     = Counter()
    total_reps       = 0
    converged        = False
    consecutive_ok   = 0
    prev_support_vec = None

    try:
        while total_reps < max_replicates:
            # ── Generate multiplicity vectors for this batch ──────────────────
            col_draws  = rng.integers(0, n_nodes, size=(batch_size, n_nodes))
            mult_vecs  = np.vstack([
                np.bincount(col_draws[r], minlength=n_nodes).astype(np.float64)
                for r in range(batch_size)
            ])                    # shape (batch_size, n_nodes)

            # ── Run batch ─────────────────────────────────────────────────────
            if n_workers > 1:
                batch_results = Parallel(n_jobs=n_workers, backend='loky')(
                    delayed(_ufboot_worker_mmap)(
                        m, x_tmp, xt_tmp, n_paths, half, all_leaves
                    )
                    for m in mult_vecs
                )
            else:
                batch_results = [
                    _ufboot_worker_mmap(m, x_tmp, xt_tmp, n_paths, half, all_leaves)
                    for m in mult_vecs
                ]

            for biparts in batch_results:
                boot_counter.update(biparts)
            total_reps += batch_size

            # ── Current support vector (same order as orig_biparts) ───────────
            curr_support_vec = np.array([
                boot_counter.get(bp, 0) / total_reps for bp in orig_biparts
            ])

            print(f"[ufboot-pav] {total_reps} replicates...")

            # ── UFBoot2 convergence criterion ──────────────────────────────────
            if total_reps >= min_replicates and prev_support_vec is not None:
                std_prev = prev_support_vec.std()
                std_curr = curr_support_vec.std()
                if std_prev < 1e-10 or std_curr < 1e-10:
                    # All support values identical → trivially converged
                    corr = 1.0
                else:
                    corr = float(np.corrcoef(prev_support_vec, curr_support_vec)[0, 1])
                    if np.isnan(corr):
                        corr = 1.0

                if corr >= convergence_threshold:
                    consecutive_ok += 1
                    print(f"[ufboot-pav] Convergence round {consecutive_ok}/{n_consecutive_converged}: "
                          f"r={corr:.4f}")
                    if consecutive_ok >= n_consecutive_converged:
                        print(f"[ufboot-pav] Converged after {total_reps} replicates.")
                        converged = True
                        break
                else:
                    if consecutive_ok > 0:
                        print(f"[ufboot-pav] Convergence reset (r={corr:.4f} < {convergence_threshold})")
                    consecutive_ok = 0

            prev_support_vec = curr_support_vec

        if not converged:
            print(f"[ufboot-pav] Max replicates ({max_replicates}) reached without convergence.")

    finally:
        for _f in (x_tmp, xt_tmp):
            try:
                os.unlink(_f)
            except OSError:
                pass

    # ── Map support percentages back to UPGMA linkage nodes ──────────────────
    bp_to_pct = {
        bp: int(100 * boot_counter.get(bp, 0) / total_reps)
        for bp in orig_biparts
    }
    support_values  = {}
    subtree_final   = {i: frozenset([i]) for i in range(n_paths)}
    for idx, (c1, c2, _, _) in enumerate(Z, start=n_paths):
        merged              = subtree_final[int(c1)] | subtree_final[int(c2)]
        subtree_final[idx]  = merged
        canonical           = merged if len(merged) <= half else (all_leaves - merged)
        support_values[idx] = bp_to_pct.get(canonical, 0)

    return {
        "Z":              Z,
        "support_values": support_values,
        "n_replicates":   total_reps,
        "converged":      converged,
    }


def _linkage_to_newick_with_support(Z, labels, support_values):
    """
    Convert a SciPy linkage matrix to Newick format with bootstrap support
    values on internal nodes and branch lengths.

    Internal node labels are written as ``{support}`` (integer 0–100) before
    the colon, following the standard convention used by RAxML and IQ-TREE.

    Parameters
    ----------
    Z              : linkage matrix from scipy.cluster.hierarchy.linkage
    labels         : list[str] — leaf names in the same order as rows of the distance matrix
    support_values : dict[int, int] — from _compute_bootstrap_support

    Returns
    -------
    str — Newick string
    """
    from collections import defaultdict
    n = len(labels)
    children = {}
    heights  = {i: 0.0 for i in range(n)}

    for idx, (c1, c2, dist, _) in enumerate(Z, start=n):
        children[idx] = [int(c1), int(c2)]
        heights[idx]  = float(dist)

    def sanitize(name):
        for ch in ("(", ")", ":", ",", "#", ";"):
            name = name.replace(ch, "_")
        return name

    def build_newick(node, parent_height):
        branch_len = parent_height - heights[node]
        if node < n:
            return f"{sanitize(labels[node])}:{branch_len:.6f}"
        left, right   = children[node]
        node_height   = heights[node]
        support       = support_values.get(node, "")
        return (
            f"({build_newick(left, node_height)},"
            f"{build_newick(right, node_height)})"
            f"{support}:{branch_len:.6f}"
        )

    root = len(Z) + n - 1
    root_height = heights[root]
    left, right = children[root]
    return f"({build_newick(left, root_height)},{build_newick(right, root_height)});"


def generate_graph_based_tree(
    odgi_matrix_tsv,
    output_prefix="graph_phylo_tree",
    n_bootstrap=100,
    n_threads=1,
    bootstrap_method="felsenstein",
    ufboot_min_rep=1000,
    ufboot_max_rep=10000,
    ufboot_convergence=0.99,
    ufboot_batch=100,
    tmp_dir=None,
):
    """
    Build a UPGMA phylogenetic tree from an ODGI PAV matrix with bootstrap support.

    Workflow
    --------
    1. Load the PAV matrix as a binary numpy array.
    2. Filter uninformative columns (all-absent, singletons, all-present).
    3. Compute Jaccard distances on the informative nodes.
    4. Cluster with UPGMA (average linkage).
    5. Run bootstrap (method selected by ``bootstrap_method``).
    6. Write ``{output_prefix}.tree`` (Newick with support) and a support TSV.

    Bootstrap methods
    -----------------
    ``"felsenstein"`` (default)
        Standard Felsenstein column-resampling bootstrap.  Runs exactly
        ``n_bootstrap`` replicates.  Each replicate uses ``scipy.spatial.distance.pdist``
        with the 'jaccard' metric (C-backed), parallelised across replicates via
        ``joblib.Parallel`` (loky backend, mmap matrix sharing via /tmp tmpfs).

    ``"ufboot"``
        UFBoot-PAV — ultrafast adaptive bootstrap (see ``_ufboot_pav``).
        Adapts two core ideas from UFBoot2 (Hoang et al. 2018, *MBE* 35:518):

        • **RELL-weighted Jaccard**: bootstrap distances are computed from a
          column-multiplicity vector via BLAS DGEMM (``(X * m) @ X.T``),
          which is *mathematically identical* to building an explicit resampled
          matrix but uses numpy's SIMD / multi-thread BLAS path instead of the
          pdist C loop.

        • **Convergence-based stopping**: replicates run in batches; the
          Pearson correlation between consecutive split-support vectors is
          tested after each batch.  Stops when r ≥ ``ufboot_convergence``
          (default 0.99, same as UFBoot2) for ``2`` consecutive batches AND
          at least ``ufboot_min_rep`` replicates have been done.  Avoids
          over-bootstrapping stable datasets and under-bootstrapping unstable ones.

        ``n_bootstrap`` has no effect when ``bootstrap_method='ufboot'``;
        use ``ufboot_min_rep`` / ``ufboot_max_rep`` instead.

    Parameters
    ----------
    odgi_matrix_tsv    : str   — ODGI TSV produced by write_odgi_matrix()
    output_prefix      : str   — path prefix for output files
    n_bootstrap        : int   — replicates for 'felsenstein' (0 = skip bootstrap)
    n_threads          : int   — worker threads
    bootstrap_method   : str   — 'felsenstein' | 'ufboot'
    ufboot_min_rep     : int   — UFBoot: min replicates before convergence testing (1000)
    ufboot_max_rep     : int   — UFBoot: hard cap on total replicates (10000)
    ufboot_convergence : float — UFBoot: Pearson r convergence threshold (0.99)
    ufboot_batch       : int   — UFBoot: replicates per convergence-check round (100)

    Returns
    -------
    dict with key ``tree_file`` pointing to the written Newick file.
    """
    from scipy.spatial.distance import pdist as _pdist

    print(f"[graph-tree] Reading ODGI matrix: {odgi_matrix_tsv}")
    paths, binary_matrix = read_odgi_as_binary_matrix(odgi_matrix_tsv)
    n_paths, n_nodes_total = binary_matrix.shape
    print(f"[graph-tree] Loaded {n_paths} paths × {n_nodes_total} nodes")

    # ── Step 1: Filter uninformative nodes ───────────────────────────────────
    filtered_matrix, _, filt_stats = _filter_informative_nodes(binary_matrix)
    print(
        f"[graph-tree] Node filtering: {filt_stats['total_nodes']} total → "
        f"{filt_stats['informative_nodes']} informative "
        f"(removed {filt_stats['removed_all_zero']} all-absent, "
        f"{filt_stats['removed_singleton']} singletons, "
        f"{filt_stats['removed_universal']} universal/core)"
    )
    if filtered_matrix.shape[1] == 0:
        raise ValueError(
            "[graph-tree] No phylogenetically informative nodes remain after filtering."
        )

    # ── Step 2 + 3: Original tree (always needed) ────────────────────────────
    print(
        f"[graph-tree] Computing Jaccard distances for {n_paths} paths × "
        f"{filtered_matrix.shape[1]} informative nodes..."
    )
    mat_orig  = filtered_matrix.astype(np.float64)
    cond_orig = _pdist(mat_orig, metric="jaccard")
    np.nan_to_num(cond_orig, copy=False, nan=1.0)
    print("[graph-tree] Building UPGMA tree (average linkage)...")
    Z = linkage(cond_orig, method="average")

    # ── Step 4: Bootstrap ────────────────────────────────────────────────────
    # Derive tmp_dir from output_prefix when not explicitly supplied so that
    # mmap temp files land on the already-mounted /output disk rather than /tmp.
    _tmp_dir = tmp_dir if tmp_dir is not None else os.path.dirname(output_prefix) or "."

    support_values  = {}
    n_reps_used     = 0
    converged_flag  = None   # None = not applicable (no ufboot)

    if bootstrap_method == "ufboot":
        print(
            f"[graph-tree] UFBoot-PAV: adaptive bootstrap "
            f"(min={ufboot_min_rep}, max={ufboot_max_rep}, "
            f"convergence r≥{ufboot_convergence}, batch={ufboot_batch}, "
            f"threads={n_threads})"
        )
        ufboot_result  = _ufboot_pav(
            filtered_matrix,
            n_threads=n_threads,
            min_replicates=ufboot_min_rep,
            max_replicates=ufboot_max_rep,
            batch_size=ufboot_batch,
            convergence_threshold=ufboot_convergence,
            tmp_dir=_tmp_dir,
        )
        support_values = ufboot_result["support_values"]
        n_reps_used    = ufboot_result["n_replicates"]
        converged_flag = ufboot_result["converged"]

    elif n_bootstrap > 0:
        print(
            f"[graph-tree] Felsenstein bootstrap: {n_bootstrap} replicates, "
            f"{n_threads} thread{'s' if n_threads != 1 else ''}..."
        )
        boot_biparts   = _bootstrap_upgma(filtered_matrix, n_bootstrap, n_threads=n_threads,
                                          tmp_dir=_tmp_dir)
        support_values = _compute_bootstrap_support(Z, n_paths, boot_biparts)
        n_reps_used    = n_bootstrap

    # ── Step 5: Write outputs ─────────────────────────────────────────────────
    if support_values:
        sup_vals = list(support_values.values())
        print(
            f"[graph-tree] Bootstrap done ({n_reps_used} replicates"
            + (f", converged={converged_flag}" if converged_flag is not None else "")
            + f"). Support range: {min(sup_vals)}–{max(sup_vals)}%"
        )
        support_tsv = f"{output_prefix}_bootstrap_support.tsv"
        with open(support_tsv, "w") as fh:
            fh.write("node_index\tsupport_pct\n")
            for node_idx, sup in sorted(support_values.items()):
                fh.write(f"{node_idx}\t{sup}\n")
        print(f"[graph-tree] Bootstrap support table: {support_tsv}")
        newick = _linkage_to_newick_with_support(Z, labels=paths,
                                                  support_values=support_values)
    else:
        newick = _linkage_to_newick(Z, labels=paths)

    tree_file = f"{output_prefix}.tree"
    with open(tree_file, "w") as fh:
        fh.write(newick)
    print(f"[graph-tree] Saved tree: {tree_file}")
    return {"tree_file": tree_file}