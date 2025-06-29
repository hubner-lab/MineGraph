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
    print("Generating data for nodes counts ...")
    my_hist = {}
    for key in paths.keys():
        for node in paths[key][0]:
            node_id = node[:-1]  # Node identifier without last character
            direction = node[-1]  # Last character, expected to be "+" or "-"

            if node_id not in my_hist:
                # Initialize with [forward, reverse] visits and a count of 1
                my_hist[node_id] = [[1, 0], 1] if direction == "+" else [[0, 1], 1]
            else:
                # Update forward/reverse visit count and total occurrence count
                if direction == "+":
                    my_hist[node_id][0][0] += 1
                else:
                    my_hist[node_id][0][1] += 1
                my_hist[node_id][1] += 1  # Increment total occurrence count

    # Find intersection of nodes across all paths
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
    plt.gcf().text(0.5, -0.15, description, wrap=True, ha='center', fontsize=10, alpha=0.8)

    # Save the plot
    output_path = f"{graph_name}_Node_Distribution.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
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
    save_all_nodes_info(largest_nodes_indices, nodes, node_sizes, node_count_hist, graph_name)

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


def extract_graph_stats(node_count_hist, edges_count, nodes_num, edges_num, avg_degree, num_paths, name, vcf_file, c_dist):
    """
    Extracts and saves graph statistics, details of nodes and edges, and VCF polymorphism stats.
    Graph stats, VCF stats, and substitution types are saved in one sheet, while nodes and edges are saved in separate Excel files.
    """
    print("Extracting and saving graph statistics...")

    # Prepare graph statistics
    graph_stats = {
        "Statistic": [
            "Number of nodes",
            "Number of unique nodes",
            "Number of shared nodes",
            "Number of core nodes",
            "Number of edges",
            "Average degree",
            "Graph density",
            "Number of paths",
            "Average nucleotide diversity"
        ],
        "Value": [
            nodes_num,
            c_dist[1],
            sum(c_dist[i] for i in c_dist.keys() if  2 <= i <= num_paths - 1),
            sum(c_dist[i] for i in c_dist.keys() if  num_paths <= i <= max(c_dist.values())),
            edges_num,
            f"{avg_degree:.2f}",
            f"{calculate_graph_density(nodes_num, edges_num):.12f}",
            num_paths,
            calculate_genomewide_pi('.'.join(vcf_file.split('.')[:-1]))
        ],
        "Explanation": [
            "Total number of nodes (sequences segments) in the graph.",
            "Total number of unique (private) nodes in the graph ,nodes with count == 1.",
            "Total number of shared nodes in the graph ,nodes with count in range [2, #paths - 1].",
            "Total number of core nodes in the graph ,nodes with count in range [#paths, MaxNodeCount].",
            "Total number of edges (connections) between the nodes.",
            "Average number of connections per node.",
            "Graph density: observed vs. maximum possible edges in the graph D = |E| / (|V|, 2). the result is multiplied by 1000 for clear interpret",
            "Total number of unique paths traversing the graph.",
            "The Average Nucleotide Diversity (π) based on biallelic SNPs from the obtained VCF file."
        ]
    }

    # Extract nodes information
    nodes_data = pd.DataFrame({
        'Node ID': list(node_count_hist.keys()),
        'Forward Count': [node_count_hist[node][0][0] for node in node_count_hist],
        'Reverse Count': [node_count_hist[node][0][1] for node in node_count_hist],
        'Total Count': [node_count_hist[node][1] for node in node_count_hist]
    })

    # Extract edges information
    edges_data = pd.DataFrame({
        'From Node': [edge[0] for edge in edges_count.keys()],
        'To Node': [edge[1] for edge in edges_count.keys()],
        'Weight': list(edges_count.values())
    })

    # Parse VCF stats
    vcf_stats = {
        "Samples": None,
        "Total Records": None,
        "No ALTs": None,
        "SNP Count": None,
        "MNP Count": None,
        "Indel Count": None,
        "Other Count": None,
        "Multiallelic Sites": None,
        "Multiallelic SNP Sites": None,
        "Transitions": None,
        "Transversions": None,
        "Transition/Transversion Ratio": None,
        "Substitution Types": {}
    }

    with open(vcf_file, "r") as f:
        for line in f:
            if line.startswith("SN"):
                parts = line.strip().split("\t")
                key = parts[2]
                value = int(parts[3])
                if key == "number of samples:":
                    vcf_stats["Samples"] = value
                elif key == "number of records:":
                    vcf_stats["Total Records"] = value
                elif key == "number of no-ALTs:":
                    vcf_stats["No ALTs"] = value
                elif key == "number of SNPs:":
                    vcf_stats["SNP Count"] = value
                elif key == "number of MNPs:":
                    vcf_stats["MNP Count"] = value
                elif key == "number of indels:":
                    vcf_stats["Indel Count"] = value
                elif key == "number of others:":
                    vcf_stats["Other Count"] = value
                elif key == "number of multiallelic sites:":
                    vcf_stats["Multiallelic Sites"] = value
                elif key == "number of multiallelic SNP sites:":
                    vcf_stats["Multiallelic SNP Sites"] = value
            elif line.startswith("TSTV"):
                parts = line.strip().split("\t")
                vcf_stats["Transitions"] = int(parts[2])
                vcf_stats["Transversions"] = int(parts[3])
                vcf_stats["Transition/Transversion Ratio"] = float(parts[4])
            elif line.startswith("ST"):
                parts = line.strip().split("\t")
                substitution_type = parts[2]
                count = int(parts[3])
                vcf_stats["Substitution Types"][substitution_type] = count

    # Combine all stats into a single sheet
    stats_combined = pd.DataFrame(graph_stats)
    vcf_data = pd.DataFrame({
        "Statistic": [
            "Samples", "Total Records", "No ALTs", "SNP Count", "MNP Count",
            "Indel Count", "Other Count", "Multiallelic Sites", "Multiallelic SNP Sites",
            "Transitions", "Transversions", "Transition/Transversion Ratio"
        ],
        "Value": [
            vcf_stats["Samples"], vcf_stats["Total Records"], vcf_stats["No ALTs"],
            vcf_stats["SNP Count"], vcf_stats["MNP Count"], vcf_stats["Indel Count"],
            vcf_stats["Other Count"], vcf_stats["Multiallelic Sites"], vcf_stats["Multiallelic SNP Sites"],
            vcf_stats["Transitions"], vcf_stats["Transversions"], vcf_stats["Transition/Transversion Ratio"]
        ]
    })

    substitution_types_df = pd.DataFrame(vcf_stats["Substitution Types"].items(),
                                         columns=["Substitution Type", "Count"])

    with pd.ExcelWriter(f"{name}.xlsx") as writer:
        # Write graph stats at the top
        stats_combined.to_excel(writer, sheet_name="All Stats", startrow=0, index=False)

        # Write VCF data below
        vcf_data.to_excel(writer, sheet_name="All Stats", startrow=len(stats_combined) + 2, index=False)

        # Write substitution types further below
        substitution_types_df.to_excel(writer, sheet_name="All Stats", startrow=len(stats_combined) + len(vcf_data) + 4,
                                       index=False)

    # Save nodes data to a separate file
    nodes_data.to_excel(f"{name}_nodes.xlsx", index=False)

    # Save edges data to a separate file
    edges_data.to_excel(f"{name}_edges.xlsx", index=False)

    print(f"Graph stats saved to {name}_stats.xlsx")
    print(f"Nodes data saved to {name}_nodes.xlsx")
    print(f"Edges data saved to {name}_edges.xlsx")

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




def extract_consensus_path(nodes, edges, node_count_hist, edges_count, min_node_count=2):
    """
    Extracts the consensus path and sequence from the graph.

    Args:
        nodes (dict): Dictionary of node IDs and their sequences.
        edges (list): List of edges in the graph.
        node_count_hist (dict): Dictionary with node visit counts, structured as [[forward_count, reverse_count], total_count].
        edges_count (dict): Dictionary of edge coverage counts.
        min_node_count (int): Minimum count threshold for including a node in the consensus.

    Returns:
        tuple: Consensus path (list of node IDs) and consensus sequence (str).
    """
    # # Preprocess edges for efficient access
    # outgoing_edges_map = preprocess_edges(edges)

    # Step 1: Initialize variables
    consensus_path = []
    consensus_sequence = ""

    # Traverse nodes in topological order
    for current_node_id in nodes.keys():
        # Extract forward and reverse counts from node_count_hist
        count_data = node_count_hist.get(current_node_id, [[0, 0], 0])
        if isinstance(count_data, list) and len(count_data) == 2:
            forward_count, reverse_count = count_data[0]
        else:
            raise ValueError(f"Invalid node count format for node {current_node_id}: {count_data}")

        # Skip nodes that do not meet the minimum count threshold
        if forward_count + reverse_count < min_node_count:
            continue

        # Add the current node to the consensus path
        consensus_path.append(current_node_id)

        # Determine the most visited strand for the current node
        if forward_count >= reverse_count:
            consensus_sequence += nodes[current_node_id].strip()  # Forward strand
            current_strand = "+"
        else:
            consensus_sequence += reverse_complement(nodes[current_node_id].strip())  # Reverse strand
            current_strand = "-"

        # # Get the outgoing edges for the current node
        # if current_node_id not in outgoing_edges_map:
        #     continue  # Skip if no outgoing edges
        # outgoing_edges = outgoing_edges_map[current_node_id]
        #
        # # Select the best outgoing edge based on edge coverage
        # best_edge = max(outgoing_edges, key=lambda edge: edges_count.get((f"{current_node_id}{current_strand}", edge), 0), default=None)
        #
        # # Add the target node of the best edge to the consensus path if it exists
        # if best_edge:
        #     next_node = best_edge[:-1]  # Remove the strand identifier (+/-) to get the node ID
        #     if next_node not in consensus_path:  # Avoid loops
        #         consensus_path.append(next_node)

    return consensus_path, consensus_sequence


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
    sequence_lengths = [len(nodes[node]) for node in sorted_nodes]
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
        seq_len = len(nodes[node])
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


def plot_node_histogram_by_paths(node_histogram, num_paths, output_file="node_histogram_by_paths.png"):
    """
    Plots a histogram of nodes categorized as private, shared, and core based on the number of paths (haplotypes),
    with a log scale on the y-axis and updated categorization logic.

    :param node_histogram: dict, A dictionary with node counts as keys and frequencies as values.
    :param num_paths: int, Total number of paths (haplotypes).
    :param output_file: str, Path to save the output plot.
    """
    # Categorize the data
    core_counts = []
    shared_counts = []
    private_counts = []

    for count, freq in node_histogram.items():
        if count == 1:  # Private nodes
            private_counts.append((count, freq))
        elif 1 < count < num_paths:  # Shared nodes
            shared_counts.append((count, freq))
        elif count >= num_paths:  # Core nodes
            core_counts.append((count, freq))

    # Split categorized data for plotting
    core_x, core_y = zip(*core_counts) if core_counts else ([], [])
    shared_x, shared_y = zip(*shared_counts) if shared_counts else ([], [])
    private_x, private_y = zip(*private_counts) if private_counts else ([], [])

    # Plot
    plt.figure(figsize=(12, 12))

    # Plot each category
    plt.bar(private_x, private_y, color="blue", label="Private nodes (count = 1)")
    plt.bar(shared_x, shared_y, color="green", label="Shared nodes (1 < count < num_paths - 1)")
    plt.bar(core_x, core_y, color="red", label="Core nodes (count >= num_paths)")

    # Add labels and legend
    plt.xlabel("Node count")
    plt.ylabel("Frequency (log scale)")
    plt.title("Histogram of Node Categories")
    plt.yscale("log")  # Set the y-axis to log scale
    plt.legend()

    # Save and show
    plt.tight_layout()
    plt.savefig(output_file, dpi=600)


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
        with open(msa_out_dir / 'odgi_matrix.tsv', 'w') as fout:
            fout.write(odgi_res.stdout)
    except Exception:
        print("error writing odgi_matrix.tsv")
        raise Exception