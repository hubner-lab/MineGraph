#!/usr/bin/env python3
"""
Convert odgi_matrix.tsv (PAV coverage matrix) to binary FASTA alignment.
Each graph node = one alignment site (column).
A = present (coverage > 0), C = absent (coverage = 0).
Only variable sites (1 to N-1 taxa having the node) are retained.
Taxon names are cleaned to match the tree labels (strip #1#0 suffix).
"""

import sys, re, numpy as np

MATRIX = "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready/graphs/192_graph/MineGraph_output/phylogenetics_msa/odgi_matrix.tsv"
OUT    = "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready/phylogenetic_trees/pav_binary_variable.fasta"

# Name corrections: odgi_matrix raw → renamed tree label
# (mirrors the corrections applied when renaming 192_graph_minegraph_bs_tree)
RENAME_MAP = {
    'Coix_lacryma_jobi':            'Coix_lacrymajobi',
    'Echinochloa_crus_galli_var__pr':'Echinochloa_crusgalli_var__pra',
    'Fargesia_sp__YYZ_2019':        'Fargesia_sp__YYZ2019',
    'Lindmania_sp__AMF_2023a':      'Lindmania_sp__AMF2023a',
    'Neomicrocalamus_sp__CYZ_2023a':'Neomicrocalamus_sp__CYZ2023a',
    'Taeniatherum_caput_medusae_sub':'Taeniatherum_caputmedusae_subs',
    'Temochloa_sp__CYZ_2023b':      'Temochloa_sp__CYZ2023b',
    'Diheteropogon_amplectens_var': 'Diheteropogon_amplectens_var__',
    'Hordeum_brevisubulatum_subsp': 'Hordeum_brevisubulatum_subsp__',
    'Saccharum_hybrid_cultivar_NCo':'Saccharum_hybrid_cultivar_NCo_',
    'Schizachyrium_sanguineum_var': 'Schizachyrium_sanguineum_var__',
    'x_Triticosecale_sp':           'x_Triticosecale_sp_',
}

print("Reading odgi_matrix.tsv ...", flush=True)
names = []
rows  = []

with open(MATRIX) as f:
    header = f.readline()  # skip header
    for line in f:
        parts = line.rstrip('\n').split('\t')
        # Clean taxon name: strip #1#0 suffix, replace spaces
        raw_name = parts[0]
        clean    = raw_name.split('#')[0]          # remove #1#0
        # Preserve all underscores as-is (tree labels also use __ for subspecies)
        clean    = re.sub(r'[^A-Za-z0-9_]', '_', clean).strip('_')
        # Apply same name corrections used when renaming the MineGraph tree
        clean = RENAME_MAP.get(clean, clean)
        names.append(clean)
        # Binary: >0 → 1 (A), 0 → 0 (C)
        rows.append([1 if int(x) > 0 else 0 for x in parts[3:]])

mat      = np.array(rows, dtype=np.uint8)   # (192, 288104)
n_taxa   = mat.shape[0]
n_nodes  = mat.shape[1]
print(f"  Matrix: {n_taxa} taxa × {n_nodes} nodes", flush=True)

# Keep only variable sites (at least 1 but not all taxa have the node)
col_sums = mat.sum(axis=0)
variable = (col_sums > 0) & (col_sums < n_taxa)
mat_var  = mat[:, variable]
print(f"  Variable sites: {mat_var.shape[1]}", flush=True)

# Write FASTA — encode 1→A, 0→C
CHUNK = 60  # line wrap
chars = np.array(['C', 'A'], dtype='U1')

print(f"Writing FASTA to {OUT} ...", flush=True)
with open(OUT, 'w') as out:
    for i, name in enumerate(names):
        seq = ''.join(chars[mat_var[i]])
        out.write(f'>{name}\n')
        # Write wrapped
        for start in range(0, len(seq), CHUNK):
            out.write(seq[start:start+CHUNK] + '\n')
        if (i+1) % 20 == 0:
            print(f"  Written {i+1}/{n_taxa} sequences ...", flush=True)

print(f"Done. FASTA saved to {OUT}", flush=True)
print(f"Alignment: {n_taxa} taxa × {mat_var.shape[1]} sites (binary, A=present C=absent)", flush=True)
