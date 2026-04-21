#!/usr/bin/env python3
"""
Process the 709-taxon MineGraph bootstrapped tree:
  - Strip _1_0 panSN suffix from tip names
  - Write named tree to phylogenetic_trees/709_graph_minegraph_bs_tree_named.nwk
  - Compute and print bootstrap support statistics
  - Write bootstrap support TSV to phylogenetic_trees/709_bs_support_summary.csv
"""

import re, sys
import numpy as np

BASE  = "/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
TREE_IN  = f"{BASE}/graphs/709_graph/MineGraph_output/phylogenetics_msa/graph_phylo_tree.tree"
TREE_OUT = f"{BASE}/phylogenetic_trees/709_graph_minegraph_bs_tree_named.nwk"
CSV_OUT  = f"{BASE}/phylogenetic_trees/709_bs_support_summary.csv"

print("Reading 709-taxon MineGraph tree ...", flush=True)
with open(TREE_IN) as f:
    newick = f.read().strip()

# --------------------------------------------------------------------------
# 1. Strip _1_0 suffix from tip names
#    Pattern: word characters ending in _1_0 immediately before a colon
# --------------------------------------------------------------------------
# Replace  SomeName_1_0:branch  →  SomeName:branch
newick_named = re.sub(r'(_1_0)(:[0-9\-\.])', r'\2', newick)

# Verify tip count unchanged
original_tips  = re.findall(r'[A-Za-z][A-Za-z0-9_]+_1_0(?=:[0-9\-])', newick)
remaining_1_0  = re.findall(r'[A-Za-z][A-Za-z0-9_]+_1_0(?=:[0-9\-])', newick_named)
clean_tips     = re.findall(r'([A-Za-z][A-Za-z0-9_]+)(?=:[0-9\-])', newick_named)
# Unique tip names (exclude bootstrap integers like 100)
tip_names = [t for t in clean_tips if not t.isdigit() and len(t) > 3 and t[0].isupper()]

print(f"  Original _1_0 tip labels:  {len(original_tips)}")
print(f"  Remaining _1_0 after sub:  {len(remaining_1_0)}")
print(f"  Clean tip labels found:    {len(set(tip_names))}")

with open(TREE_OUT, 'w') as f:
    f.write(newick_named + '\n')
print(f"Saved named tree: {TREE_OUT}", flush=True)

# --------------------------------------------------------------------------
# 2. Extract bootstrap values from internal nodes
#    In IQ-TREE Newick: bootstrap appears as node label: )100:0.xxx or )95:0.xxx
# --------------------------------------------------------------------------
bs_values = [int(m) for m in re.findall(r'\)(\d+):', newick_named) if m.strip().isdigit()]
bs_arr = np.array(bs_values, dtype=float)

print(f"\n--- Bootstrap Support Statistics (709-taxon MineGraph tree) ---")
print(f"  Internal nodes with BS:   {len(bs_arr)}")
print(f"  Mean BS:                  {bs_arr.mean():.1f}%")
print(f"  Median BS:                {np.median(bs_arr):.0f}%")
print(f"  SD BS:                    {bs_arr.std():.1f}%")
print(f"  % >= 50:                  {100*np.mean(bs_arr >= 50):.1f}%")
print(f"  % >= 70:                  {100*np.mean(bs_arr >= 70):.1f}%")
print(f"  % >= 90:                  {100*np.mean(bs_arr >= 90):.1f}%")
print(f"  % >= 95:                  {100*np.mean(bs_arr >= 95):.1f}%")
print(f"  % = 100:                  {100*np.mean(bs_arr == 100):.1f}%")

# Write CSV
with open(CSV_OUT, 'w') as f:
    f.write("tree,n_internal,mean_BS,median_BS,sd_BS,pct_50,pct_70,pct_90,pct_95,pct_100\n")
    f.write(
        f'"Graph_MineGraph_BS_709 (1000 replicates, Felsenstein)",'
        f'{len(bs_arr)},{bs_arr.mean():.1f},{np.median(bs_arr):.0f},{bs_arr.std():.1f},'
        f'{100*np.mean(bs_arr >= 50):.1f},{100*np.mean(bs_arr >= 70):.1f},'
        f'{100*np.mean(bs_arr >= 90):.1f},{100*np.mean(bs_arr >= 95):.1f},'
        f'{100*np.mean(bs_arr == 100):.1f}\n'
    )
print(f"\nSaved BS statistics: {CSV_OUT}", flush=True)
print("=== process_709_tree.py DONE ===")
