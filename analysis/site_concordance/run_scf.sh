#!/bin/bash
# =============================================================================
# A4+A5: PAV Node Concordance Factor (structural sCF)
# Uses the binary PAV FASTA (graph nodes as sites, A=present C=absent)
# with IQ-TREE --scf to compute concordance factors for two reference trees:
#   1. Graph-PAV MineGraph bootstrapped tree (primary)
#   2. Graph-PAV Ward tree (original, no bootstrap)
# =============================================================================

IQTREE="/Users/rakanhaib/miniconda3/bin/iqtree"
BASE="/Users/rakanhaib/GenomeCenter/revise/claude_poales_revision_ready"
TREES="$BASE/phylogenetic_trees"
FASTA="$TREES/pav_binary_variable.fasta"
THREADS=4
SCF_REPS=1000

echo "=== PAV Node Concordance Factor Analysis ==="
echo "Alignment: $FASTA"
echo "Sites: 285728 binary (A=present, C=absent)"
echo "Model: MK (Mk Lewis 2001, 2-state binary)"
echo "sCF replicates: $SCF_REPS"
echo ""

# ---- Run 1: MineGraph bootstrapped PAV tree as reference ----
REF1="$TREES/192_graph_minegraph_bs_tree_named.nwk"
OUT1="$TREES/pav_scf_minegraph_bs"
echo "--- Run 1: Reference = MineGraph_BS tree ---"
$IQTREE \
  -s  "$FASTA" \
  -m  MK+G4 \
  -te "$REF1" \
  --scf $SCF_REPS \
  --seed 42 \
  --prefix "$OUT1" \
  -T $THREADS \
  --redo \
  2>&1 | tee "${OUT1}.log"

echo ""
echo "--- Run 2: Reference = PAV Ward tree ---"
REF2="$TREES/192_graph_pav_tree_named.nwk"
OUT2="$TREES/pav_scf_pav_ward"
$IQTREE \
  -s  "$FASTA" \
  -m  MK+G4 \
  -te "$REF2" \
  --scf $SCF_REPS \
  --seed 42 \
  --prefix "$OUT2" \
  -T $THREADS \
  --redo \
  2>&1 | tee "${OUT2}.log"

echo ""
echo "=== sCF runs complete ==="
echo "Output files:"
echo "  MineGraph BS: ${OUT1}.cf.branch"
echo "  PAV Ward:     ${OUT2}.cf.branch"
