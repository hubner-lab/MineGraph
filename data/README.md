# Example Input Data

Minimal single-file examples for every input type that MineGraph accepts. These are **demonstrations of format**, not a full runnable dataset — full inputs for the 192-genome and 709-genome Poales graphs are archived with the manuscript (see the Data availability statement).

## Files

| File | Used by | Format |
|------|---------|--------|
| `metadata.csv` | `construct --metadata` | CSV, one FASTA filename per row, first column |
| `genes.bed` | `inject --bed` | BED6 — `path_name \t start \t end \t gene_id \t score \t strand` |
| `representatives.txt` | `inject --genomes` | Plain text, one representative `path_name` per line |
| `subgraph.gfa.gz` | `tubemap -i`, `extract -- subgraph -i` | GFAv1 (gzipped to keep the repo small; decompress before use) |

## Decompressing the example subgraph

```bash
gunzip -k data/subgraph.gfa.gz   # keeps the .gz, produces data/subgraph.gfa
```

## Pipeline demos against these examples

```bash
# Interactive SequenceTubeMap view of the example subgraph
python MineGraph.py tubemap -i data/subgraph.gfa
# → open http://localhost:3210

# Extract a smaller sub-subgraph anchored on a chosen path
python MineGraph.py extract -- subgraph \
  -i data/subgraph.gfa \
  -w <(printf "Triticum_aestivum#1#1\nElymus_hystrix#1#1\n") \
  -o triticum_elymus_mini.gfa

# Inject gene arrows (requires an .og — convert first if needed)
python MineGraph.py convert vg -i data/subgraph.gfa -o subgraph.vg   # or build .og from .gfa
python MineGraph.py inject \
  --graph subgraph.og \
  --bed data/genes.bed \
  --genomes data/representatives.txt \
  --output inject_demo --threads 4
```

## Provenance

- `metadata.csv` — first 10 entries of the 192-taxon Poales metadata list (one representative per tribe).
- `genes.bed` — inverted-repeat (IR) gene coordinates projected onto four family anchor paths (*Bromeliaceae*, *Eriocaulaceae*, *Cyperaceae*, *Typhaceae*).
- `representatives.txt` — one representative path per Poales family used for gene-arrow rendering.
- `subgraph.gfa.gz` — `Triticum × Elymus` 10-vs-10 subgraph extracted from the 192-genome chloroplast graph; used as the backing graph for Figure 4 of the manuscript.
