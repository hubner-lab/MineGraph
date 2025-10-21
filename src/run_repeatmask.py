import re
import yaml
from collections import defaultdict
from itertools import groupby
import math

PAREN = re.compile(r'[()]')

def to_int(x):
    return int(PAREN.sub('', x))

def parse_repeatmasker_out(out_file):
    """Parse RepeatMasker .out into list of entries with typed fields."""
    entries = []
    with open(out_file, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith(("SW", "score", "=")):
                continue
            cols = s.split()
            # Expected columns (RepeatMasker .out):
            # 0 score, 1 div, 2 del, 3 ins, 4 qname,
            # 5 q_beg, 6 q_end, 7 q_left, 8 strand,
            # 9 rname, 10 classfam, 11 r_beg, 12 r_end, 13 r_left, 14 id
            if len(cols) < 14:
                continue
            try:
                e = {
                    "score": int(cols[0]),
                    "qname": cols[4],
                    "q_beg": to_int(cols[5]),
                    "q_end": to_int(cols[6]),
                    "q_left": to_int(cols[7]),
                    "strand": cols[8],
                    "rname": cols[9],
                    "classfam": cols[10],
                    "r_beg": to_int(cols[11]),
                    "r_end": to_int(cols[12]),
                    # ID may be missing or non-numeric; ignore
                }
                entries.append(e)
            except Exception:
                # Skip malformed lines
                continue
    return entries

def is_interspersed(e):
    """Transposon-like classes have a slash, e.g. LTR/Gypsy, DNA/hAT, LINE/L1, SINE/*, RC/Helitron."""
    cf = e.get("classfam", "")
    return "/" in cf

def is_simple(e):
    """Simple_repeat or Satellite classes."""
    cf = e.get("classfam", "")
    return cf.startswith("Simple_repeat") or cf.startswith("Satellite")

def is_rna(e):
    """rRNA/tRNA/etc."""
    cf = e.get("classfam", "")
    head = cf.split("/", 1)[0]
    return head in {"rRNA", "tRNA", "scRNA", "snRNA", "srpRNA"}

def longest_span_on_query(entries, pred=lambda e: True):
    """Max contiguous q_end - q_beg + 1 among entries that satisfy pred()."""
    L = 0
    for e in entries:
        if pred(e):
            L = max(L, abs(e["q_end"] - e["q_beg"]) + 1)
    return L

def longest_merged_on_query(entries, pred=lambda e: True):
    """Merge overlapping intervals per qname and return longest merged span among entries that satisfy pred()."""
    data = [e for e in entries if pred(e)]
    data.sort(key=lambda e: (e["qname"], min(e["q_beg"], e["q_end"]), max(e["q_beg"], e["q_end"])))
    best = 0
    for qname, grp in groupby(data, key=lambda e: e["qname"]):
        merged = []
        for e in grp:
            s = min(e["q_beg"], e["q_end"])
            t = max(e["q_beg"], e["q_end"])
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, t])
            else:
                merged[-1][1] = max(merged[-1][1], t)
        for s, t in merged:
            best = max(best, t - s + 1)
    return best

def compute_segment_length(out_file, yaml_file="/data/params.yaml", include_rna=True, multiplier=1.2, round_to=10):
    entries = parse_repeatmasker_out(out_file)

    # Longest segments by category on query coords
    L_inter = longest_merged_on_query(entries, is_interspersed)  # transposons
    L_simple = longest_merged_on_query(entries, is_simple)       # simple/satellite
    L_rna    = longest_merged_on_query(entries, is_rna)          # rRNA/tRNA
    L_all    = longest_merged_on_query(entries, lambda e: True)

    # Choose base length
    if include_rna:
        base = max(L_inter, L_simple, L_rna)
    else:
        base = max(L_inter, L_simple)

    seg = int(math.ceil((base * multiplier) / round_to) * round_to)

    # Update YAML
    try:
        with open(yaml_file, "r") as yf:
            data = yaml.safe_load(yf) or {}
    except FileNotFoundError:
        data = {}

    data["segment_length"] = seg
    data["repeat_lengths"] = {
        "interspersed_merged_max": L_inter,
        "simple_satellite_merged_max": L_simple,
        "rna_merged_max": L_rna,
        "all_repeats_merged_max": L_all,
        "multiplier": multiplier,
    }

    with open(yaml_file, "w") as yf:
        yaml.dump(data, yf, default_flow_style=False)

    print(f"[info] interspersed={L_inter} | simple/satellite={L_simple} | rRNA/tRNA={L_rna} | all={L_all}")
    print(f"[info] segment_length set to {seg} (include_rna={include_rna}, ×{multiplier})")

if __name__ == "__main__":
    compute_segment_length("/data/downsampled_panSN_output.fasta.out", yaml_file="/data/params.yaml", include_rna=True)
