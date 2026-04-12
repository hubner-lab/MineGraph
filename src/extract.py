#!/usr/bin/env python3
import argparse
import os
import sys
import subprocess
from typing import List, Set, Optional


def run(cmd: List[str]) -> str:
    res = subprocess.run(cmd, check=True, text=True, capture_output=True)
    return res.stdout


def require_file(path: str, what: str) -> None:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{what} not found: {path}")


def read_nonempty_lines(path_file: str) -> List[str]:
    with open(path_file, "r") as f:
        return [ln.strip() for ln in f if ln.strip()]


def write_lines(out_path: str, lines: List[str]) -> None:
    with open(out_path, "w") as f:
        for ln in lines:
            f.write(ln + "\n")


def write_text(out_path: str, text: str) -> None:
    with open(out_path, "w") as f:
        f.write(text)


# ---------------- MODES ---------------- #

def cmd_paths(graph: str, out: str, prefix: Optional[str]) -> None:
    require_file(graph, "Graph")
    paths = run(["odgi", "paths", "-i", graph, "-L"]).splitlines()
    paths = [p for p in paths if p]
    if prefix:
        paths = [p for p in paths if p.startswith(prefix)]
    if not paths:
        raise ValueError("No paths matched.")
    write_lines(out, paths)


def cmd_stats(graph: str, out: str) -> None:
    require_file(graph, "Graph")
    stats = run(["odgi", "stats", "-i", graph, "-S"])
    paths = run(["odgi", "paths", "-i", graph, "-L"])
    write_text(out, stats + "\nPATHS:\n" + paths)


def cmd_subgraph(graph: str, wanted: str, out_og: str, check_subset: bool) -> None:
    require_file(graph, "Graph")
    require_file(wanted, "Wanted paths file")

    wanted_paths = read_nonempty_lines(wanted)
    if not wanted_paths:
        raise ValueError("Wanted paths file is empty.")

    # Check graph readability
    run(["odgi", "stats", "-i", graph, "-S"])

    if check_subset:
        graph_paths = set(run(["odgi", "paths", "-i", graph, "-L"]).splitlines())
        missing = [p for p in wanted_paths if p not in graph_paths]
        if missing:
            raise ValueError(
                f"{len(missing)} wanted paths not found in graph.\n"
                f"First missing:\n" + "\n".join(missing[:10])
            )

    anchor = wanted_paths[0]
    cmd = [
        "odgi", "extract",
        "-i", graph,
        "-r", anchor,
        "-p", wanted,
        "-R", wanted,
        "-o", out_og
    ]

    print("Running:", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)

    #print(run(["odgi", "stats", "-i", out_og, "-S"]), file=sys.stderr)


# ---------------- CLI ---------------- #

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="""
ODGI utility script with three modes:

  paths     → save path names from a graph
  stats     → save graph stats + path names
  subgraph  → extract a subgraph from a wanted path list

All commands rely purely on odgi.
""",
        epilog="""
Examples:

  Save all paths:
    odgi_tool.py paths -i graph.og -o all_paths.txt

  Save Elymus paths only:
    odgi_tool.py paths -i graph.og -o elymus.txt --prefix Elymus_

  Save stats + paths:
    odgi_tool.py stats -i graph.og -o graph_info.txt

  Extract a 10v10 subgraph:
    odgi_tool.py subgraph -i graph.og -w wanted.txt -o subgraph.og
""",
        formatter_class=argparse.RawTextHelpFormatter
    )

    sub = parser.add_subparsers(dest="mode")

    p_paths = sub.add_parser("paths", help="Save path names (odgi paths -L)")
    p_paths.add_argument("-i", "--graph", required=True)
    p_paths.add_argument("-o", "--out", required=True)
    p_paths.add_argument("--prefix", help="Optional path prefix filter")

    p_stats = sub.add_parser("stats", help="Save stats (-S) and path names")
    p_stats.add_argument("-i", "--graph", required=True)
    p_stats.add_argument("-o", "--out", required=True)

    p_sub = sub.add_parser("subgraph", help="Extract subgraph from path list")
    p_sub.add_argument("-i", "--graph", required=True)
    p_sub.add_argument("-w", "--wanted", required=True)
    p_sub.add_argument("-o", "--out_og", required=True)
    p_sub.add_argument("--no_subset_check", action="store_true")

    return parser


def main():
    parser = build_parser()

    # 👇 PRINT HELP IF NO ARGUMENTS
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.mode == "paths":
        cmd_paths(args.graph, args.out, args.prefix)
    elif args.mode == "stats":
        cmd_stats(args.graph, args.out)
    elif args.mode == "subgraph":
        cmd_subgraph(args.graph, args.wanted, args.out_og, not args.no_subset_check)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()