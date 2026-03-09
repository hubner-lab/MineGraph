#!/usr/bin/env python3
"""
MineGraph.py – Host-side CLI wrapper.

Orchestrates all pipeline steps by launching Docker containers on the host.
No step runs docker from inside another container.

Pipeline (construct command):
  Step 1 – prepare_and_mash_input.py  → rakanhaib/opggb
  Step 2 – RepeatMasker               → pegi3s/repeat_masker
  Step 3 – run_repeatmask.py          → rakanhaib/opggb
  Step 4 – run_pggb.py                → rakanhaib/opggb
  Step 5 – run_stats.py               → rakanhaib/opggb
"""

import argparse
import os
import platform
import shlex
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

# =============================================================================
# Constants
# =============================================================================

OPGGB_IMAGE_DEFAULT  = "rakanhaib/opggb:latest"
TUBEMAP_IMAGE_DEFAULT = "rakanhaib/sequencetubemap:latest"
REPEATMASKER_IMAGE   = "pegi3s/repeat_masker"

TUBEMAP_HOST_PORT      = 3210
TUBEMAP_CONTAINER_PORT = 3000

# Scripts inside the opggb container
_PREPARE_SCRIPT      = "/src/prepare_and_mash_input.py"
_EXTRACT_TR_SCRIPT   = "/src/run_repeatmask.py"
_PGGB_SCRIPT         = "/src/run_pggb.py"
_STATS_SCRIPT        = "/src/run_stats.py"
_EXTRACT_SCRIPT      = "/src/extract.py"
_GFA2FASTA_SCRIPT    = "/src/gfa2fasta.py"

# Path flags used by collect_paths_and_rewrite_args for extract/convert
INPUT_PATH_FLAGS = {
    "extract": {"-i", "--input", "-w", "--wanted"},
    "convert": {"-i", "--input"},
}
OUTPUT_PATH_FLAGS = {
    "extract": {"-o", "--output"},
    "convert": {"-o", "--output"},
}


# =============================================================================
# Helpers
# =============================================================================

def _platform_flags() -> List[str]:
    m = platform.machine()
    if m == "x86_64":
        return ["--platform", "linux/amd64"]
    if m in ("arm64", "aarch64"):
        return ["--platform", "linux/arm64"]
    return []


def _looks_like_flag(s: str) -> bool:
    return s.startswith("-") and s != "-"


def _run(cmd: List[str], dry_run: bool, label: str = "") -> int:
    if label:
        print(f"[INFO] {label}")
    print("[CMD] " + " ".join(shlex.quote(x) for x in cmd))
    if dry_run:
        print("[DRY-RUN] skipped.")
        return 0
    return subprocess.run(cmd).returncode


# =============================================================================
# Path rewriting for extract / convert (unchanged from original)
# =============================================================================

def collect_paths_and_rewrite_args(
    command: str,
    forwarded_args: List[str],
    mkdir: bool,
) -> Tuple[List[str], Dict[Path, str], List[str]]:
    in_flags  = INPUT_PATH_FLAGS.get(command, set())
    out_flags = OUTPUT_PATH_FLAGS.get(command, set())

    mount_map: Dict[Path, str] = {}
    messages:  List[str] = []
    new_args = forwarded_args[:]
    vol_idx  = 0

    def hostdir_to_container(host_dir: Path) -> str:
        nonlocal vol_idx
        host_dir = host_dir.resolve(strict=False)
        if host_dir not in mount_map:
            mount_map[host_dir] = f"/data/vol{vol_idx}"
            vol_idx += 1
        return mount_map[host_dir]

    def handle_path(flag: str, val: str, is_output: bool) -> str:
        host_path = Path(val).expanduser()
        if is_output:
            parent = host_path.parent.expanduser()
            if not parent.exists():
                if mkdir:
                    parent.mkdir(parents=True, exist_ok=True)
                else:
                    messages.append(
                        f"[ERROR] Output directory does not exist for {flag}: {parent}\n"
                        f"        Re-run with --mkdir to auto-create."
                    )
                    return val
            container_dir = hostdir_to_container(parent.resolve(strict=False))
            return f"{container_dir}/{host_path.name}"
        if not host_path.exists():
            messages.append(f"[ERROR] Input path not found for {flag}: {host_path}")
            return val
        container_dir = hostdir_to_container(host_path.resolve(strict=False).parent)
        return f"{container_dir}/{host_path.name}"

    i = 0
    while i < len(new_args):
        tok = new_args[i]
        if tok in in_flags | out_flags and i + 1 < len(new_args):
            val = new_args[i + 1]
            if not _looks_like_flag(val):
                new_args[i + 1] = handle_path(tok, val, tok in out_flags)
                i += 2
                continue
        if tok.startswith("--") and "=" in tok:
            flag, val = tok.split("=", 1)
            if flag in in_flags | out_flags:
                new_args[i] = f"{flag}={handle_path(flag, val, flag in out_flags)}"
        i += 1

    return new_args, mount_map, messages


# =============================================================================
# Single-step opggb runner
# =============================================================================

def _opggb(
    image: str,
    script: str,
    script_args: List[str],
    data_dir_host:   str,
    output_dir_host: str,
    extra_mounts: List[str],
    dry_run: bool,
    label: str = "",
) -> int:
    """
    Run one script inside the opggb container.

    Mounts:
      <data_dir_host>   → /data      (read-only input FASTA directory)
      <output_dir_host> → /output    (read-write; matches hardcoded paths in scripts)
    """
    # output_dir_host may not exist yet — create first, then resolve strictly
    Path(output_dir_host).expanduser().mkdir(parents=True, exist_ok=True)
    abs_data = str(Path(data_dir_host).expanduser().resolve(strict=False))
    abs_out  = str(Path(output_dir_host).expanduser().resolve())

    cmd = (
        ["docker", "run", "--rm"]
        + _platform_flags()
        + ["-v", f"{abs_data}:/data",
           "-v", f"{abs_out}:/output"]
    )
    for m in extra_mounts:
        cmd += ["-v", m]
    cmd += [image, "python", script] + script_args

    return _run(cmd, dry_run, label)


# =============================================================================
# construct — five-step host-orchestrated pipeline
# =============================================================================

def run_construct(args) -> int:
    """
    Orchestrate the full construct pipeline from the host.
    Each step is a separate docker run — no docker-in-docker.
    """
    # Strip leading '--' separator if present
    fwd = args.forward_args
    if fwd and fwd[0] == "--":
        fwd = fwd[1:]

    # Parse the forwarded construct arguments so we can extract paths and flags
    cp = argparse.ArgumentParser(add_help=False)
    cp.add_argument("--data_dir",    required=True)
    cp.add_argument("--output_dir",  required=True)
    cp.add_argument("--metadata",    required=True)
    cp.add_argument("--threads",     default="16")
    cp.add_argument("--tree_pars",   default="10")
    cp.add_argument("--tree_bs",     default="10")
    cp.add_argument("--quantile",    default="25")
    cp.add_argument("--top_n",       default="50")
    cp.add_argument("--window_size", default="1000")
    cp.add_argument("--mode",        default="all")
    cp.add_argument("--tree_type",   default="graph")
    cp.add_argument("--sc_th",       default="5")
    cp.add_argument("--phyl-tree",   action="store_true", dest="phyl_tree")
    cp.add_argument("--plots",       action="store_true")
    cp.add_argument("--only-stats",  action="store_true", dest="only_stats")
    cp.add_argument("--sq_view",     action="store_true")

    try:
        ca, _ = cp.parse_known_args(fwd)
    except SystemExit:
        print("[ERROR] Could not parse construct arguments. Run with --help.", file=sys.stderr)
        return 2

    data_dir    = ca.data_dir
    output_dir  = ca.output_dir
    metadata    = ca.metadata
    threads     = ca.threads
    mode        = ca.mode

    # Resolve host paths.
    # strict=False: without it Python 3.12 raises FileNotFoundError for any
    # path that doesn't exist yet (e.g. the output_dir we are about to create).
    abs_data   = str(Path(data_dir).expanduser().resolve(strict=False))
    abs_meta   = str(Path(metadata).expanduser().resolve(strict=False))
    meta_dir   = str(Path(abs_meta).parent)
    meta_name  = Path(abs_meta).name

    # Validate inputs exist BEFORE touching the output dir
    if not Path(abs_data).exists():
        print(f"[ERROR] data_dir not found: {abs_data}", file=sys.stderr)
        return 2
    if not Path(abs_meta).exists():
        print(f"[ERROR] metadata not found: {abs_meta}", file=sys.stderr)
        return 2

    # Create output dir first, then resolve with strict=True so the
    # resolved path is canonical and real (no symlink ambiguity).
    Path(output_dir).expanduser().mkdir(parents=True, exist_ok=True)
    abs_out = str(Path(output_dir).expanduser().resolve())

    dry_run = args.dry_run
    image   = args.image

    # ------------------------------------------------------------------
    # STEP 1 – Prepare FASTA input (opggb container)
    # ------------------------------------------------------------------
    if mode in ("all", "extract-tr", "construct-graph"):
        rc = _opggb(
            image=image,
            script=_PREPARE_SCRIPT,
            script_args=["/data", "/output", "--output_yaml", "/output/params.yaml",
                         "--metadata", f"/meta/{meta_name}"],
            data_dir_host=abs_data,
            output_dir_host=abs_out,
            extra_mounts=[f"{meta_dir}:/meta:ro"],
            dry_run=dry_run,
            label="Step 1 — Prepare FASTA input",
        )
        if rc != 0:
            return rc

        # ------------------------------------------------------------------
        # STEP 2 – RepeatMasker (dedicated pegi3s image — only docker call
        #          that is not opggb, and it runs on the HOST so no nesting)
        # ------------------------------------------------------------------
        print("[INFO] Step 2 — RepeatMasker")
        rm_cmd = (
            ["docker", "run", "--rm"]
            + _platform_flags()
            + [
                "-v", f"{abs_out}:/data",
                REPEATMASKER_IMAGE,
                "bash", "-c",
                f"RepeatMasker /data/downsampled_panSN_output.fasta "
                f"-no_is -pa {threads} -s -gff",
            ]
        )
        rc = _run(rm_cmd, dry_run, "")
        if rc != 0:
            return rc

        # ------------------------------------------------------------------
        # STEP 3 – Parse RepeatMasker output → update params.yaml (opggb)
        # ------------------------------------------------------------------
        rc = _opggb(
            image=image,
            script=_EXTRACT_TR_SCRIPT,
            script_args=[],
            data_dir_host=abs_out,   # /data = output_dir (has the .out file)
            output_dir_host=abs_out,
            extra_mounts=[],
            dry_run=dry_run,
            label="Step 3 — Extract longest TR / update params.yaml",
        )
        if rc != 0:
            return rc

    if mode == "extract-tr":
        print("[INFO] Extract-only mode complete.")
        return 0

    # ------------------------------------------------------------------
    # STEP 4 – PGGB graph construction (opggb)
    # ------------------------------------------------------------------
    rc = _opggb(
        image=image,
        script=_PGGB_SCRIPT,
        script_args=[threads],
        data_dir_host=abs_out,   # /data = output_dir (has params.yaml)
        output_dir_host=abs_out,
        extra_mounts=[],
        dry_run=dry_run,
        label="Step 4 — PGGB graph construction",
    )
    if rc != 0:
        return rc

    if mode == "construct-graph":
        print("[INFO] Graph construction complete. Skipping statistics.")
        return 0

    # ------------------------------------------------------------------
    # STEP 5 – Graph statistics (opggb)
    # ------------------------------------------------------------------
    stats_args = [
        "--input_dir",  "/output/pggb_output",
        "--output_dir", "/output",
        "--threads",    threads,
        "--tree_pars",  ca.tree_pars,
        "--tree_bs",    ca.tree_bs,
        "--quantile",   str(float(ca.quantile) / 100.0),
        "--top_n",      ca.top_n,
        "--window_size", ca.window_size,
        "--tree_type",  ca.tree_type,
        "--softcore_threshold", ca.sc_th,
    ]
    if ca.phyl_tree:   stats_args.append("--phyl-tree")
    if ca.plots:       stats_args.append("--plots")
    if ca.only_stats:  stats_args.append("--only-stats")

    rc = _opggb(
        image=image,
        script=_STATS_SCRIPT,
        script_args=stats_args,
        data_dir_host=abs_out,
        output_dir_host=abs_out,
        extra_mounts=[],
        dry_run=dry_run,
        label="Step 5 — Graph statistics",
    )
    if rc != 0:
        return rc

    # ------------------------------------------------------------------
    # Optional: SequenceTubeMap viewer
    # ------------------------------------------------------------------
    if ca.sq_view and not dry_run:
        return run_tubemap(
            image=args.tubemap_image,
            input_graph_host=str(Path(abs_out) / "MineGraph_output" / "gfa_convert" / "gfa_to_vg.vg"),
            interactive=False,
            dry_run=dry_run,
        )

    return 0


# =============================================================================
# extract / convert  (unchanged behaviour, kept simple)
# =============================================================================

def _opggb_simple(
    image: str,
    script: str,
    forwarded_args: List[str],
    mounts: List[str],
    dry_run: bool,
) -> int:
    cwd = Path.cwd().resolve(strict=False)
    cmd = (
        ["docker", "run", "--rm"]
        + _platform_flags()
        + ["-v", f"{cwd}:/work", "-w", "/work"]
    )
    for m in mounts:
        cmd += ["-v", m]
    cmd += [image, "python", script] + forwarded_args
    return _run(cmd, dry_run, "Running:")


def run_extract(args) -> int:
    fwd = args.forward_args
    if fwd and fwd[0] == "--":
        fwd = fwd[1:]
    rewritten, mount_map, msgs = collect_paths_and_rewrite_args("extract", fwd, args.mkdir)
    for m in msgs:
        print(m, file=sys.stderr)
    if any(m.startswith("[ERROR]") for m in msgs):
        return 2
    auto_mounts = [f"{h}:{c}" for h, c in mount_map.items()]
    return _opggb_simple(args.image, _EXTRACT_SCRIPT, rewritten,
                         args.mount + auto_mounts, args.dry_run)


def run_convert(args) -> int:
    in_path  = Path(args.input).expanduser().resolve(strict=False)
    out_path = Path(args.output).expanduser().resolve(strict=False)

    if not in_path.exists():
        print(f"[ERROR] Input not found: {in_path}", file=sys.stderr)
        return 2
    if not out_path.parent.exists():
        if args.mkdir:
            out_path.parent.mkdir(parents=True, exist_ok=True)
        else:
            print(f"[ERROR] Output dir missing: {out_path.parent}  (use --mkdir)", file=sys.stderr)
            return 2

    in_dir  = in_path.parent
    out_dir = out_path.parent
    mounts  = [f"{in_dir}:/data/in"] + args.mount
    if out_dir != in_dir:
        mounts += [f"{out_dir}:/data/out"]
        out_c = f"/data/out/{out_path.name}"
    else:
        out_c = f"/data/in/{out_path.name}"
    in_c = f"/data/in/{in_path.name}"

    cwd = Path.cwd().resolve(strict=False)
    cmd = (
        ["docker", "run", "--rm"]
        + _platform_flags()
        + ["-v", f"{cwd}:/work", "-w", "/work"]
    )
    for m in mounts:
        cmd += ["-v", m]

    if args.mode == "fasta":
        cmd += [args.image, "python", _GFA2FASTA_SCRIPT, in_c, out_c]
    else:
        cmd += [args.image, "bash", "-lc",
                f"set -euo pipefail; vg convert -g {shlex.quote(in_c)} > {shlex.quote(out_c)}"]

    return _run(cmd, args.dry_run, "Converting:")


# =============================================================================
# tubemap
# =============================================================================

def run_tubemap(
    image: str,
    input_graph_host: str,
    interactive: bool,
    dry_run: bool,
) -> int:
    graph = Path(input_graph_host).expanduser().resolve(strict=False)
    if not graph.exists():
        print(f"[ERROR] Input graph not found: {graph}", file=sys.stderr)
        return 2

    cname = f"minegraph_tubemap_{os.getpid()}_{int(time.time())}"
    cmd = ["docker", "run", "--rm", "--name", cname]
    if interactive:
        cmd.append("-it")
    cmd += [
        "-v", f"{graph.parent}:/in:ro",
        "-p", f"{TUBEMAP_HOST_PORT}:{TUBEMAP_CONTAINER_PORT}",
        image,
        "bash", "-lc",
        f"""set -euo pipefail
mkdir -p /data
cp "/in/{shlex.quote(graph.name)}" /data/input_graph
./scripts/prepare_vg.sh /data/input_graph
XG=$(find /data -name '*.xg' | head -n1)
[ -z "$XG" ] && {{ echo '[ERROR] No .xg produced'; exit 1; }}
mkdir -p exampleData && cp "$XG" exampleData/
echo '[INFO] Open http://localhost:{TUBEMAP_HOST_PORT}'
npm start serve""",
    ]

    if dry_run:
        print("[DRY-RUN] " + " ".join(shlex.quote(x) for x in cmd))
        return 0

    print(f"[INFO] TubeMap container: {cname}")
    print(" ".join(shlex.quote(x) for x in cmd))
    proc = None
    try:
        proc = subprocess.Popen(cmd)
        return proc.wait()
    except KeyboardInterrupt:
        print("\n[INFO] Stopping TubeMap...", file=sys.stderr)
        return 130
    finally:
        for cleanup in (["docker", "stop", "-t", "1", cname],
                        ["docker", "rm",   "-f",     cname]):
            subprocess.run(cleanup, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if proc and proc.poll() is None:
            proc.terminate()


# =============================================================================
# CLI
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="MineGraph.py",
        description="MineGraph CLI — host-side Docker orchestrator",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=f"""
EXAMPLES
========
Construct (full pipeline):
  MineGraph.py construct -- --data_dir ./mito_Data --output_dir ./out \\
    --metadata ./samples.csv --threads 32 --phyl-tree --plots

Extract subgraph:
  MineGraph.py extract -- subgraph -i graph.gfa -w wanted.txt -o sub.gfa

Convert GFA → FASTA:
  MineGraph.py convert fasta -i graph.gfa -o graph.fasta

TubeMap (http://localhost:{TUBEMAP_HOST_PORT}):
  MineGraph.py tubemap -i out/MineGraph_output/gfa_convert/gfa_to_vg.vg
""",
    )

    parser.add_argument("--image",         default=OPGGB_IMAGE_DEFAULT,
                        help=f"opggb image (default: {OPGGB_IMAGE_DEFAULT})")
    parser.add_argument("--tubemap-image", default=TUBEMAP_IMAGE_DEFAULT,
                        help=f"TubeMap image (default: {TUBEMAP_IMAGE_DEFAULT})")
    parser.add_argument("-m", "--mount",   action="append", default=[],
                        help="Extra host:container mount (repeatable)")
    parser.add_argument("--mkdir",         action="store_true",
                        help="Auto-create missing output directories")
    parser.add_argument("--dry-run",       action="store_true",
                        help="Print docker commands without running")

    sub = parser.add_subparsers(dest="command", required=True)

    # construct
    pc = sub.add_parser("construct", help="Run full pangenome construction pipeline")
    pc.add_argument("forward_args", nargs=argparse.REMAINDER,
                    help="Arguments forwarded to construct steps (use '--' separator)")

    # extract
    pe = sub.add_parser("extract", help="Subgraph extraction / graph inspection")
    pe.add_argument("forward_args", nargs=argparse.REMAINDER,
                    help="Arguments forwarded to extract.py (use '--' separator)")

    # convert
    pconv = sub.add_parser("convert", help="Convert GFA to FASTA or VG")
    pconv.add_argument("mode", choices=["fasta", "vg"])
    pconv.add_argument("-i", "--input",  required=True, help="Input GFA file")
    pconv.add_argument("-o", "--output", required=True, help="Output file")

    # tubemap
    pt = sub.add_parser("tubemap",
                        help=f"Interactive graph browser (http://localhost:{TUBEMAP_HOST_PORT})")
    pt.add_argument("-i", "--input",       required=True, help="Input .gfa or .vg")
    pt.add_argument("--interactive",       action="store_true", help="Run with -it (debug)")

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args   = parser.parse_args(argv)

    if args.command == "construct":
        return run_construct(args)
    if args.command == "extract":
        return run_extract(args)
    if args.command == "convert":
        return run_convert(args)
    if args.command == "tubemap":
        return run_tubemap(
            image=args.tubemap_image,
            input_graph_host=args.input,
            interactive=args.interactive,
            dry_run=args.dry_run,
        )
    parser.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
