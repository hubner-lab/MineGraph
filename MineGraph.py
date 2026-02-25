#!/usr/bin/env python3
import argparse
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

# =========================
# Docker images
# =========================
OPGGB_IMAGE_DEFAULT = "rakanhaib/opggb:latest"
TUBEMAP_IMAGE_DEFAULT = "rakanhaib/sequencetubemap:latest"

# scripts INSIDE the opggb container
CONSTRUCT_SCRIPT = "/src/construct.py"
EXTRACT_SCRIPT   = "/src/extract.py"
GFA2FASTA_SCRIPT = "/src/gfa2fasta.py"   # <-- added

# TubeMap fixed ports
TUBEMAP_HOST_PORT = 3210
TUBEMAP_CONTAINER_PORT = 3000

# =========================
# Auto-mount path handling (opggb commands)
# =========================
INPUT_PATH_FLAGS = {
    "extract": {"-i", "--input", "-w", "--wanted"},
    "construct": {"--data_dir", "--metadata"},
    "convert": {"-i", "--input"},  # <-- added
}
OUTPUT_PATH_FLAGS = {
    "extract": {"-o", "--output"},
    "construct": {"--output_dir"},
    "convert": {"-o", "--output"},  # <-- added
}

def _looks_like_flag(s: str) -> bool:
    return s.startswith("-") and s != "-"

def collect_paths_and_rewrite_args(
    command: str,
    forwarded_args: List[str],
    mkdir: bool,
) -> Tuple[List[str], Dict[Path, str], List[str]]:
    """
    For known path flags:
      - input paths must exist (mount their parent)
      - output paths mount their parent; optionally create parent (--mkdir)
    Rewrites args to /data/volN/<basename>
    """
    in_flags = INPUT_PATH_FLAGS.get(command, set())
    out_flags = OUTPUT_PATH_FLAGS.get(command, set())

    mount_map: Dict[Path, str] = {}
    messages: List[str] = []
    new_args = forwarded_args[:]
    vol_idx = 0

    def hostdir_to_container(host_dir: Path) -> str:
        nonlocal vol_idx
        host_dir = host_dir.resolve()
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
                        f"        Create it, or re-run with --mkdir to auto-create."
                    )
                    return val
            container_dir = hostdir_to_container(parent.resolve())
            return f"{container_dir}/{host_path.name}"

        if not host_path.exists():
            messages.append(f"[ERROR] Input path not found for {flag}: {host_path}")
            return val

        container_dir = hostdir_to_container(host_path.resolve().parent)
        return f"{container_dir}/{host_path.name}"

    i = 0
    while i < len(new_args):
        tok = new_args[i]

        # --flag value
        if tok in in_flags.union(out_flags) and i + 1 < len(new_args):
            val = new_args[i + 1]
            if _looks_like_flag(val):
                i += 1
                continue
            is_output = tok in out_flags
            new_args[i + 1] = handle_path(tok, val, is_output)
            i += 2
            continue

        # --flag=value
        if tok.startswith("--") and "=" in tok:
            flag, val = tok.split("=", 1)
            if flag in in_flags or flag in out_flags:
                is_output = flag in out_flags
                new_args[i] = f"{flag}={handle_path(flag, val, is_output)}"
            i += 1
            continue

        i += 1

    return new_args, mount_map, messages

def docker_run_opggb(
    image: str,
    script: str,
    forwarded_args: List[str],
    mounts: List[str],
    dry_run: bool,
) -> int:
    cmd: List[str] = ["docker", "run", "--rm"]

    # mount current directory as /work
    cwd = Path.cwd().resolve()
    cmd += ["-v", f"{cwd}:/work", "-w", "/work"]

    for m in mounts:
        cmd += ["-v", m]

    cmd += [image, "python", script] + forwarded_args

    if dry_run:
        print("[DRY-RUN]")
        print(" ".join(shlex.quote(x) for x in cmd))
        return 0

    print("[INFO] Running:")
    print(" ".join(shlex.quote(x) for x in cmd))
    return subprocess.run(cmd).returncode

def docker_run_tubemap(
    image: str,
    input_graph_host: str,
    interactive: bool,
    dry_run: bool,
) -> int:
    """
    Run SequenceTubeMap image directly (no tubemap.py, no docker-in-docker).
    Ensures cleanup on Ctrl+C / exit by stopping + removing the container explicitly.
    """
    graph_path = Path(input_graph_host).expanduser().resolve()
    if not graph_path.exists():
        print(f"[ERROR] Input graph not found: {graph_path}", file=sys.stderr)
        return 2

    graph_dir = graph_path.parent
    graph_name = graph_path.name

    # Unique container name for guaranteed cleanup
    cname = f"minegraph_tubemap_{os.getpid()}_{int(time.time())}"

    cmd: List[str] = ["docker", "run", "--rm", "--name", cname]
    if interactive:
        cmd.append("-it")

    cmd += [
        "-v", f"{graph_dir}:/in:ro",
        "-p", f"{TUBEMAP_HOST_PORT}:{TUBEMAP_CONTAINER_PORT}",
        image,
        "bash", "-lc",
    ]

    inner = f"""
set -euo pipefail
echo "[INFO] Preparing SequenceTubeMap data"
mkdir -p /data
cp "/in/{shlex.quote(graph_name)}" /data/input_graph

./scripts/prepare_vg.sh /data/input_graph

echo "[INFO] Looking for generated .xg file"
XG_FILE=$(find /data -type f -name "*.xg" | head -n 1)
if [ -z "$XG_FILE" ]; then
  echo "[ERROR] No .xg file produced by prepare_vg.sh"
  exit 1
fi

mkdir -p exampleData
cp "$XG_FILE" exampleData/

echo "[INFO] Starting SequenceTubeMap"
echo "[INFO] Open http://localhost:{TUBEMAP_HOST_PORT}"
npm start serve
""".strip()

    cmd.append(inner)

    if dry_run:
        print("[DRY-RUN]")
        print(" ".join(shlex.quote(x) for x in cmd))
        return 0

    print(f"[INFO] Running SequenceTubeMap (container name: {cname})")
    print(" ".join(shlex.quote(x) for x in cmd))

    proc = None
    try:
        proc = subprocess.Popen(cmd)
        return proc.wait()

    except KeyboardInterrupt:
        print("\n[INFO] Ctrl+C received — stopping TubeMap container...", file=sys.stderr)
        return 130

    finally:
        try:
            subprocess.run(["docker", "stop", "-t", "1", cname],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except Exception:
            pass
        try:
            subprocess.run(["docker", "rm", "-f", cname],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except Exception:
            pass
        if proc and proc.poll() is None:
            try:
                proc.terminate()
            except Exception:
                pass

# =========================
# NEW: convert runner
# =========================
def docker_run_convert(
    image: str,
    mode: str,
    input_gfa_host: str,
    output_host: str,
    extra_mounts: List[str],
    mkdir: bool,
    dry_run: bool,
) -> int:
    """
    Convert:
      - fasta: runs /src/gfa2fasta.py <in.gfa> <out.fasta>
      - vg:    runs vg convert -g <in.gfa> -v <out.vg>
    Runs inside opggb image.
    """
    in_path = Path(input_gfa_host).expanduser().resolve()
    if not in_path.exists():
        print(f"[ERROR] Input graph not found: {in_path}", file=sys.stderr)
        return 2

    out_path = Path(output_host).expanduser().resolve()
    if not out_path.parent.exists():
        if mkdir:
            out_path.parent.mkdir(parents=True, exist_ok=True)
        else:
            print(f"[ERROR] Output directory does not exist: {out_path.parent}\n"
                  f"        Create it, or re-run with --mkdir to auto-create.",
                  file=sys.stderr)
            return 2

    # If input and output are in different dirs, mount both
    mounts: List[str] = []
    mounts += extra_mounts

    in_dir = in_path.parent
    out_dir = out_path.parent

    # mount dirs as /data/in and /data/out (avoid conflict)
    mounts += [f"{in_dir}:/data/in"]
    if out_dir != in_dir:
        mounts += [f"{out_dir}:/data/out"]
        out_container = f"/data/out/{out_path.name}"
    else:
        out_container = f"/data/in/{out_path.name}"

    in_container = f"/data/in/{in_path.name}"

    cmd: List[str] = ["docker", "run", "--rm"]

    cwd = Path.cwd().resolve()
    cmd += ["-v", f"{cwd}:/work", "-w", "/work"]

    for m in mounts:
        cmd += ["-v", m]

    if mode == "fasta":
        cmd += [image, "python", GFA2FASTA_SCRIPT, in_container, out_container]
    else:  # mode == "vg"
        cmd += [
            image, "bash", "-lc",
            f"set -euo pipefail; vg convert -g {shlex.quote(in_container)} > {shlex.quote(out_container)}"
        ]

    if dry_run:
        print("[DRY-RUN]")
        print(" ".join(shlex.quote(x) for x in cmd))
        return 0

    print("[INFO] Running:")
    print(" ".join(shlex.quote(x) for x in cmd))
    return subprocess.run(cmd).returncode

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="MineGraph.py",
        description="MineGraph CLI wrapper (Docker-based): construct/extract/convert/tubemap",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=f"""
EXAMPLES
========
Construct:
  MineGraph.py construct -- --data_dir ./mito_Data --output_dir ./out --threads 120

Extract:
  MineGraph.py extract -- subgraph -i graph.gfa -w wanted.txt -o out.gfa

Convert (GFA -> FASTA):
  MineGraph.py convert fasta -i graph.gfa -o graph.fasta

Convert (GFA -> VG):
  MineGraph.py convert vg -i graph.gfa -o graph.vg

TubeMap (fixed port http://localhost:{TUBEMAP_HOST_PORT}):
  MineGraph.py tubemap -i out.gfa
"""
    )

    parser.add_argument("--image", default=OPGGB_IMAGE_DEFAULT,
                        help=f"Image for construct/extract/convert (default: {OPGGB_IMAGE_DEFAULT})")
    parser.add_argument("--tubemap-image", default=TUBEMAP_IMAGE_DEFAULT,
                        help=f"Image for tubemap (default: {TUBEMAP_IMAGE_DEFAULT})")
    parser.add_argument("-m", "--mount", action="append", default=[],
                        help="Extra mounts host:container (repeatable) for construct/extract/convert")
    parser.add_argument("--mkdir", action="store_true",
                        help="Auto-create missing output directories on host (for output flags)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print docker command without running")

    sub = parser.add_subparsers(dest="command", required=True)

    pc = sub.add_parser("construct", help="Run construction pipeline (inside opggb)")
    pc.add_argument("forward_args", nargs=argparse.REMAINDER,
                    help="Args forwarded to /src/construct.py (use `--`)")
    pc.set_defaults(script=CONSTRUCT_SCRIPT)

    pe = sub.add_parser("extract", help="Run subgraph extraction (inside opggb)")
    pe.add_argument("forward_args", nargs=argparse.REMAINDER,
                    help="Args forwarded to /src/extract.py (use `--`)")
    pe.set_defaults(script=EXTRACT_SCRIPT)

    # NEW: convert
    pconv = sub.add_parser("convert", help="Convert graph formats (inside opggb)")
    pconv.add_argument("mode", choices=["fasta", "vg"], help="Conversion mode")
    pconv.add_argument("-i", "--input", required=True, help="Input GFA file on host")
    pconv.add_argument("-o", "--output", required=True, help="Output file on host (.fasta/.vg)")

    pt = sub.add_parser("tubemap", help=f"Explore a graph in SequenceTubeMap (http://localhost:{TUBEMAP_HOST_PORT})")
    pt.add_argument("-i", "--input", required=True, help="Input graph file (.gfa/.vg) on the host")
    pt.add_argument("--interactive", action="store_true", help="Run tubemap docker with -it (debug)")

    return parser

def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    # tubemap runs its own image directly
    if args.command == "tubemap":
        return docker_run_tubemap(
            image=args.tubemap_image,
            input_graph_host=args.input,
            interactive=args.interactive,
            dry_run=args.dry_run,
        )

    # NEW: convert is handled separately (not forwarding to python script)
    if args.command == "convert":
        return docker_run_convert(
            image=args.image,
            mode=args.mode,
            input_gfa_host=args.input,
            output_host=args.output,
            extra_mounts=args.mount,
            mkdir=args.mkdir,
            dry_run=args.dry_run,
        )

    # construct/extract:
    if args.forward_args and args.forward_args[0] == "--":
        args.forward_args = args.forward_args[1:]

    if not args.forward_args:
        parser.print_help()
        return 2

    rewritten_args, mount_map, msgs = collect_paths_and_rewrite_args(
        args.command, args.forward_args, mkdir=args.mkdir
    )
    for m in msgs:
        print(m, file=sys.stderr)
    if any(m.startswith("[ERROR]") for m in msgs):
        return 2

    auto_mounts = [f"{host_dir}:{container_dir}" for host_dir, container_dir in mount_map.items()]

    return docker_run_opggb(
        image=args.image,
        script=args.script,
        forwarded_args=rewritten_args,
        mounts=args.mount + auto_mounts,
        dry_run=args.dry_run,
    )

if __name__ == "__main__":
    raise SystemExit(main())