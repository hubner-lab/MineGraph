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
import json
import os
import yaml
import platform
import re
import shlex
import subprocess
import sys
import threading
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# =============================================================================
# Benchmark profiler (host-side; polls docker stats per pipeline step)
# =============================================================================

_MEM_RE = re.compile(r"^([\d.]+)\s*(B|kB|KiB|MB|MiB|GB|GiB|TB|TiB)$", re.IGNORECASE)
_TO_MIB = {
    "b":   1 / 1024**2,
    "kb":  1 / 1024,
    "kib": 1 / 1024,
    "mb":  1 / 1.048576,
    "mib": 1.0,
    "gb":  1024 / 1.048576,
    "gib": 1024.0,
    "tb":  1024**2 / 1.048576,
    "tib": 1024**2,
}


def _parse_mem_mib(s: str) -> float:
    m = _MEM_RE.match(s.strip())
    if not m:
        return 0.0
    return float(m.group(1)) * _TO_MIB.get(m.group(2).lower(), 0.0)


@dataclass
class StepMetrics:
    step: str
    wall_seconds: float = 0.0
    cpu_pct_peak: float = 0.0
    cpu_pct_mean: float = 0.0
    mem_mib_peak: float = 0.0
    mem_mib_mean: float = 0.0
    poll_samples: int = 0
    exit_code: int = -1
    error: Optional[str] = None


def _inject_container_name(cmd: List[str], name: str) -> List[str]:
    cmd = list(cmd)
    try:
        run_idx = cmd.index("run") + 1
        cmd.insert(run_idx, "--name")
        cmd.insert(run_idx + 1, name)
    except ValueError:
        pass
    return cmd


def _print_step_inline(m: StepMetrics) -> None:
    status = "OK" if m.exit_code == 0 else f"EXIT={m.exit_code}"
    print(
        f"[BENCH] {m.step:<42}  "
        f"wall={m.wall_seconds:>8.1f}s  "
        f"cpu_peak={m.cpu_pct_peak:>6.1f}%  "
        f"mem_peak={m.mem_mib_peak:>8.1f} MiB  "
        f"[{status}]"
    )


def _print_summary_table(report: dict, output_path: Path) -> None:
    W = 78
    print()
    print("=" * W)
    print("  MINEGRAPH PIPELINE BENCHMARK REPORT")
    print("=" * W)
    print(
        f"{'Step':<44} {'Wall(s)':>8} {'CPU%pk':>7} {'CPU%avg':>8} "
        f"{'MEMpk(MiB)':>11} {'MEMavg(MiB)':>12}"
    )
    print("-" * W)
    for s in report["steps"]:
        status = "" if s["exit_code"] == 0 else f" [EXIT={s['exit_code']}]"
        print(
            f"{s['step'] + status:<44} {s['wall_seconds']:>8.1f}"
            f" {s['cpu_pct_peak']:>7.1f} {s['cpu_pct_mean']:>8.1f}"
            f" {s['mem_mib_peak']:>11.1f} {s['mem_mib_mean']:>12.1f}"
        )
    print("-" * W)
    print(f"{'TOTAL WALL TIME':<44} {report['pipeline_wall_seconds']:>8.1f}")
    print("=" * W)
    print(f"[BENCH] Full JSON report → {output_path}")
    print()


def _emergency_stop(proc: Optional[subprocess.Popen], container_name: Optional[str]) -> None:
    """
    Best-effort cleanup on KeyboardInterrupt or unexpected exit.

    Sends ``docker stop`` (graceful, 5 s timeout) followed by ``docker rm -f``
    to ensure the container is gone even if the stop times out.  Then
    terminates the host-side docker-client process if it is still running.
    Called from both ``_run()`` and ``BenchmarkProfiler.run()``.
    """
    print("\n[INTERRUPT] Stopping container and cleaning up...", file=sys.stderr)
    if container_name:
        for cleanup in (
            ["docker", "stop", "-t", "5", container_name],
            ["docker", "rm",   "-f",     container_name],
        ):
            subprocess.run(cleanup, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if proc and proc.poll() is None:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()


class BenchmarkProfiler:
    def __init__(self, output_path: str, poll_interval: float = 1.0):
        self.output_path = Path(output_path)
        self.poll_interval = poll_interval
        self.steps: List[StepMetrics] = []
        self._wall_start = time.perf_counter()

    def run(self, cmd: List[str], step_name: str, dry_run: bool = False) -> int:
        print(f"[BENCH] Starting: {step_name}")
        if dry_run:
            print(f"[DRY-RUN][BENCH] skipped profiling: {step_name}")
            self.steps.append(StepMetrics(step=step_name, exit_code=0))
            return 0

        container_name = f"mg_bench_{int(time.time() * 1000) % 10 ** 9}"
        cmd = _inject_container_name(cmd, container_name)

        metrics = StepMetrics(step=step_name)
        cpu_samples: List[float] = []
        mem_samples: List[float] = []
        stop_event = threading.Event()

        def _poll_loop():
            while not stop_event.is_set():
                try:
                    result = subprocess.run(
                        ["docker", "stats", "--no-stream",
                         "--format", "{{.CPUPerc}}\t{{.MemUsage}}",
                         container_name],
                        capture_output=True, text=True, timeout=5,
                    )
                    if result.returncode == 0 and result.stdout.strip():
                        for line in result.stdout.strip().splitlines():
                            parts = line.split("\t")
                            if len(parts) != 2:
                                continue
                            cpu_str, mem_str = parts
                            try:
                                cpu = float(cpu_str.strip().rstrip("%") or "0")
                            except ValueError:
                                cpu = 0.0
                            mem = _parse_mem_mib(mem_str.split("/")[0].strip())
                            cpu_samples.append(cpu)
                            mem_samples.append(mem)
                except Exception:
                    pass
                stop_event.wait(self.poll_interval)

        poll_thread = threading.Thread(target=_poll_loop, daemon=True)
        t0 = time.perf_counter()
        try:
            proc = subprocess.Popen(cmd)
        except Exception as exc:
            metrics.error = str(exc)
            metrics.exit_code = -1
            self.steps.append(metrics)
            print(f"[BENCH][ERROR] Could not start container: {exc}")
            return -1

        poll_thread.start()
        try:
            exit_code = proc.wait()
        except KeyboardInterrupt:
            stop_event.set()
            poll_thread.join(timeout=2)
            _emergency_stop(proc, container_name)
            raise
        stop_event.set()
        poll_thread.join(timeout=5)

        metrics.wall_seconds = round(time.perf_counter() - t0, 3)
        metrics.exit_code = exit_code
        metrics.poll_samples = len(cpu_samples)
        if cpu_samples:
            metrics.cpu_pct_peak = round(max(cpu_samples), 2)
            metrics.cpu_pct_mean = round(sum(cpu_samples) / len(cpu_samples), 2)
        if mem_samples:
            metrics.mem_mib_peak = round(max(mem_samples), 2)
            metrics.mem_mib_mean = round(sum(mem_samples) / len(mem_samples), 2)

        self.steps.append(metrics)
        _print_step_inline(metrics)
        return exit_code

    def finalize(self) -> dict:
        total_wall = round(time.perf_counter() - self._wall_start, 3)
        report = {
            "pipeline_wall_seconds": total_wall,
            "steps": [asdict(m) for m in self.steps],
        }
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.output_path, "w") as fh:
            json.dump(report, fh, indent=2)
        _print_summary_table(report, self.output_path)
        return report

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
_INJECT_SCRIPT       = "/src/inject.py"

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

    # Inject a unique container name into docker-run commands so we can
    # stop the container on interrupt even though the process handle is local.
    cname    = None
    run_cmd  = list(cmd)
    if len(cmd) >= 2 and cmd[0] == "docker" and cmd[1] == "run":
        cname   = f"mg_{os.getpid()}_{int(time.time() * 1000) % 10 ** 9}"
        run_cmd = _inject_container_name(run_cmd, cname)

    proc = subprocess.Popen(run_cmd)
    try:
        return proc.wait()
    except KeyboardInterrupt:
        _emergency_stop(proc, cname)
        raise


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
    profiler: Optional["BenchmarkProfiler"] = None,
) -> int:
    """
    Run one script inside the opggb container.

    Mounts:
      <data_dir_host>   → /data      (read-only input FASTA directory)
      <output_dir_host> → /output    (read-write; matches hardcoded paths in scripts)

    When *profiler* is provided the step is measured (wall time, CPU, RAM).
    """
    # output_dir_host may not exist yet — create first, then resolve strictly
    Path(output_dir_host).expanduser().mkdir(parents=True, exist_ok=True)
    abs_data = str(Path(data_dir_host).expanduser().resolve(strict=False))
    abs_out  = str(Path(output_dir_host).expanduser().resolve())

    cmd = (
        ["docker", "run", "--rm", "--init"]
        + _platform_flags()
        + ["-v", f"{abs_data}:/data",
           "-v", f"{abs_out}:/output"]
    )
    for m in extra_mounts:
        cmd += ["-v", m]
    cmd += [image, "python", script] + script_args

    if profiler is not None and label:
        return profiler.run(cmd, step_name=label, dry_run=dry_run)
    return _run(cmd, dry_run, label)



# =============================================================================
# Resume checkpoint detection (host-side, no Docker needed)
# =============================================================================

def _detect_resume_point(abs_out: str) -> str:
    """
    Inspect *abs_out* on the host and return the earliest step that still needs
    to run.  Return values:

      'all'          – nothing done; run every step
      'tr'           – FASTA ready (Step 1 done); run Steps 2-3-4-5
      'pggb'         – FASTA + TR ready (Steps 1-3 done); run full PGGB + stats
      'pggb_resume'  – pggb_output/ exists but incomplete; run PGGB --resume + stats
      'stats'         – PGGB complete; no Step 5 outputs yet; run full Step 5
      'stats_partial' – PGGB complete; some Step 5 outputs exist; run missing sub-steps
      'done'          – all steps complete; nothing to do
    """
    out = Path(abs_out)

    # ── Step 1 ──────────────────────────────────────────────────────────────
    fasta_gz    = out / "panSN_output.fasta.gz"
    downsampled = out / "downsampled_panSN_output.fasta"
    params_yaml = out / "params.yaml"

    if not (fasta_gz.exists() and downsampled.exists() and params_yaml.exists()):
        return "all"

    # ── Steps 2 + 3 ─────────────────────────────────────────────────────────
    try:
        with open(params_yaml) as _f:
            _params = yaml.safe_load(_f) or {}
        step2_3_done = "segment_length" in _params
    except Exception:
        step2_3_done = False

    if not step2_3_done:
        return "tr"

    # ── Step 4 ──────────────────────────────────────────────────────────────
    pggb_dir = out / "pggb_output"
    if not pggb_dir.exists():
        return "pggb"

    gfa_done = any(pggb_dir.glob("*.smooth.final.gfa"))
    paf_done = any(pggb_dir.glob("*.alignments.wfmash.paf"))
    maf_done = any(pggb_dir.glob("*.smooth.maf"))

    if not (gfa_done and paf_done and maf_done):
        return "pggb_resume"

    # ── Step 5 ──────────────────────────────────────────────────────────────
    # "done" only when the two heaviest outputs both exist (XLSX + tree).
    # If either is missing we fall through to "stats" and let
    # _detect_stats_substeps() choose which sub-steps to skip.
    xlsx  = out / "MineGraph_output" / "statistics"  / "graph_stats.xlsx"
    odgi  = out / "MineGraph_output" / "phylogenetics_msa" / "odgi_matrix.tsv"
    tree  = out / "MineGraph_output" / "phylogenetics_msa" / "graph_phylo_tree.nwk"
    tree_png = out / "MineGraph_output" / "phylogenetics_msa" / "graph_phylo_tree.png"

    if xlsx.exists() and odgi.exists() and (tree.exists() or tree_png.exists()):
        return "done"
    if xlsx.exists():
        return "stats_partial"

    return "stats"


def _detect_stats_substeps(abs_out: str) -> List[str]:
    """
    Inspect existing Step 5 output files and return the list of --skip-*
    flags to forward to run_stats.py so only missing sub-steps are executed.

    Sub-step mapping (sentinel → flag):
      graph_stats.xlsx                    → --skip-core-stats
      plots/node_histogram_by_paths.png   → --skip-plots
      phylogenetics_msa/odgi_matrix.tsv   → --skip-odgi
      phylogenetics_msa/graph_phylo_tree.nwk (or .png) → --skip-tree
    """
    out   = Path(abs_out)
    flags: List[str] = []

    xlsx     = out / "MineGraph_output" / "statistics"         / "graph_stats.xlsx"
    hist     = out / "MineGraph_output" / "plots"              / "node_histogram_by_paths.png"
    odgi     = out / "MineGraph_output" / "phylogenetics_msa"  / "odgi_matrix.tsv"
    tree_nwk = out / "MineGraph_output" / "phylogenetics_msa"  / "graph_phylo_tree.nwk"
    tree_png = out / "MineGraph_output" / "phylogenetics_msa"  / "graph_phylo_tree.png"

    if xlsx.exists():
        flags.append("--skip-core-stats")
        print(f"[RESUME]   Sub-step 5a (core stats)   — found {xlsx.name}, skipping.")
    if hist.exists():
        flags.append("--skip-plots")
        print(f"[RESUME]   Sub-step 5b (plots)         — found {hist.name}, skipping.")
    if odgi.exists():
        flags.append("--skip-odgi")
        print(f"[RESUME]   Sub-step 5c (ODGI matrix)   — found {odgi.name}, skipping.")
    if tree_nwk.exists() or tree_png.exists():
        flags.append("--skip-tree")
        found = tree_nwk.name if tree_nwk.exists() else tree_png.name
        print(f"[RESUME]   Sub-step 5d (phylo tree)    — found {found}, skipping.")

    return flags


# =============================================================================
# construct — five-step host-orchestrated pipeline
# =============================================================================

def run_construct(args, profiler: Optional["BenchmarkProfiler"] = None) -> int:
    """
    Orchestrate the full construct pipeline from the host.
    Each step is a separate docker run — no docker-in-docker.

    When *profiler* is provided every docker step is measured.
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
    cp.add_argument("--tree_bs",     default="1000",
                        help="Bootstrap replicates (graph tree: PAV Felsenstein bootstrap; "
                             "MSA tree: RAxML --bs-trees). Default: 1000.")
    cp.add_argument("--quantile",    default="25")
    cp.add_argument("--top_n",       default="50")
    cp.add_argument("--window_size", default="1000")
    cp.add_argument("--mode",        default="all")
    cp.add_argument("--tree_type",   default="graph")
    cp.add_argument("--sc_th",       default="5")
    cp.add_argument("--phyl-tree",   action="store_true", dest="phyl_tree",  default=True)
    cp.add_argument("--plots",       action="store_true",                    default=True)
    cp.add_argument("--only-stats",  action="store_true", dest="only_stats", default=False)
    cp.add_argument("--resume",      action="store_true", default=False,
                        help="Resume from the last completed step")
    cp.add_argument("--convert-gfa", action="store_true", dest="convert_gfa", default=False,
                        help="Convert the final GFA to VG and FASTA (default: off).")
    cp.add_argument("--bootstrap-method", dest="bootstrap_method", default="felsenstein",
                        choices=["felsenstein", "ufboot"],
                        help="Bootstrap algorithm for the graph-based UPGMA tree: "
                             "'felsenstein' (default) or 'ufboot' (UFBoot-PAV adaptive).")
    cp.add_argument("--ufboot-min-rep", type=int, dest="ufboot_min_rep", default=1000,
                        help="UFBoot: minimum replicates before convergence testing (default: 1000).")
    cp.add_argument("--ufboot-max-rep", type=int, dest="ufboot_max_rep", default=10000,
                        help="UFBoot: hard cap on total replicates (default: 10000).")
    cp.add_argument("--ufboot-convergence", type=float, dest="ufboot_convergence", default=0.99,
                        help="UFBoot: Pearson r convergence threshold (default: 0.99).")
    cp.add_argument("--ufboot-batch", type=int, dest="ufboot_batch", default=100,
                        help="UFBoot: replicates per convergence-check round (default: 100).")
    cp.add_argument("--species",     default="viridiplantae",
                        help="Species/clade name passed to RepeatMasker -species "
                             "(default: viridiplantae — highest available plant clade). "
                             "Examples: embryophyta, magnoliophyta, monocotyledons, poaceae.")
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
    rm_species  = ca.species

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
    # Resume: detect which step to start from
    # ------------------------------------------------------------------
    checkpoint = _detect_resume_point(abs_out) if ca.resume else "all"
    if ca.resume:
        if checkpoint == "done":
            print("[RESUME] All steps already complete. Nothing to do.")
            return 0
        print(f"[RESUME] Detected checkpoint: '{checkpoint}' — resuming from there.")

    skip_fasta    = ca.resume and checkpoint not in ("all",)
    skip_tr       = ca.resume and checkpoint not in ("all", "tr")
    skip_pggb     = ca.resume and checkpoint in ("stats", "stats_partial", "done")
    pggb_resume   = ca.resume and checkpoint == "pggb_resume"
    skip_stats    = ca.resume and checkpoint == "done"

    # Intra-step 5 resume: detect which sub-steps already have output files
    stats_skip_flags: List[str] = []
    if ca.resume and checkpoint in ("stats_partial",):
        print("[RESUME] Step 5 partially complete — detecting finished sub-steps...")
        stats_skip_flags = _detect_stats_substeps(abs_out)

    # ------------------------------------------------------------------
    # STEP 1 – Prepare FASTA input (opggb container)
    # ------------------------------------------------------------------
    if mode in ("all", "extract-tr", "construct-graph"):
        if skip_fasta:
            print("[RESUME] Step 1 (FASTA prep) — already done, skipping.")
        else:
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
                profiler=profiler,
            )
            if rc != 0:
                return rc

        # ------------------------------------------------------------------
        # STEP 2 – RepeatMasker
        # ------------------------------------------------------------------
        if skip_tr:
            print("[RESUME] Steps 2+3 (RepeatMasker + TR extract) — already done, skipping.")
        else:
            rm_label = "Step 2 — RepeatMasker"
            print(f"[INFO] {rm_label}")
            rm_cmd = (
                ["docker", "run", "--rm"]
                + _platform_flags()
                + [
                    "-v", f"{abs_out}:/data",
                    REPEATMASKER_IMAGE,
                    "bash", "-c",
                    f"RepeatMasker /data/downsampled_panSN_output.fasta "
                    f"-no_is -pa {threads} -s -gff -species {rm_species}",
                ]
            )
            if profiler is not None:
                rc = profiler.run(rm_cmd, step_name=rm_label, dry_run=dry_run)
            else:
                rc = _run(rm_cmd, dry_run, "")
            if rc != 0:
                return rc

            # ------------------------------------------------------------------
            # STEP 3 – Parse RepeatMasker output → update params.yaml
            # ------------------------------------------------------------------
            rc = _opggb(
                image=image,
                script=_EXTRACT_TR_SCRIPT,
                script_args=[],
                data_dir_host=abs_out,
                output_dir_host=abs_out,
                extra_mounts=[],
                dry_run=dry_run,
                label="Step 3 — Extract longest TR / update params.yaml",
                profiler=profiler,
            )
            if rc != 0:
                return rc

    if mode == "extract-tr":
        print("[INFO] Extract-only mode complete.")
        return 0

    # ------------------------------------------------------------------
    # STEP 4 – PGGB graph construction (opggb)
    # ------------------------------------------------------------------
    if skip_pggb:
        print("[RESUME] Step 4 (PGGB) — already done, skipping.")
    else:
        pggb_args = [threads]
        if pggb_resume:
            pggb_args.append("--resume")
            print("[RESUME] Step 4 (PGGB) — partial output found, running with --resume.")
        rc = _opggb(
            image=image,
            script=_PGGB_SCRIPT,
            script_args=pggb_args,
            data_dir_host=abs_out,
            output_dir_host=abs_out,
            extra_mounts=[],
            dry_run=dry_run,
            label="Step 4 — PGGB graph construction",
            profiler=profiler,
        )
        if rc != 0:
            return rc

    if mode == "construct-graph":
        print("[INFO] Graph construction complete. Skipping statistics.")
        return 0

    # ------------------------------------------------------------------
    # STEP 5 – Graph statistics (opggb)
    # ------------------------------------------------------------------
    if skip_stats:
        print("[RESUME] Step 5 (stats) — already done, skipping.")
        return 0
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
    if ca.phyl_tree:    stats_args.append("--phyl-tree")
    if ca.plots:        stats_args.append("--plots")
    if ca.only_stats:   stats_args.append("--only-stats")
    if ca.convert_gfa:  stats_args.append("--convert-gfa")
    stats_args += [
        "--bootstrap-method", ca.bootstrap_method,
        "--ufboot-min-rep",   str(ca.ufboot_min_rep),
        "--ufboot-max-rep",   str(ca.ufboot_max_rep),
        "--ufboot-convergence", str(ca.ufboot_convergence),
        "--ufboot-batch",     str(ca.ufboot_batch),
    ]
    # Append intra-step resume skip flags (empty list when not resuming)
    stats_args += stats_skip_flags

    rc = _opggb(
        image=image,
        script=_STATS_SCRIPT,
        script_args=stats_args,
        data_dir_host=abs_out,
        output_dir_host=abs_out,
        extra_mounts=[],
        dry_run=dry_run,
        label="Step 5 — Graph statistics",
        profiler=profiler,
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
# inject — gene annotation injection + gggenes visualisation
# =============================================================================

def run_inject(args) -> int:
    """
    Mount all input files into the container and run inject.py.

    Each unique parent directory is mounted to a separate /data/volN path so
    the user can keep graph, BED, and genomes files in different directories.
    The output directory is mounted read-write at its own /data/volN path.
    """
    # ── Collect path-bearing arguments ───────────────────────────────────
    path_flags: List[Tuple[str, str, bool]] = [
        ("--graph",    args.graph,   False),
        ("--bed",      args.bed,     False),
        ("--genomes",  args.genomes, False),
        ("--output",   args.output,  True),
    ]
    if args.label_map:
        path_flags.append(("--label-map", args.label_map, False))
    if args.color_map:
        path_flags.append(("--color-map", args.color_map, False))

    # ── Build mount map: host_dir → /data/volN ────────────────────────────
    mount_map: Dict[Path, str] = {}
    vol_idx = 0
    errors: List[str] = []

    def register(host_str: str, is_output: bool) -> str:
        nonlocal vol_idx
        p = Path(host_str).expanduser()
        if is_output:
            p.mkdir(parents=True, exist_ok=True)
            parent = p.resolve(strict=False)
        else:
            p = p.resolve(strict=False)
            if not p.exists():
                errors.append(f"[ERROR] Input not found: {p}")
                return str(p)
            parent = p.parent
        if parent not in mount_map:
            mount_map[parent] = f"/data/vol{vol_idx}"
            vol_idx += 1
        if is_output:
            return mount_map[parent]
        return f"{mount_map[parent]}/{p.name}"

    # Build container-path rewrites for each flag
    rewrites: Dict[str, str] = {}
    for flag, host_val, is_out in path_flags:
        rewrites[flag] = register(host_val, is_out)

    for err in errors:
        print(err, file=sys.stderr)
    if errors:
        return 2

    # ── Build forwarded args list with container paths ────────────────────
    container_args: List[str] = [
        "--graph",   rewrites["--graph"],
        "--bed",     rewrites["--bed"],
        "--genomes", rewrites["--genomes"],
        "--output",  rewrites["--output"],
        "--threads", args.threads,
        "--bins",    args.bins,
        "--gap",     args.gap,
    ]
    if args.anchor:
        container_args += ["--anchor", args.anchor]
    if args.label_map:
        container_args += ["--label-map", rewrites["--label-map"]]
    if args.color_map:
        container_args += ["--color-map", rewrites["--color-map"]]
    if args.no_subgraph:
        container_args.append("--no-subgraph")
    if args.skip_viz:
        container_args.append("--skip-viz")

    # ── Docker command ────────────────────────────────────────────────────
    auto_mounts = [f"{h}:{c}" for h, c in mount_map.items()]
    return _opggb_simple(
        args.image,
        _INJECT_SCRIPT,
        container_args,
        args.mount + auto_mounts,
        args.dry_run,
    )


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

  # Use monocot-specific RepeatMasker database (better sensitivity for Poales):
  MineGraph.py construct -- --data_dir ./mito_Data --output_dir ./out \\
    --metadata ./samples.csv --threads 32 --phyl-tree --plots --species monocotyledons

Resume interrupted run:
  MineGraph.py construct -- --data_dir ./mito_Data --output_dir ./out \\
    --metadata ./samples.csv --threads 32 --phyl-tree --plots --resume

Extract subgraph:
  MineGraph.py extract -- subgraph -i graph.gfa -w wanted.txt -o sub.gfa

Inject gene annotations and visualise:
  MineGraph.py inject --graph graph.og --bed genes.bed \\
    --genomes reps.txt --output ./inject_out --threads 8

  # With custom molecule labels and gene colours:
  MineGraph.py inject --graph graph.og --bed genes.bed \\
    --genomes reps.txt --output ./inject_out \\
    --label-map labels.tsv --color-map colors.tsv

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
    parser.add_argument("--benchmark",     action="store_true",
                        help="Profile each pipeline step: wall time, CPU%%, and RAM (MiB). "
                             "Writes benchmark_report.json to the output directory.")

    sub = parser.add_subparsers(dest="command", required=True)

    # construct
    pc = sub.add_parser("construct", help="Run full pangenome construction pipeline")
    pc.add_argument("forward_args", nargs=argparse.REMAINDER,
                    help="Arguments forwarded to construct steps (use '--' separator)")

    # extract
    pe = sub.add_parser("extract", help="Subgraph extraction / graph inspection")
    pe.add_argument("forward_args", nargs=argparse.REMAINDER,
                    help="Arguments forwarded to extract.py (use '--' separator)")

    # inject
    pi = sub.add_parser(
        "inject",
        help="Inject gene BED annotations into graph and produce gggenes visualisation",
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
inject — Gene annotation injection and gggenes visualisation.

Embeds gene annotations from a BED file into a pangenome graph as traversal
paths, then resolves their positions in each representative genome using
odgi untangle, and outputs gene-arrow plots and bin-coverage heatmaps.

BED format (tab-separated):
  path_name  start  end  gene_name  score  strand
  - path_name must match a path already present in the graph.
  - gene_name (col 4) becomes the injected path name.

Outputs written to --output/:
  inject_gene_arrows.png / .pdf   — gene-arrow visualisation
  inject_bin_coverage.png / .pdf  — per-path bin-coverage heatmap
  inject_gene_presence.tsv        — gene copy counts per genome
  untangle.tsv                    — raw odgi untangle positions
  bin_coverage.tsv                — raw odgi bin coverage
  inject_stats.txt                — graph statistics summary
""",
    )
    pi.add_argument("--graph",    "-g", required=True,
                    help="Input pangenome graph (.gfa, .gfa.gz, or .og)")
    pi.add_argument("--bed",      "-b", required=True,
                    help="Gene annotations in BED6 format.\n"
                         "Column 1 must be a path name present in the graph.")
    pi.add_argument("--genomes",  "-r", required=True,
                    help="Text file: one representative genome path name per line.")
    pi.add_argument("--output",   "-o", required=True,
                    help="Output directory (created if it does not exist).")
    pi.add_argument("--threads",  "-t", default="4",
                    help="Number of threads (default: 4)")
    pi.add_argument("--anchor",   default=None,
                    help="Anchor path for subgraph extraction.\n"
                         "Default: first line of --genomes.")
    pi.add_argument("--bins",     default="500",
                    help="Number of bins for coverage plot (default: 500).")
    pi.add_argument("--gap",      default="5000",
                    help="Gap threshold (bp) for merging fragmented gene\n"
                         "segments within the same copy (default: 5000).")
    pi.add_argument("--label-map", dest="label_map", default=None,
                    help="Optional TSV (path_prefix<TAB>display_label).")
    pi.add_argument("--color-map", dest="color_map", default=None,
                    help="Optional TSV (gene_name<TAB>hex_color).")
    pi.add_argument("--no-subgraph", action="store_true", dest="no_subgraph",
                    help="Skip subgraph extraction; run on full injected graph.")
    pi.add_argument("--skip-viz",    action="store_true", dest="skip_viz",
                    help="Skip visualisation; output TSV files only.")

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
        profiler = None
        if getattr(args, "benchmark", False):
            if BenchmarkProfiler is None:
                print("[WARN] benchmark.py not found — running without profiling.", file=sys.stderr)
            else:
                # Derive output_dir from forwarded args to know where to write the report
                import argparse as _ap
                _cp = _ap.ArgumentParser(add_help=False)
                _cp.add_argument("--output_dir", default=".")
                _fwd = args.forward_args[1:] if args.forward_args and args.forward_args[0] == "--" else args.forward_args
                _ca, _ = _cp.parse_known_args(_fwd)
                _report_path = str(Path(_ca.output_dir).expanduser() / "benchmark_report.json")
                profiler = BenchmarkProfiler(output_path=_report_path, poll_interval=1.0)
                print(f"[BENCH] Benchmark mode ON — report will be written to: {_report_path}")

        rc = run_construct(args, profiler=profiler)

        if profiler is not None:
            profiler.finalize()

        return rc
    if args.command == "extract":
        return run_extract(args)
    if args.command == "inject":
        return run_inject(args)
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
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\n[INFO] Pipeline interrupted by user.", file=sys.stderr)
        raise SystemExit(130)
