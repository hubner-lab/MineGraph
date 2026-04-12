#!/usr/bin/env python3
"""
MineGraph pipeline benchmarking.

Measures per-step resource usage by:
  - Injecting --name <id> into every docker run command
  - Polling docker stats --no-stream in a background thread
  - Recording wall clock, peak/mean CPU%, peak/mean RAM (MiB)

Output: JSON report + console table.
"""

import json
import re
import subprocess
import threading
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Optional


# =============================================================================
# Memory unit converter
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
    """Convert a docker stats memory string (e.g. '1.23GiB', '456MiB') to MiB."""
    m = _MEM_RE.match(s.strip())
    if not m:
        return 0.0
    return float(m.group(1)) * _TO_MIB.get(m.group(2).lower(), 0.0)


# =============================================================================
# Per-step result dataclass
# =============================================================================

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


# =============================================================================
# Profiler
# =============================================================================

class BenchmarkProfiler:
    """
    Profile each MineGraph pipeline step.

    Usage in MineGraph.py:

        profiler = BenchmarkProfiler(output_path="out/benchmark_report.json")
        rc = profiler.run(cmd, step_name="Step 1 — Prepare FASTA", dry_run=False)
        ...
        profiler.finalize()
    """

    def __init__(self, output_path: str, poll_interval: float = 1.0):
        """
        Parameters
        ----------
        output_path   : where to write the JSON report
        poll_interval : seconds between docker stats polls (default 1.0)
        """
        self.output_path = Path(output_path)
        self.poll_interval = poll_interval
        self.steps: List[StepMetrics] = []
        self._wall_start = time.perf_counter()

    # -------------------------------------------------------------------------
    # Public API
    # -------------------------------------------------------------------------

    def run(self, cmd: List[str], step_name: str, dry_run: bool = False) -> int:
        """
        Execute *cmd* (a docker run command) and profile it.

        Automatically injects ``--name <unique_id>`` so docker stats can target
        the container. Returns the container exit code.
        """
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
                        [
                            "docker", "stats", "--no-stream",
                            "--format", "{{.CPUPerc}}\t{{.MemUsage}}",
                            container_name,
                        ],
                        capture_output=True,
                        text=True,
                        timeout=5,
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
                            mem_used = mem_str.split("/")[0].strip()
                            mem = _parse_mem_mib(mem_used)
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
        exit_code = proc.wait()
        stop_event.set()
        poll_thread.join(timeout=5)

        elapsed = time.perf_counter() - t0

        metrics.wall_seconds = round(elapsed, 3)
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
        """Write JSON report and print summary table. Call once after pipeline ends."""
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
# Internal helpers
# =============================================================================

def _inject_container_name(cmd: List[str], name: str) -> List[str]:
    """
    Insert ``--name <name>`` immediately after ``docker run`` in *cmd*.
    Leaves the list unchanged if 'run' is not found.
    """
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
