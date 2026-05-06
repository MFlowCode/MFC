"""
Floating-point stability test suite using Verrou.

Each test case runs the simulation N times under random IEEE-754 rounding
and measures the max deviation from a nearest-rounding reference.  Cases
that exceed their threshold are reported as FAIL.

Requires:
  - Verrou-enabled Valgrind at $VERROU_HOME/bin/valgrind
    (default: $HOME/.local/verrou)
  - A serial (no-MPI, no-GPU) simulation binary
  - A serial pre_process binary (to generate initial conditions)

Usage:
  ./mfc.sh fp-stability --sim-binary PATH --pre-binary PATH
  ./mfc.sh fp-stability  # auto-discovers newest installed binaries
"""

import glob
import os
import shutil
import subprocess
import sys
import tempfile
import time

from .common import MFC_ROOT_DIR, MFCException
from .printer import cons
from .state import ARG

CASES_DIR = os.path.join(MFC_ROOT_DIR, "tests", "fp_stability", "cases")

# Each case:
#   name         - subdirectory under CASES_DIR
#   description  - human-readable purpose
#   compare      - list of D/ filenames to compare (relative to save step)
#   threshold    - max L∞ deviation allowed (in conserved-variable units)
#   ill_cond     - description of the known ill-conditioning
CASES = [
    {
        "name": "sod_strong",
        "description": "1-D Sod, p_L/p_R=100,000, ideal gas",
        "compare": ["cons.1.00.000005.dat", "cons.3.00.000005.dat"],
        "threshold": 1e-10,
        "ill_cond": "HLLC xi factor: (s_L - vel_L)/(s_L - s_S) cancels near sonic contact",
    },
    {
        "name": "water_stiffened",
        "description": "1-D water shock, stiffened EOS (pi_inf=4046)",
        "compare": ["cons.1.00.000005.dat", "prim.3.00.000005.dat"],
        "threshold": 1e-10,
        "ill_cond": "Pressure recovery: p = (E - pi_inf)/gamma loses ~4 digits (pi_inf/p_right ~ 40,000)",
    },
]


def _find_verrou() -> str:
    verrou_home = os.environ.get("VERROU_HOME", os.path.join(os.path.expanduser("~"), ".local", "verrou"))
    candidate = os.path.join(verrou_home, "bin", "valgrind")
    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    fallback = shutil.which("valgrind")
    return fallback or ""


def _find_binary(name: str) -> str:
    install_dir = os.path.join(os.getcwd(), "build", "install")
    candidates = glob.glob(os.path.join(install_dir, "*", "bin", name))
    if not candidates:
        return ""
    return max(candidates, key=os.path.getmtime)


def _exclude_args(exclude_file: str) -> list:
    """Return --exclude flag if a verrou exclusion file is provided."""
    if exclude_file and os.path.isfile(exclude_file):
        return [f"--exclude={exclude_file}"]
    return []


def _run_preprocess(pp_bin: str, case_dir: str, work_dir: str):
    """Generate initial conditions in work_dir."""
    shutil.copy2(os.path.join(case_dir, "pre_process.inp"), work_dir)
    with open(os.path.join(work_dir, "pre.log"), "w") as f:
        result = subprocess.run(
            [pp_bin],
            cwd=work_dir,
            stdout=f,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode != 0:
        raise MFCException(f"pre_process failed (rc={result.returncode}). See {work_dir}/pre.log")


def _run_simulation_verrou(
    verrou_bin: str,
    sim_bin: str,
    work_dir: str,
    run_dir: str,
    rounding_mode: str,
    exclude_args: list,
):
    """Copy IC files into a fresh tmpdir, run simulation under verrou, collect D/ output."""
    with tempfile.TemporaryDirectory(prefix="mfc-fps-") as tmpdir:
        # Copy static inputs
        for fname in ["simulation.inp", "indices.dat", "pre_time_data.dat", "io_time_data.dat"]:
            src = os.path.join(work_dir, fname)
            if os.path.exists(src):
                shutil.copy2(src, tmpdir)

        # Copy ICs (p_all_ic → p_all)
        shutil.copytree(
            os.path.join(work_dir, "p_all"),
            os.path.join(tmpdir, "p_all"),
        )
        os.makedirs(os.path.join(tmpdir, "D"))

        log_path = os.path.join(run_dir, "verrou.log")
        cmd = [
            verrou_bin,
            "--tool=verrou",
            f"--rounding-mode={rounding_mode}",
            "--error-limit=no",
            f"--log-file={log_path}",
            *exclude_args,
            sim_bin,
        ]

        with open(os.path.join(run_dir, "sim.out"), "w") as f:
            result = subprocess.run(cmd, cwd=tmpdir, stdout=f, stderr=subprocess.STDOUT, check=False)

        if result.returncode != 0:
            raise MFCException(f"simulation (rounding={rounding_mode}) exited {result.returncode}. See {run_dir}/sim.out")

        # Collect D/ outputs
        os.makedirs(run_dir, exist_ok=True)
        d_src = os.path.join(tmpdir, "D")
        for fn in os.listdir(d_src):
            shutil.copy2(os.path.join(d_src, fn), run_dir)


def _max_diff(ref_dir: str, run_dir: str, compare_files: list) -> float:
    try:
        import numpy as np
    except ImportError:
        raise MFCException("numpy is required for fp-stability comparison")

    total = 0.0
    for fname in compare_files:
        ref_path = os.path.join(ref_dir, fname)
        run_path = os.path.join(run_dir, fname)
        if not os.path.exists(ref_path):
            raise MFCException(f"Reference output missing: {ref_path}")
        if not os.path.exists(run_path):
            raise MFCException(f"Run output missing: {run_path}")
        ref = np.loadtxt(ref_path)[:, 1]
        run = np.loadtxt(run_path)[:, 1]
        total = max(total, float(np.max(np.abs(ref - run))))
    return total


def _run_case(case: dict, verrou_bin: str, sim_bin: str, pp_bin: str, n_samples: int) -> bool:
    """Run one stability test case. Returns True if the case passes."""
    name = case["name"]
    threshold = case["threshold"]
    compare = case["compare"]
    case_dir = os.path.join(CASES_DIR, name)

    cons.print(f"[bold]{name}[/bold]: {case['description']}")
    cons.indent()
    cons.print(f"  ill-conditioning: {case['ill_cond']}")
    cons.print(f"  threshold: {threshold:.0e}")

    work_dir = tempfile.mkdtemp(prefix=f"mfc-fps-{name}-")
    try:
        # Generate ICs
        cons.print("  [dim]running pre_process...[/dim]")
        shutil.copy2(os.path.join(case_dir, "simulation.inp"), work_dir)
        _run_preprocess(pp_bin, case_dir, work_dir)

        # Reference run (nearest rounding — deterministic baseline)
        ref_dir = os.path.join(work_dir, "ref")
        os.makedirs(ref_dir)
        cons.print("  [dim]reference run (rounding=nearest)...[/dim]")
        _run_simulation_verrou(verrou_bin, sim_bin, work_dir, ref_dir, "nearest", [])

        # Random-rounding samples
        max_dev = 0.0
        cons.print(f"  [dim]random-rounding runs (N={n_samples})...[/dim]")
        for i in range(n_samples):
            run_dir = os.path.join(work_dir, f"run_{i:02d}")
            os.makedirs(run_dir)
            _run_simulation_verrou(verrou_bin, sim_bin, work_dir, run_dir, "random", [])
            dev = _max_diff(ref_dir, run_dir, compare)
            max_dev = max(max_dev, dev)

        passed = max_dev <= threshold

        tag = "[bold green]PASS[/bold green]" if passed else "[bold red]FAIL[/bold red]"
        cons.print(f"  {tag}  max_dev={max_dev:.3e}  threshold={threshold:.0e}")
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    cons.unindent()
    cons.print()
    return passed


def fp_stability():
    verrou_bin = ARG("verrou_binary") or _find_verrou()
    if not verrou_bin or not os.path.isfile(verrou_bin):
        cons.print("[bold yellow]SKIP[/bold yellow]: verrou not found. Install at $HOME/.local/verrou or set VERROU_HOME.")
        sys.exit(0)

    sim_bin = ARG("sim_binary") or _find_binary("simulation")
    if not sim_bin or not os.path.isfile(sim_bin):
        raise MFCException("simulation binary not found. Build with --debug --no-mpi, or pass --sim-binary.")

    pp_bin = ARG("pre_binary") or _find_binary("pre_process")
    if not pp_bin or not os.path.isfile(pp_bin):
        raise MFCException("pre_process binary not found. Build with --no-mpi, or pass --pre-binary.")

    n_samples = ARG("samples") or 5

    cons.print()
    cons.print("[bold]MFC Floating-Point Stability Suite[/bold]")
    cons.print(f"  verrou:      {verrou_bin}")
    cons.print(f"  simulation:  {sim_bin}")
    cons.print(f"  pre_process: {pp_bin}")
    cons.print(f"  samples:     {n_samples}")
    cons.print()

    start = time.time()
    results = []
    for case in CASES:
        try:
            ok = _run_case(case, verrou_bin, sim_bin, pp_bin, n_samples)
        except MFCException as exc:
            cons.print(f"  [bold red]ERROR[/bold red]: {exc}")
            ok = False
        results.append((case["name"], ok))

    elapsed = time.time() - start
    n_pass = sum(1 for _, ok in results if ok)
    n_fail = len(results) - n_pass

    cons.print(f"[bold]Results[/bold] ({elapsed:.0f}s):  [green]{n_pass} passed[/green]  [red]{n_fail} failed[/red]")
    for name, ok in results:
        mark = "[green]✓[/green]" if ok else "[red]✗[/red]"
        cons.print(f"  {mark} {name}")

    sys.exit(0 if n_fail == 0 else 1)
