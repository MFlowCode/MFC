"""
Floating-point stability test suite using Verrou.

Each test case runs the simulation N times under random IEEE-754 rounding
and measures the max deviation from a nearest-rounding reference.  Cases
that exceed their threshold are reported as FAIL.

On failure, verrou_dd_sym is run to identify the minimal set of functions
responsible for the instability.  Logs are saved to fp-stability-logs/.

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
import stat
import subprocess
import sys
import tempfile
import textwrap
import time

from .common import MFC_ROOT_DIR, MFCException
from .printer import cons
from .state import ARG

CASES_DIR = os.path.join(MFC_ROOT_DIR, "tests", "fp_stability", "cases")

# Each case:
#   name         - subdirectory under CASES_DIR
#   description  - human-readable purpose
#   compare      - list of D/ filenames to compare
#   threshold    - max L∞ deviation allowed (in conserved-variable units)
#   ill_cond     - description of the known ill-conditioning (empty = none expected)
CASES = [
    {
        "name": "sod_standard",
        "description": "1-D standard Sod, p_L/p_R=10, ideal gas (well-conditioned baseline)",
        "compare": ["cons.1.00.000005.dat", "cons.3.00.000005.dat"],
        "threshold": 1e-13,
        "ill_cond": "",
    },
    {
        "name": "sod_strong",
        "description": "1-D Sod, p_L/p_R=100,000, ideal gas",
        "compare": ["cons.1.00.000010.dat", "cons.3.00.000010.dat"],
        "threshold": 1e-10,
        "ill_cond": "HLLC xi factor: (s_L - vel_L)/(s_L - s_S) cancels near sonic contact",
    },
    {
        "name": "water_stiffened",
        "description": "1-D water shock, stiffened EOS (pi_inf=4046)",
        "compare": ["cons.1.00.000010.dat", "prim.3.00.000010.dat"],
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


def _find_dd_sym(verrou_bin: str) -> str:
    candidate = os.path.join(os.path.dirname(verrou_bin), "verrou_dd_sym")
    return candidate if os.path.isfile(candidate) else ""


def _verrou_pythonpath(verrou_bin: str) -> str:
    """Return the path that must be on PYTHONPATH for verrou_dd_sym's imports."""
    verrou_home = os.path.dirname(os.path.dirname(verrou_bin))
    matches = glob.glob(os.path.join(verrou_home, "lib", "python*", "site-packages", "valgrind"))
    return matches[0] if matches else ""


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
):
    """Copy IC files into a fresh tmpdir, run simulation under verrou, collect D/ output."""
    with tempfile.TemporaryDirectory(prefix="mfc-fps-") as tmpdir:
        for fname in ["simulation.inp", "indices.dat", "pre_time_data.dat", "io_time_data.dat"]:
            src = os.path.join(work_dir, fname)
            if os.path.exists(src):
                shutil.copy2(src, tmpdir)

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
            sim_bin,
        ]

        with open(os.path.join(run_dir, "sim.out"), "w") as f:
            result = subprocess.run(cmd, cwd=tmpdir, stdout=f, stderr=subprocess.STDOUT, check=False)

        if result.returncode != 0:
            raise MFCException(f"simulation (rounding={rounding_mode}) exited {result.returncode}. See {run_dir}/sim.out")

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


def _write_dd_run_sh(path: str, verrou_bin: str, sim_bin: str, ic_dir: str):
    """Write a dd_run.sh script for verrou_dd_sym.

    verrou_dd_sym calls: dd_run.sh RUNDIR
    It sets VERROU_ROUNDING_MODE (reference) or VERROU_EXCLUDE (sample runs)
    in the environment.  The script must place D/ output files into RUNDIR.
    """
    content = textwrap.dedent(f"""\
        #!/usr/bin/env bash
        # Generated by mfc.sh fp-stability — do not edit by hand.
        VERROU_BIN={verrou_bin!r}
        SIM_BIN={sim_bin!r}
        IC_DIR={ic_dir!r}

        RUNDIR="$1"
        TMPDIR_RUN=$(mktemp -d)
        trap 'rm -rf "$TMPDIR_RUN"' EXIT

        cp -r "$IC_DIR/p_all" "$TMPDIR_RUN/p_all"
        cp "$IC_DIR/simulation.inp" "$TMPDIR_RUN/simulation.inp"
        for fname in indices.dat pre_time_data.dat io_time_data.dat; do
            [ -f "$IC_DIR/$fname" ] && cp "$IC_DIR/$fname" "$TMPDIR_RUN/"
        done
        mkdir -p "$TMPDIR_RUN/D"

        cd "$TMPDIR_RUN"
        "$VERROU_BIN" --tool=verrou --error-limit=no "$SIM_BIN"
        rc=$?

        [ -d "$TMPDIR_RUN/D" ] && cp -a "$TMPDIR_RUN/D/." "$RUNDIR/"
        exit $rc
    """)
    with open(path, "w") as f:
        f.write(content)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _write_dd_cmp_py(path: str, compare_files: list, threshold: float):
    """Write a dd_cmp.py script for verrou_dd_sym.

    verrou_dd_sym calls: dd_cmp.py REF_DIR RUN_DIR
    Exits 0 if max L∞ deviation is within threshold, 1 otherwise.
    """
    files_repr = repr(compare_files)
    content = textwrap.dedent(f"""\
        #!/usr/bin/env python3
        # Generated by mfc.sh fp-stability — do not edit by hand.
        import sys, os, numpy as np

        COMPARE_FILES = {files_repr}
        THRESHOLD = {threshold!r}

        ref_dir, run_dir = sys.argv[1], sys.argv[2]
        max_dev = 0.0
        for fname in COMPARE_FILES:
            ref_p = os.path.join(ref_dir, fname)
            run_p = os.path.join(run_dir, fname)
            if not os.path.exists(ref_p) or not os.path.exists(run_p):
                print(f"MISSING: {{fname}}")
                sys.exit(1)
            ref = np.loadtxt(ref_p)[:, 1]
            run = np.loadtxt(run_p)[:, 1]
            dev = float(np.max(np.abs(ref - run)))
            max_dev = max(max_dev, dev)

        print(f"max_dev={{max_dev:.3e}}  threshold={{THRESHOLD:.0e}}")
        sys.exit(0 if max_dev <= THRESHOLD else 1)
    """)
    with open(path, "w") as f:
        f.write(content)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _run_dd_sym(case: dict, verrou_bin: str, sim_bin: str, work_dir: str, log_dir: str):
    """Run verrou_dd_sym to find the minimal set of symbols causing instability."""
    name = case["name"]
    dd_sym_bin = _find_dd_sym(verrou_bin)
    if not dd_sym_bin:
        cons.print("  [dim]verrou_dd_sym not found; skipping delta-debug[/dim]")
        return

    dd_dir = os.path.join(log_dir, name)
    os.makedirs(dd_dir, exist_ok=True)

    dd_run_sh = os.path.join(dd_dir, "dd_run.sh")
    dd_cmp_py = os.path.join(dd_dir, "dd_cmp.py")
    _write_dd_run_sh(dd_run_sh, verrou_bin, sim_bin, work_dir)
    _write_dd_cmp_py(dd_cmp_py, case["compare"], case["threshold"])

    # PYTHONPATH: verrou_dd_sym imports dd_config etc. from the valgrind/ subpackage
    py_pkg = _verrou_pythonpath(verrou_bin)
    env = os.environ.copy()
    if py_pkg:
        existing = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = ":".join(filter(None, [py_pkg, existing]))

    log_file = os.path.join(dd_dir, "dd_sym.log")
    cmd = [
        dd_sym_bin,
        "--nruns=10",
        "--rddmin=d",
        "--reference-rounding=nearest",
        dd_run_sh,
        dd_cmp_py,
    ]

    cons.print("  [dim]running verrou_dd_sym (--nruns=10 --rddmin=d)...[/dim]")
    with open(log_file, "w") as f:
        result = subprocess.run(cmd, cwd=dd_dir, env=env, stdout=f, stderr=subprocess.STDOUT, check=False)

    if result.returncode == 0:
        summary_path = os.path.join(dd_dir, "dd.sym", "rddmin_summary")
        if os.path.isfile(summary_path):
            cons.print("  [bold yellow]dd_sym result[/bold yellow]:")
            with open(summary_path) as f:
                for line in f:
                    cons.print(f"    {line.rstrip()}")
        else:
            cons.print(f"  [dim]dd_sym done; see {log_file}[/dim]")
    else:
        cons.print(f"  [bold yellow]dd_sym exited {result.returncode}[/bold yellow] (see {log_file})")

    cons.print(f"  [dim]dd_sym logs: {dd_dir}[/dim]")


def _run_case(
    case: dict,
    verrou_bin: str,
    sim_bin: str,
    pp_bin: str,
    n_samples: int,
    log_dir: str,
) -> bool:
    """Run one stability test case.  Returns True if the case passes."""
    name = case["name"]
    threshold = case["threshold"]
    compare = case["compare"]
    case_dir = os.path.join(CASES_DIR, name)

    cons.print(f"[bold]{name}[/bold]: {case['description']}")
    cons.indent()
    if case["ill_cond"]:
        cons.print(f"  ill-conditioning: {case['ill_cond']}")
    cons.print(f"  threshold: {threshold:.0e}")

    work_dir = tempfile.mkdtemp(prefix=f"mfc-fps-{name}-")
    passed = False
    try:
        cons.print("  [dim]running pre_process...[/dim]")
        shutil.copy2(os.path.join(case_dir, "simulation.inp"), work_dir)
        _run_preprocess(pp_bin, case_dir, work_dir)

        ref_dir = os.path.join(work_dir, "ref")
        os.makedirs(ref_dir)
        cons.print("  [dim]reference run (rounding=nearest)...[/dim]")
        _run_simulation_verrou(verrou_bin, sim_bin, work_dir, ref_dir, "nearest")

        max_dev = 0.0
        cons.print(f"  [dim]random-rounding runs (N={n_samples})...[/dim]")
        for i in range(n_samples):
            run_dir = os.path.join(work_dir, f"run_{i:02d}")
            os.makedirs(run_dir)
            _run_simulation_verrou(verrou_bin, sim_bin, work_dir, run_dir, "random")
            dev = _max_diff(ref_dir, run_dir, compare)
            max_dev = max(max_dev, dev)

        passed = max_dev <= threshold
        tag = "[bold green]PASS[/bold green]" if passed else "[bold red]FAIL[/bold red]"
        cons.print(f"  {tag}  max_dev={max_dev:.3e}  threshold={threshold:.0e}")

        if not passed:
            try:
                _run_dd_sym(case, verrou_bin, sim_bin, work_dir, log_dir)
            except Exception as exc:
                cons.print(f"  [bold yellow]dd_sym error[/bold yellow]: {exc}")
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

    log_dir = os.path.join(os.getcwd(), "fp-stability-logs")
    os.makedirs(log_dir, exist_ok=True)

    cons.print()
    cons.print("[bold]MFC Floating-Point Stability Suite[/bold]")
    cons.print(f"  verrou:      {verrou_bin}")
    cons.print(f"  simulation:  {sim_bin}")
    cons.print(f"  pre_process: {pp_bin}")
    cons.print(f"  samples:     {n_samples}")
    cons.print(f"  logs:        {log_dir}")
    cons.print()

    start = time.time()
    results = []
    for case in CASES:
        try:
            ok = _run_case(case, verrou_bin, sim_bin, pp_bin, n_samples, log_dir)
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

    if n_fail > 0:
        cons.print(f"\n  dd_sym logs in: {log_dir}")

    sys.exit(0 if n_fail == 0 else 1)
