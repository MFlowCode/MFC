"""
Floating-point stability test suite using Verrou.

Features
--------
A. Stability suite (always)
   N random-rounding samples per case, threshold-based PASS/FAIL.

B. Float proxy (--no-float-proxy to skip)
   One run with --rounding-mode=float — deterministic proxy for
   single-precision sensitivity without recompiling.

C. VPREC precision sweep (--no-vprec to skip)
   One run per mantissa-bit level [52,23,16,10] with
   --backend=vprec --vprec-mode=full; shows where each case breaks.

D. verrou_dd_sym on failure (--no-dd-sym to skip)
   Delta-debug bisection isolates the minimal set of *functions* causing
   instability.

E. verrou_dd_line on failure, after dd_sym (--no-dd-line to skip)
   Further bisects to exact *source lines* within the responsible functions.

Logs are saved to fp-stability-logs/ and uploaded as CI artifacts.

Requires:
  - Verrou-enabled Valgrind at $VERROU_HOME/bin/valgrind
    (default: $HOME/.local/verrou)
  - A serial (no-MPI, no-GPU) simulation binary
  - A serial pre_process binary (to generate initial conditions)

Usage:
  ./mfc.sh fp-stability
  ./mfc.sh fp-stability --no-vprec --no-dd-line
  ./mfc.sh fp-stability --sim-binary PATH --pre-binary PATH
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

# Mantissa-bit levels for the VPREC sweep (C).
# 52 = full double, 23 = single, 16 = half-ish, 10 = ultra-low.
VPREC_MANTISSA_BITS = [52, 23, 16, 10]

# Each case:
#   name         - subdirectory under CASES_DIR
#   description  - human-readable purpose
#   compare      - list of D/ filenames to compare
#   threshold    - max L∞ deviation allowed (conserved-variable units)
#   ill_cond     - known ill-conditioning (empty string = none expected)
CASES = [
    {
        "name": "sod_standard",
        "description": "1-D standard Sod, p_L/p_R=10, ideal gas (well-conditioned baseline)",
        "compare": ["cons.1.00.000050.dat", "cons.3.00.000050.dat"],
        "threshold": 1e-13,
        "ill_cond": "",
    },
    {
        "name": "sod_strong",
        "description": "1-D Sod, p_L/p_R=100,000, ideal gas",
        "compare": ["cons.1.00.000050.dat", "cons.3.00.000050.dat"],
        "threshold": 1e-10,
        "ill_cond": "HLLC xi factor: (s_L - vel_L)/(s_L - s_S) cancels near sonic contact",
    },
    {
        "name": "water_stiffened",
        "description": "1-D water shock, stiffened EOS (pi_inf=4046)",
        "compare": ["cons.1.00.000050.dat", "prim.3.00.000050.dat"],
        "threshold": 1e-8,
        "ill_cond": "Pressure recovery: p=(E-pi_inf)/gamma loses ~4 digits (pi_inf/p_right~40,000) [threshold loosened until reduced-energy (Etilde) scheme is merged]",
    },
    {
        "name": "air_water_interface",
        "description": "1-D air/water isobaric contact (two-fluid, pi_inf=4046)",
        "compare": ["cons.1.00.000050.dat", "cons.4.00.000050.dat", "cons.5.00.000050.dat"],
        "threshold": 1e-10,
        "ill_cond": "Mixed-cell pressure recovery: E-alpha_w*gamma_w*pi_inf cancels when alpha_w<<1",
    },
]


def _find_verrou() -> str:
    verrou_home = os.environ.get("VERROU_HOME", os.path.join(os.path.expanduser("~"), ".local", "verrou"))
    candidate = os.path.join(verrou_home, "bin", "valgrind")
    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    return shutil.which("valgrind") or ""


def _find_binary(name: str) -> str:
    install_dir = os.path.join(os.getcwd(), "build", "install")
    candidates = glob.glob(os.path.join(install_dir, "*", "bin", name))
    return max(candidates, key=os.path.getmtime) if candidates else ""


def _find_dd_sym(verrou_bin: str) -> str:
    c = os.path.join(os.path.dirname(verrou_bin), "verrou_dd_sym")
    return c if os.path.isfile(c) else ""


def _find_dd_line(verrou_bin: str) -> str:
    c = os.path.join(os.path.dirname(verrou_bin), "verrou_dd_line")
    return c if os.path.isfile(c) else ""


def _verrou_pythonpath(verrou_bin: str) -> str:
    """Path that must be on PYTHONPATH for verrou_dd_* imports (valgrind/ subdir)."""
    verrou_home = os.path.dirname(os.path.dirname(verrou_bin))
    matches = glob.glob(os.path.join(verrou_home, "lib", "python*", "site-packages", "valgrind"))
    return matches[0] if matches else ""


def _run_preprocess(pp_bin: str, case_dir: str, work_dir: str):
    shutil.copy2(os.path.join(case_dir, "pre_process.inp"), work_dir)
    with open(os.path.join(work_dir, "pre.log"), "w") as f:
        result = subprocess.run([pp_bin], cwd=work_dir, stdout=f, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        raise MFCException(f"pre_process failed (rc={result.returncode}). See {work_dir}/pre.log")


def _run_simulation_verrou(
    verrou_bin: str,
    sim_bin: str,
    work_dir: str,
    run_dir: str,
    rounding_mode: str = None,
    extra_flags: list = None,
):
    """Copy ICs into a fresh tmpdir, run simulation under verrou, collect D/ output.

    rounding_mode is passed as --rounding-mode=<mode> when not None.
    extra_flags are appended before the binary (e.g. --backend=vprec ...).
    """
    with tempfile.TemporaryDirectory(prefix="mfc-fps-") as tmpdir:
        for fname in ["simulation.inp", "indices.dat", "pre_time_data.dat", "io_time_data.dat"]:
            src = os.path.join(work_dir, fname)
            if os.path.exists(src):
                shutil.copy2(src, tmpdir)
        shutil.copytree(os.path.join(work_dir, "p_all"), os.path.join(tmpdir, "p_all"))
        os.makedirs(os.path.join(tmpdir, "D"))

        log_path = os.path.join(run_dir, "verrou.log")
        cmd = [verrou_bin, "--tool=verrou", "--error-limit=no", f"--log-file={log_path}"]
        if rounding_mode:
            cmd.append(f"--rounding-mode={rounding_mode}")
        cmd.extend(extra_flags or [])
        cmd.append(sim_bin)

        with open(os.path.join(run_dir, "sim.out"), "w") as f:
            result = subprocess.run(cmd, cwd=tmpdir, stdout=f, stderr=subprocess.STDOUT, check=False)

        if result.returncode != 0:
            tag = rounding_mode or "vprec"
            raise MFCException(f"simulation ({tag}) exited {result.returncode}. See {run_dir}/sim.out")

        os.makedirs(run_dir, exist_ok=True)
        for fn in os.listdir(os.path.join(tmpdir, "D")):
            shutil.copy2(os.path.join(tmpdir, "D", fn), run_dir)


def _max_diff_np(ref_dir: str, run_dir: str, compare_files: list) -> float:
    import numpy as np

    total = 0.0
    for fname in compare_files:
        ref_p, run_p = os.path.join(ref_dir, fname), os.path.join(run_dir, fname)
        if not os.path.exists(ref_p) or not os.path.exists(run_p):
            return float("inf")
        ref = np.loadtxt(ref_p)[:, 1]
        run = np.loadtxt(run_p)[:, 1]
        total = max(total, float(np.max(np.abs(ref - run))))
    return total


def _run_float_proxy(case: dict, verrou_bin: str, sim_bin: str, work_dir: str, ref_dir: str) -> float:
    """One run with --rounding-mode=float; returns L∞ deviation from nearest-ref."""
    run_dir = os.path.join(work_dir, "float_proxy")
    os.makedirs(run_dir)
    _run_simulation_verrou(verrou_bin, sim_bin, work_dir, run_dir, rounding_mode="float")
    return _max_diff_np(ref_dir, run_dir, case["compare"])


def _run_vprec_sweep(case: dict, verrou_bin: str, sim_bin: str, work_dir: str, ref_dir: str) -> list:
    """Run at each mantissa-bit level. Returns [(bits, dev), ...]."""
    results = []
    for bits in VPREC_MANTISSA_BITS:
        run_dir = os.path.join(work_dir, f"vprec_{bits}")
        os.makedirs(run_dir)
        flags = [
            "--backend=vprec",
            "--vprec-mode=full",
            f"--vprec-precision-binary64={bits}",
            "--vprec-range-binary64=11",
        ]
        try:
            _run_simulation_verrou(verrou_bin, sim_bin, work_dir, run_dir, extra_flags=flags)
            dev = _max_diff_np(ref_dir, run_dir, case["compare"])
        except MFCException:
            dev = float("inf")
        results.append((bits, dev))
    return results


def _write_dd_run_sh(path: str, verrou_bin: str, sim_bin: str, ic_dir: str):
    """Generate dd_run.sh for verrou_dd_sym / verrou_dd_line.

    verrou_dd_* calls: dd_run.sh RUNDIR
    It sets VERROU_ROUNDING_MODE / VERROU_EXCLUDE / VERROU_SOURCE etc. in the
    environment; the script just invokes verrou and copies D/ output to RUNDIR.
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
    """Generate dd_cmp.py for verrou_dd_sym / verrou_dd_line.

    verrou_dd_* calls: dd_cmp.py REF_DIR RUN_DIR
    Exits 0 (stable) or 1 (unstable) based on threshold.
    """
    content = textwrap.dedent(f"""\
        #!/usr/bin/env python3
        # Generated by mfc.sh fp-stability — do not edit by hand.
        import sys, os, numpy as np

        COMPARE_FILES = {compare_files!r}
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


def _dd_env(verrou_bin: str) -> dict:
    """Environment with PYTHONPATH set for verrou_dd_* imports."""
    py_pkg = _verrou_pythonpath(verrou_bin)
    env = os.environ.copy()
    if py_pkg:
        existing = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = ":".join(filter(None, [py_pkg, existing]))
    return env


def _run_dd_tool(
    dd_bin: str,
    dd_dir: str,
    dd_run_sh: str,
    dd_cmp_py: str,
    env: dict,
    log_name: str,
    summary_subdir: str,
    label: str,
):
    """Generic runner for verrou_dd_sym / verrou_dd_line."""
    log_file = os.path.join(dd_dir, log_name)
    cmd = [dd_bin, "--nruns=10", "--rddmin=d", "--reference-rounding=nearest", dd_run_sh, dd_cmp_py]
    cons.print(f"  [dim]running {label} (--nruns=10 --rddmin=d)...[/dim]")
    with open(log_file, "w") as f:
        result = subprocess.run(cmd, cwd=dd_dir, env=env, stdout=f, stderr=subprocess.STDOUT, check=False)
    if result.returncode == 0:
        summary_path = os.path.join(dd_dir, summary_subdir, "rddmin_summary")
        if os.path.isfile(summary_path):
            cons.print(f"  [bold yellow]{label} result[/bold yellow]:")
            with open(summary_path) as f:
                for line in f:
                    cons.print(f"    {line.rstrip()}")
        else:
            cons.print(f"  [dim]{label} done; see {log_file}[/dim]")
    else:
        cons.print(f"  [bold yellow]{label} exited {result.returncode}[/bold yellow] (see {log_file})")


def _run_dd_sym(case: dict, verrou_bin: str, sim_bin: str, work_dir: str, log_dir: str):
    dd_bin = _find_dd_sym(verrou_bin)
    if not dd_bin:
        cons.print("  [dim]verrou_dd_sym not found; skipping delta-debug[/dim]")
        return

    dd_dir = os.path.join(log_dir, case["name"])
    os.makedirs(dd_dir, exist_ok=True)
    dd_run_sh = os.path.join(dd_dir, "dd_run.sh")
    dd_cmp_py = os.path.join(dd_dir, "dd_cmp.py")
    _write_dd_run_sh(dd_run_sh, verrou_bin, sim_bin, work_dir)
    _write_dd_cmp_py(dd_cmp_py, case["compare"], case["threshold"])
    _run_dd_tool(dd_bin, dd_dir, dd_run_sh, dd_cmp_py, _dd_env(verrou_bin), "dd_sym.log", "dd.sym", "verrou_dd_sym")
    cons.print(f"  [dim]dd_sym logs: {dd_dir}[/dim]")


def _run_dd_line(case: dict, verrou_bin: str, sim_bin: str, work_dir: str, log_dir: str):
    dd_bin = _find_dd_line(verrou_bin)
    if not dd_bin:
        cons.print("  [dim]verrou_dd_line not found; skipping line-level debug[/dim]")
        return

    # Reuse scripts written by dd_sym (or write them now if dd_sym was skipped).
    dd_dir = os.path.join(log_dir, case["name"])
    os.makedirs(dd_dir, exist_ok=True)
    dd_run_sh = os.path.join(dd_dir, "dd_run.sh")
    dd_cmp_py = os.path.join(dd_dir, "dd_cmp.py")
    if not os.path.isfile(dd_run_sh):
        _write_dd_run_sh(dd_run_sh, verrou_bin, sim_bin, work_dir)
        _write_dd_cmp_py(dd_cmp_py, case["compare"], case["threshold"])
    _run_dd_tool(dd_bin, dd_dir, dd_run_sh, dd_cmp_py, _dd_env(verrou_bin), "dd_line.log", "dd.line", "verrou_dd_line")


def _run_case(
    case: dict,
    verrou_bin: str,
    sim_bin: str,
    pp_bin: str,
    n_samples: int,
    log_dir: str,
    run_float: bool,
    run_vprec: bool,
    run_dd_sym: bool,
    run_dd_line: bool,
) -> bool:
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
        _run_simulation_verrou(verrou_bin, sim_bin, work_dir, ref_dir, rounding_mode="nearest")

        # --- A: random-rounding stability samples ---
        max_dev = 0.0
        cons.print(f"  [dim]random-rounding runs (N={n_samples})...[/dim]")
        for i in range(n_samples):
            run_dir = os.path.join(work_dir, f"run_{i:02d}")
            os.makedirs(run_dir)
            _run_simulation_verrou(verrou_bin, sim_bin, work_dir, run_dir, rounding_mode="random")
            max_dev = max(max_dev, _max_diff_np(ref_dir, run_dir, compare))

        passed = max_dev <= threshold
        tag = "[bold green]PASS[/bold green]" if passed else "[bold red]FAIL[/bold red]"
        cons.print(f"  {tag}  max_dev={max_dev:.3e}  threshold={threshold:.0e}")

        # --- B: float proxy ---
        if run_float:
            try:
                fdev = _run_float_proxy(case, verrou_bin, sim_bin, work_dir, ref_dir)
                cons.print(f"  float proxy: dev={fdev:.3e}  (single-precision sensitivity)")
            except MFCException as exc:
                cons.print(f"  [dim]float proxy error: {exc}[/dim]")

        # --- C: VPREC sweep ---
        if run_vprec:
            cons.print("  VPREC precision sweep:")
            vprec_results = _run_vprec_sweep(case, verrou_bin, sim_bin, work_dir, ref_dir)
            labels = {52: "double", 23: "single", 16: "~half", 10: "ultra-low"}
            for bits, dev in vprec_results:
                label = labels.get(bits, "")
                label_str = f" ({label})" if label else ""
                marker = ""
                if dev == float("inf"):
                    marker = "  [red]crashed[/red]"
                elif dev > threshold:
                    marker = "  [red]FAIL[/red]"
                cons.print(f"    {bits:2d} bits{label_str}: dev={dev:.3e}{marker}")

        # --- D/E: delta-debug on failure ---
        if not passed:
            if run_dd_sym:
                try:
                    _run_dd_sym(case, verrou_bin, sim_bin, work_dir, log_dir)
                except Exception as exc:
                    cons.print(f"  [bold yellow]dd_sym error[/bold yellow]: {exc}")
            if run_dd_line:
                try:
                    _run_dd_line(case, verrou_bin, sim_bin, work_dir, log_dir)
                except Exception as exc:
                    cons.print(f"  [bold yellow]dd_line error[/bold yellow]: {exc}")
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
    run_float = not ARG("no_float_proxy")
    run_vprec = not ARG("no_vprec")
    run_dd_sym = not ARG("no_dd_sym")
    run_dd_line = not ARG("no_dd_line")

    log_dir = os.path.join(os.getcwd(), "fp-stability-logs")
    os.makedirs(log_dir, exist_ok=True)

    cons.print()
    cons.print("[bold]MFC Floating-Point Stability Suite[/bold]")
    cons.print(f"  verrou:      {verrou_bin}")
    cons.print(f"  simulation:  {sim_bin}")
    cons.print(f"  pre_process: {pp_bin}")
    cons.print(f"  samples:     {n_samples}")
    features = []
    if run_float:
        features.append("float-proxy")
    if run_vprec:
        features.append("vprec-sweep")
    if run_dd_sym:
        features.append("dd_sym")
    if run_dd_line:
        features.append("dd_line")
    cons.print(f"  features:    {', '.join(features) if features else 'stability only'}")
    cons.print(f"  logs:        {log_dir}")
    cons.print()

    start = time.time()
    results = []
    for case in CASES:
        try:
            ok = _run_case(
                case,
                verrou_bin,
                sim_bin,
                pp_bin,
                n_samples,
                log_dir,
                run_float,
                run_vprec,
                run_dd_sym,
                run_dd_line,
            )
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
        cons.print(f"\n  dd_sym/dd_line logs in: {log_dir}")

    sys.exit(0 if n_fail == 0 else 1)
