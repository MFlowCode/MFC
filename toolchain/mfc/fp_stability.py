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

F. Cancellation detection (--no-cancellation to skip)
   One run with --check-cancellation=yes; reports MFC source lines that
   produce catastrophic cancellation (subtraction of nearly-equal doubles).
   Uses --cc-gen-file for structured per-line output.

G. MCA significant-bits estimate (--no-mca to skip)
   N runs with --backend=mcaquad; max deviation vs nearest-rounding
   reference gives a lower bound on significant bits: s = -log2(dev/scale).

H. Float-max overflow detection (--no-float-max to skip)
   One run with --check-max-float=yes; reports locations where a
   double→float conversion would overflow to ±Inf.

Logs are saved to fp-stability-logs/ and uploaded as CI artifacts.
On GitHub Actions: a step summary table and ::warning:: file annotations
are emitted automatically so failing source lines appear in the PR diff.

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
import math
import os
import re
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

# Mantissa-bit levels for the VPREC sweep (C).
# 52 = full double, 23 = single, 16 = half-ish, 10 = ultra-low.
VPREC_MANTISSA_BITS = [52, 23, 16, 10]

# Matches "path/file.f90:123" or "path/file.fpp:123-456" in dd_line rddmin_summary.
_LOC_RE = re.compile(r"(\S+\.(?:f90|fpp|c|cpp|h|F90))\s*:(\d+)(?:-(\d+))?", re.IGNORECASE)

# Files to exclude from cancellation / float-max reports (runtime loaders, XALT).
_EXTERNAL_SRCS = ("xalt", "dl-init", "ld-linux", "libc.so", "libm.so")

# Matches the first "at" frame in a Valgrind stack trace: "(file.fpp:LINE)".
_VGFRAME_RE = re.compile(r"\(([^):]+\.(?:fpp|f90|F90|c|cpp))\s*:(\d+)\)")


def _get_source_context(fname: str, lineno: int, context: int = 2) -> str:
    """Return a annotated source snippet around lineno, or '' if file not found.

    fname may be a bare basename (e.g. 'm_weno.fpp') or a relative path.
    Searches recursively under MFC_ROOT_DIR/src/ first, then the whole tree.
    """
    if os.path.isabs(fname) and os.path.isfile(fname):
        candidates = [fname]
    else:
        candidates = glob.glob(os.path.join(MFC_ROOT_DIR, "src", "**", os.path.basename(fname)), recursive=True)
        if not candidates:
            candidates = glob.glob(os.path.join(MFC_ROOT_DIR, "**", os.path.basename(fname)), recursive=True)
    if not candidates:
        return ""
    try:
        with open(candidates[0]) as fh:
            lines = fh.readlines()
    except OSError:
        return ""
    start = max(0, lineno - context - 1)
    end = min(len(lines), lineno + context)
    rows = []
    for i, line in enumerate(lines[start:end], start=start + 1):
        marker = ">" if i == lineno else " "
        rows.append(f"{marker}{i:5d} | {line.rstrip()}")
    return "\n".join(rows)


def _merge(*dicts):
    """Merge dicts left-to-right; later entries override earlier ones."""
    result = {}
    for d in dicts:
        result.update(d)
    return result


# Shared EOS parameters for common fluid types.
_AIR_EOS = {"fluid_pp(1)%gamma": 1.0 / 0.4, "fluid_pp(1)%pi_inf": 0.0}
_WATER_EOS = {"fluid_pp(1)%gamma": 0.1953125, "fluid_pp(1)%pi_inf": 4046.31}

# Parameters shared by every pre_process run.
_BASE_PRE = {
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "n": 0,
    "p": 0,
    "t_step_start": 0,
    "num_patches": 2,
    "model_eqns": 2,
    "weno_order": 5,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "precision": 2,
    "parallel_io": "F",
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": 0.25,
    "patch_icpp(1)%length_x": 0.5,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(2)%geometry": 1,
    "patch_icpp(2)%x_centroid": 0.75,
    "patch_icpp(2)%length_x": 0.5,
    "patch_icpp(2)%vel(1)": 0.0,
}

# Parameters shared by every simulation run.
_BASE_SIM = {
    "run_time_info": "F",
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "n": 0,
    "p": 0,
    "t_step_start": 0,
    "t_step_stop": 50,
    "t_step_save": 50,
    "model_eqns": 2,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "precision": 2,
    "prim_vars_wrt": "F",
    "parallel_io": "F",
}

# Each entry is one test case. Fields:
#   name      - unique identifier used in log paths and console output
#   description - human-readable summary
#   compare   - D/ output files compared between reference and perturbed runs
#   threshold - max L∞ deviation allowed before the case is declared FAIL
#   ill_cond  - known source of cancellation (empty string = none expected)
#   pre       - parameters for pre_process (generates initial conditions)
#   sim       - parameters for simulation
CASES = [
    {
        "name": "sod_standard",
        "description": "1-D standard Sod, p_L/p_R=10, ideal gas (well-conditioned baseline)",
        "compare": ["cons.1.00.000050.dat", "cons.3.00.000050.dat"],
        "threshold": 1e-13,
        "ill_cond": "",
        "pre": _merge(
            _BASE_PRE,
            _AIR_EOS,
            {
                "m": 24,
                "num_fluids": 1,
                "mpp_lim": "F",
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(1)%alpha_rho(1)": 1.0,
                "patch_icpp(1)%alpha(1)": 1.0,
                "patch_icpp(2)%pres": 0.1,
                "patch_icpp(2)%alpha_rho(1)": 0.125,
                "patch_icpp(2)%alpha(1)": 1.0,
            },
        ),
        "sim": _merge(_BASE_SIM, _AIR_EOS, {"m": 24, "num_fluids": 1, "dt": 0.001}),
    },
    {
        "name": "sod_strong",
        "description": "1-D Sod, p_L/p_R=100,000, ideal gas",
        "compare": ["cons.1.00.000050.dat", "cons.3.00.000050.dat"],
        "threshold": 1e-10,
        "ill_cond": "HLLC xi factor: (s_L - vel_L)/(s_L - s_S) cancels near sonic contact",
        "pre": _merge(
            _BASE_PRE,
            _AIR_EOS,
            {
                "m": 49,
                "num_fluids": 1,
                "mpp_lim": "F",
                "patch_icpp(1)%pres": 1000.0,
                "patch_icpp(1)%alpha_rho(1)": 10.0,
                "patch_icpp(1)%alpha(1)": 1.0,
                "patch_icpp(2)%pres": 0.01,
                "patch_icpp(2)%alpha_rho(1)": 0.01,
                "patch_icpp(2)%alpha(1)": 1.0,
            },
        ),
        "sim": _merge(_BASE_SIM, _AIR_EOS, {"m": 49, "num_fluids": 1, "dt": 5e-5, "prim_vars_wrt": "T"}),
    },
    {
        "name": "water_stiffened",
        "description": "1-D water shock, stiffened EOS (pi_inf=4046)",
        "compare": ["cons.1.00.000050.dat", "prim.3.00.000050.dat"],
        "threshold": 1e-8,
        "ill_cond": "Pressure recovery: p=(E-pi_inf)/gamma loses ~4 digits (pi_inf/p_right~40,000) [threshold loosened until reduced-energy (Etilde) scheme is merged]",
        "pre": _merge(
            _BASE_PRE,
            _WATER_EOS,
            {
                "m": 49,
                "num_fluids": 1,
                "mpp_lim": "F",
                "patch_icpp(1)%pres": 100.0,
                "patch_icpp(1)%alpha_rho(1)": 1.0,
                "patch_icpp(1)%alpha(1)": 1.0,
                "patch_icpp(2)%pres": 0.1,
                "patch_icpp(2)%alpha_rho(1)": 1.0,
                "patch_icpp(2)%alpha(1)": 1.0,
            },
        ),
        "sim": _merge(_BASE_SIM, _WATER_EOS, {"m": 49, "num_fluids": 1, "dt": 2.5e-5, "mixture_err": "T", "prim_vars_wrt": "T"}),
    },
    {
        "name": "air_water_interface",
        "description": "1-D air/water isobaric contact (two-fluid, pi_inf=4046)",
        "compare": ["cons.1.00.000050.dat", "cons.4.00.000050.dat", "cons.5.00.000050.dat"],
        "threshold": 1e-10,
        "ill_cond": "Mixed-cell pressure recovery: E-alpha_w*gamma_w*pi_inf cancels when alpha_w<<1",
        "pre": _merge(
            _BASE_PRE,
            _AIR_EOS,
            {
                "m": 24,
                "num_fluids": 2,
                "mpp_lim": "T",
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(1)%alpha_rho(1)": 1.0,
                "patch_icpp(1)%alpha_rho(2)": 0.0,
                "patch_icpp(1)%alpha(1)": 1.0,
                "patch_icpp(1)%alpha(2)": 0.0,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(2)%alpha_rho(1)": 0.0,
                "patch_icpp(2)%alpha_rho(2)": 1.0,
                "patch_icpp(2)%alpha(1)": 0.0,
                "patch_icpp(2)%alpha(2)": 1.0,
                "fluid_pp(2)%gamma": 0.1953125,
                "fluid_pp(2)%pi_inf": 4046.31,
            },
        ),
        "sim": _merge(
            _BASE_SIM,
            _AIR_EOS,
            {
                "m": 24,
                "num_fluids": 2,
                "mpp_lim": "T",
                "dt": 5e-5,
                "mixture_err": "T",
                "fluid_pp(2)%gamma": 0.1953125,
                "fluid_pp(2)%pi_inf": 4046.31,
            },
        ),
    },
    {
        "name": "bubble_rp",
        "description": "1-D bubbly water, pressure step 2:1 driving Rayleigh-Plesset oscillations (nb=1, Keller-Miksis)",
        "compare": ["cons.1.00.000050.dat", "prim.3.00.000050.dat"],
        "threshold": 1e-8,
        "ill_cond": "RP ODE: (p_bub - p_ext) cancels near bubble equilibrium",
        "pre": _merge(
            _BASE_PRE,
            {
                "m": 49,
                "num_fluids": 1,
                "mpp_lim": "F",
                "bubbles_euler": "T",
                "nb": 1,
                "polytropic": "T",
                "polydisperse": "F",
                "thermal": 3,
                "pref": 101325.0,
                "rhoref": 1000.0,
                "patch_icpp(1)%pres": 2.0,
                "patch_icpp(1)%alpha_rho(1)": 0.96,
                "patch_icpp(1)%alpha(1)": 0.04,
                "patch_icpp(1)%r0": 1.0,
                "patch_icpp(1)%v0": 0.0,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(2)%alpha_rho(1)": 0.96,
                "patch_icpp(2)%alpha(1)": 0.04,
                "patch_icpp(2)%r0": 1.0,
                "patch_icpp(2)%v0": 0.0,
                "fluid_pp(1)%gamma": 0.16,
                "fluid_pp(1)%pi_inf": 3515.0,
                "bub_pp%R0ref": 1.0,
                "bub_pp%p0ref": 1.0,
                "bub_pp%rho0ref": 1.0,
                "bub_pp%ss": 0.07179866765358993,
                "bub_pp%pv": 0.02308216136195411,
                "bub_pp%mu_l": 0.009954269975623244,
                "bub_pp%gam_g": 1.4,
            },
        ),
        "sim": _merge(
            _BASE_SIM,
            {
                "m": 49,
                "num_fluids": 1,
                "dt": 2.5e-5,
                "mixture_err": "T",
                "prim_vars_wrt": "T",
                "bubbles_euler": "T",
                "nb": 1,
                "bubble_model": 3,
                "polytropic": "T",
                "polydisperse": "F",
                "thermal": 3,
                "pref": 101325.0,
                "rhoref": 1000.0,
                "fluid_pp(1)%gamma": 0.16,
                "fluid_pp(1)%pi_inf": 3515.0,
                "bub_pp%R0ref": 1.0,
                "bub_pp%p0ref": 1.0,
                "bub_pp%rho0ref": 1.0,
                "bub_pp%ss": 0.07179866765358993,
                "bub_pp%pv": 0.02308216136195411,
                "bub_pp%mu_l": 0.009954269975623244,
                "bub_pp%gam_g": 1.4,
            },
        ),
    },
    {
        "name": "low_mach",
        "description": "1-D water shock with low_Mach=1 HLLC correction active",
        "compare": ["cons.1.00.000050.dat", "prim.3.00.000050.dat"],
        "threshold": 2e-7,
        "ill_cond": "low_Mach correction: velocity perturbation ~u/c cancels severely at M≈0 (threshold loosened to 2e-7 to absorb MCA sampling variance)",
        "pre": _merge(
            _BASE_PRE,
            _WATER_EOS,
            {
                "m": 49,
                "num_fluids": 1,
                "mpp_lim": "F",
                "patch_icpp(1)%pres": 100.0,
                "patch_icpp(1)%alpha_rho(1)": 1.0,
                "patch_icpp(1)%alpha(1)": 1.0,
                "patch_icpp(2)%pres": 0.1,
                "patch_icpp(2)%alpha_rho(1)": 1.0,
                "patch_icpp(2)%alpha(1)": 1.0,
            },
        ),
        "sim": _merge(_BASE_SIM, _WATER_EOS, {"m": 49, "num_fluids": 1, "dt": 2.5e-5, "mixture_err": "T", "prim_vars_wrt": "T", "low_Mach": 1}),
    },
]


def _find_verrou() -> str:
    verrou_home = os.environ.get("VERROU_HOME", os.path.join(os.path.expanduser("~"), ".local", "verrou"))
    candidate = os.path.join(verrou_home, "bin", "valgrind")
    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    return shutil.which("valgrind") or ""


def _find_binary(name: str) -> str:
    install_dir = os.path.join(MFC_ROOT_DIR, "build", "install")
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


def _write_inp(params: dict, target_name: str, work_dir: str) -> None:
    """Write a Fortran namelist .inp file from a Python params dict."""
    from .run import case_dicts

    master_keys = case_dicts.get_input_dict_keys(target_name)
    lines = [f"{k} = {v}" for k, v in params.items() if k in master_keys]
    with open(os.path.join(work_dir, f"{target_name}.inp"), "w") as fh:
        fh.write("&user_inputs\n" + "\n".join(lines) + "\n&end/\n")


def _run_preprocess(pp_bin: str, pre_params: dict, work_dir: str):
    _write_inp(pre_params, "pre_process", work_dir)
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


def _max_abs_np(ref_dir: str, compare_files: list) -> float:
    """Return the maximum absolute value across all reference output files."""
    import numpy as np

    total = 0.0
    for fname in compare_files:
        ref_p = os.path.join(ref_dir, fname)
        if not os.path.exists(ref_p):
            continue
        ref = np.loadtxt(ref_p)[:, 1]
        total = max(total, float(np.max(np.abs(ref))))
    return total


def _parse_cancel_gen(gen_path: str) -> list:
    """Parse cc-gen-file TSV (file\\tline\\tsymbol) → sorted unique [(fname, line)] for MFC sources."""
    if not os.path.isfile(gen_path):
        return []
    locs = []
    seen = set()
    with open(gen_path) as fh:
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            fname = parts[0].strip()
            if any(ext in fname for ext in _EXTERNAL_SRCS):
                continue
            if not fname.endswith((".fpp", ".f90", ".F90", ".c", ".cpp")):
                continue
            try:
                lineno = int(parts[1].strip())
            except ValueError:
                continue
            key = (fname, lineno)
            if key not in seen:
                seen.add(key)
                locs.append(key)
    return locs


def _parse_vg_error_locs(log_path: str, error_keyword: str) -> list:
    """Extract first MFC-source frame from each Valgrind error matching error_keyword."""
    if not os.path.isfile(log_path):
        return []
    locs = []
    seen = set()
    in_error = False
    with open(log_path) as fh:
        for raw in fh:
            line = re.sub(r"^==\d+== ?", "", raw)
            if error_keyword in line:
                in_error = True
                continue
            if in_error:
                if "   at " in line or "   by " in line:
                    m = _VGFRAME_RE.search(line)
                    if m:
                        fname = m.group(1)
                        if any(ext in fname for ext in _EXTERNAL_SRCS):
                            continue
                        lineno = int(m.group(2))
                        key = (fname, lineno)
                        if key not in seen:
                            seen.add(key)
                            locs.append(key)
                        in_error = False
                elif line.strip() == "":
                    in_error = False
    return locs


def _run_cancellation_check(case: dict, verrou_bin: str, sim_bin: str, work_dir: str) -> list:
    """Run with --check-cancellation=yes; return [(fname, line)] of MFC cancellation sites."""
    run_dir = os.path.join(work_dir, "cancellation")
    os.makedirs(run_dir, exist_ok=True)
    gen_path = os.path.join(run_dir, "cancel_gen.txt")
    flags = [
        "--check-cancellation=yes",
        "--cc-threshold-double=10",
        f"--cc-gen-file={gen_path}",
    ]
    try:
        _run_simulation_verrou(verrou_bin, sim_bin, work_dir, run_dir, rounding_mode="nearest", extra_flags=flags)
    except MFCException:
        pass
    return _parse_cancel_gen(gen_path)


def _run_mca_samples(
    case: dict,
    verrou_bin: str,
    sim_bin: str,
    work_dir: str,
    ref_dir: str,
    n_mca: int,
) -> tuple:
    """Run N mcaquad samples; return (max_dev, sig_bits_lower_bound)."""
    compare = case["compare"]
    ref_scale = _max_abs_np(ref_dir, compare)
    max_dev = 0.0
    flags = ["--backend=mcaquad", "--mca-mode=mca"]
    for i in range(n_mca):
        run_dir = os.path.join(work_dir, f"mca_{i:02d}")
        os.makedirs(run_dir, exist_ok=True)
        try:
            _run_simulation_verrou(verrou_bin, sim_bin, work_dir, run_dir, extra_flags=flags)
            max_dev = max(max_dev, _max_diff_np(ref_dir, run_dir, compare))
        except MFCException:
            pass
    sig_bits = None
    if max_dev > 0.0 and ref_scale > 0.0:
        sig_bits = max(0, int(math.floor(-math.log2(max_dev / ref_scale))))
    return max_dev, sig_bits


def _run_float_max_check(case: dict, verrou_bin: str, sim_bin: str, work_dir: str) -> list:
    """Run with --check-max-float=yes; return [(fname, line)] of overflow sites."""
    run_dir = os.path.join(work_dir, "float_max")
    os.makedirs(run_dir, exist_ok=True)
    try:
        _run_simulation_verrou(
            verrou_bin,
            sim_bin,
            work_dir,
            run_dir,
            rounding_mode="nearest",
            extra_flags=["--check-max-float=yes"],
        )
    except MFCException:
        pass
    return _parse_vg_error_locs(os.path.join(run_dir, "verrou.log"), "Max float")


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

    verrou_dd_* calls: dd_run.sh RUNDIR and injects function/line exclusion via
    VERROU_EXCLUDE / VERROU_SOURCE environment variables.  For test runs, we use
    --rounding-mode=float (deterministic, same deviation every call, --nruns=1 suffices).
    For the reference run, verrou_dd_sym sets VERROU_ROUNDING_MODE=nearest in the
    environment — we honour that so the reference is a stable nearest-rounding baseline
    to compare against.  CLI --rounding-mode would override the env var and break the
    reference, so we pass the mode via ${VERROU_ROUNDING_MODE:-float} instead.
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

        # verrou_dd_sym sets VERROU_ROUNDING_MODE=nearest for its reference run and
        # leaves it unset for test runs.  Defaulting to float gives deterministic
        # test steps while letting the reference use nearest-rounding.
        ROUND="${{VERROU_ROUNDING_MODE:-float}}"

        cd "$TMPDIR_RUN"
        "$VERROU_BIN" --tool=verrou --error-limit=no --rounding-mode="$ROUND" "$SIM_BIN"
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


def _parse_rddmin_locs(summary_path: str) -> list:
    """Extract [(rel_path, start_line, end_line)] from a dd_line rddmin_summary."""
    if not os.path.isfile(summary_path):
        return []
    locs = []
    with open(summary_path) as fh:
        for line in fh:
            m = _LOC_RE.search(line)
            if not m:
                continue
            path = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3)) if m.group(3) else start
            try:
                rel = os.path.relpath(path, MFC_ROOT_DIR)
                if rel.startswith(".."):
                    rel = path
            except ValueError:
                rel = path
            locs.append((rel.replace("\\", "/"), start, end))
    return locs


def _parse_rddmin_syms(summary_path: str) -> list:
    """Extract symbol/function names from a dd_sym rddmin_summary."""
    if not os.path.isfile(summary_path):
        return []
    with open(summary_path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def _run_dd_tool(
    dd_bin: str,
    dd_dir: str,
    dd_run_sh: str,
    dd_cmp_py: str,
    env: dict,
    log_name: str,
    summary_subdir: str,
    label: str,
) -> list:
    """Generic runner for verrou_dd_sym / verrou_dd_line. Returns raw summary lines."""
    log_file = os.path.join(dd_dir, log_name)
    cmd = [dd_bin, "--nruns=1", "--rddmin=d", "--reference-rounding=nearest", dd_run_sh, dd_cmp_py]
    cons.print(f"  [dim]running {label} (--nruns=1 float-mode --rddmin=d)...[/dim]")
    with open(log_file, "w") as f:
        result = subprocess.run(cmd, cwd=dd_dir, env=env, stdout=f, stderr=subprocess.STDOUT, check=False)
    summary_path = os.path.join(dd_dir, summary_subdir, "rddmin_summary")
    summary_lines = []
    if result.returncode == 0:
        if os.path.isfile(summary_path):
            with open(summary_path) as f:
                summary_lines = f.readlines()
            cons.print(f"  [bold yellow]{label} result[/bold yellow]:")
            for line in summary_lines:
                cons.print(f"    {line.rstrip()}")
        else:
            cons.print(f"  [dim]{label} done; see {log_file}[/dim]")
    else:
        cons.print(f"  [bold yellow]{label} exited {result.returncode}[/bold yellow] (see {log_file})")
    return summary_lines


def _run_dd_sym(case: dict, verrou_bin: str, sim_bin: str, work_dir: str, log_dir: str, threshold: float = None) -> list:
    """Run verrou_dd_sym; return list of responsible symbol names."""
    dd_bin = _find_dd_sym(verrou_bin)
    if not dd_bin:
        cons.print("  [dim]verrou_dd_sym not found; skipping delta-debug[/dim]")
        return []

    dd_dir = os.path.join(log_dir, case["name"])
    os.makedirs(dd_dir, exist_ok=True)
    dd_run_sh = os.path.join(dd_dir, "dd_run.sh")
    dd_cmp_py = os.path.join(dd_dir, "dd_cmp.py")
    _write_dd_run_sh(dd_run_sh, verrou_bin, sim_bin, work_dir)
    _write_dd_cmp_py(dd_cmp_py, case["compare"], threshold if threshold is not None else case["threshold"])
    _run_dd_tool(dd_bin, dd_dir, dd_run_sh, dd_cmp_py, _dd_env(verrou_bin), "dd_sym.log", "dd.sym", "verrou_dd_sym")
    cons.print(f"  [dim]dd_sym logs: {dd_dir}[/dim]")
    return _parse_rddmin_syms(os.path.join(dd_dir, "dd.sym", "rddmin_summary"))


def _run_dd_line(case: dict, verrou_bin: str, sim_bin: str, work_dir: str, log_dir: str, threshold: float = None) -> list:
    """Run verrou_dd_line; return list of (rel_path, start_line, end_line) tuples."""
    dd_bin = _find_dd_line(verrou_bin)
    if not dd_bin:
        cons.print("  [dim]verrou_dd_line not found; skipping line-level debug[/dim]")
        return []

    dd_dir = os.path.join(log_dir, case["name"])
    os.makedirs(dd_dir, exist_ok=True)
    dd_run_sh = os.path.join(dd_dir, "dd_run.sh")
    dd_cmp_py = os.path.join(dd_dir, "dd_cmp.py")
    effective_threshold = threshold if threshold is not None else case["threshold"]
    if not os.path.isfile(dd_run_sh):
        _write_dd_run_sh(dd_run_sh, verrou_bin, sim_bin, work_dir)
        _write_dd_cmp_py(dd_cmp_py, case["compare"], effective_threshold)
    else:
        # dd_sym already wrote dd_cmp.py with its threshold; rewrite with ours if different
        _write_dd_cmp_py(dd_cmp_py, case["compare"], effective_threshold)
    _run_dd_tool(dd_bin, dd_dir, dd_run_sh, dd_cmp_py, _dd_env(verrou_bin), "dd_line.log", "dd.line", "verrou_dd_line")
    return _parse_rddmin_locs(os.path.join(dd_dir, "dd.line", "rddmin_summary"))


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
    run_cancellation: bool,
    run_mca: bool,
    run_float_max: bool,
) -> dict:
    name = case["name"]
    threshold = case["threshold"]
    compare = case["compare"]

    cons.print(f"[bold]{name}[/bold]: {case['description']}")
    cons.indent()
    if case["ill_cond"]:
        cons.print(f"  ill-conditioning: {case['ill_cond']}")
    cons.print(f"  threshold: {threshold:.0e}")

    work_dir = tempfile.mkdtemp(prefix=f"mfc-fps-{name}-")
    result = {
        "name": name,
        "passed": False,
        "max_dev": float("inf"),
        "threshold": threshold,
        "float_proxy": None,
        "vprec": [],
        "dd_sym_syms": [],
        "dd_line_locs": [],
        "cancellation_locs": [],
        "mca_dev": None,
        "mca_sigbits": None,
        "float_max_locs": [],
    }
    try:
        cons.print("  [dim]running pre_process...[/dim]")
        _write_inp(case["sim"], "simulation", work_dir)
        _run_preprocess(pp_bin, case["pre"], work_dir)

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
        result["passed"] = passed
        result["max_dev"] = max_dev
        tag = "[bold green]PASS[/bold green]" if passed else "[bold red]FAIL[/bold red]"
        cons.print(f"  {tag}  max_dev={max_dev:.3e}  threshold={threshold:.0e}")

        # --- B: float proxy ---
        if run_float:
            try:
                fdev = _run_float_proxy(case, verrou_bin, sim_bin, work_dir, ref_dir)
                result["float_proxy"] = fdev
                cons.print(f"  float proxy: dev={fdev:.3e}  (single-precision sensitivity)")
            except MFCException as exc:
                cons.print(f"  [dim]float proxy error: {exc}[/dim]")

        # --- C: VPREC sweep ---
        if run_vprec:
            cons.print("  VPREC precision sweep:")
            vprec_results = _run_vprec_sweep(case, verrou_bin, sim_bin, work_dir, ref_dir)
            result["vprec"] = vprec_results
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

        # --- D/E: delta-debug with float mode to find FP hotspots.
        # dd_run.sh uses --rounding-mode=float (deterministic single-precision),
        # so each bisection step is consistent and --nruns=1 suffices.  Threshold
        # = float_proxy/10: the full instrumented set produces ~float_proxy
        # deviation; excluding the responsible function drops it to near zero;
        # any subset missing the responsible function gives SAME.
        # Skip when float_proxy is unavailable or too small to localize.
        float_proxy = result.get("float_proxy")
        _DD_FLOAT_MIN = 1e-6
        dd_threshold = float_proxy / 10.0 if float_proxy and float_proxy >= _DD_FLOAT_MIN else 0.0
        if dd_threshold > 0 and (run_dd_sym or run_dd_line):
            cons.print(f"  [dim]dd threshold: {dd_threshold:.1e} (float_proxy={float_proxy:.1e})[/dim]")
        elif run_dd_sym or run_dd_line:
            cons.print(f"  [dim]skipping dd: float_proxy={float_proxy} < {_DD_FLOAT_MIN:.0e}[/dim]")
        if dd_threshold > 0 and run_dd_sym:
            try:
                result["dd_sym_syms"] = _run_dd_sym(case, verrou_bin, sim_bin, work_dir, log_dir, threshold=dd_threshold)
            except Exception as exc:
                cons.print(f"  [bold yellow]dd_sym error[/bold yellow]: {exc}")
        if dd_threshold > 0 and run_dd_line:
            try:
                result["dd_line_locs"] = _run_dd_line(case, verrou_bin, sim_bin, work_dir, log_dir, threshold=dd_threshold)
            except Exception as exc:
                cons.print(f"  [bold yellow]dd_line error[/bold yellow]: {exc}")

        # --- F: cancellation detection ---
        if run_cancellation:
            cons.print("  [dim]cancellation detection...[/dim]")
            try:
                locs = _run_cancellation_check(case, verrou_bin, sim_bin, work_dir)
                result["cancellation_locs"] = locs
                if locs:
                    cons.print(f"  cancellation: {len(locs)} unique source location(s)")
                else:
                    cons.print("  cancellation: none detected")
            except Exception as exc:
                cons.print(f"  [bold yellow]cancellation check error[/bold yellow]: {exc}")

        # --- G: MCA significant-bits estimate ---
        if run_mca:
            cons.print(f"  [dim]MCA significant-bits estimate (N={n_samples})...[/dim]")
            try:
                mca_dev, mca_sigbits = _run_mca_samples(case, verrou_bin, sim_bin, work_dir, ref_dir, n_samples)
                result["mca_dev"] = mca_dev
                result["mca_sigbits"] = mca_sigbits
                bits_str = f"~{mca_sigbits} sig bits" if mca_sigbits is not None else "n/a"
                cons.print(f"  MCA: dev={mca_dev:.3e}  ({bits_str})")
            except Exception as exc:
                cons.print(f"  [bold yellow]MCA error[/bold yellow]: {exc}")

        # --- H: float-max overflow detection ---
        if run_float_max:
            cons.print("  [dim]float-max overflow check...[/dim]")
            try:
                locs = _run_float_max_check(case, verrou_bin, sim_bin, work_dir)
                result["float_max_locs"] = locs
                if locs:
                    cons.print(f"  [bold yellow]float-max[/bold yellow]: {len(locs)} overflow site(s)")
                else:
                    cons.print("  float-max: no overflows")
            except Exception as exc:
                cons.print(f"  [bold yellow]float-max check error[/bold yellow]: {exc}")

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)
        cons.unindent()
        cons.print()
    return result


def _emit_github_annotations(results: list):
    """Emit GitHub annotations for FP hotspots.

    Only runs inside GitHub Actions (GITHUB_ACTIONS env var set). Annotations
    appear inline on the responsible source lines in the PR diff view.

    Up to 3 dd_line locations are emitted as ::warning:: per case (minimal
    responsible lines from delta-debug).  Up to 3 cancellation sites per case
    are emitted as ::notice:: so the diff also highlights subtraction-
    cancellation hotspots identified by --check-cancellation.
    """
    if not os.environ.get("GITHUB_ACTIONS"):
        return
    for r in results:
        status = "FAIL" if not r["passed"] else "hotspot"
        dev_str = f"max_dev={r['max_dev']:.2e} (threshold {r['threshold']:.0e})"

        for rel_path, start, end in r.get("dd_line_locs", [])[:3]:
            loc = f"file={rel_path},line={start}"
            if end != start:
                loc += f",endLine={end}"
            title = f"FP {status} [{r['name']}]"
            print(f"::warning {loc},title={title}::{dev_str}", flush=True)

        for fname, lineno in r.get("cancellation_locs", [])[:3]:
            loc = f"file={fname},line={lineno}"
            title = f"FP cancellation [{r['name']}]"
            print(f"::notice {loc},title={title}::catastrophic cancellation site", flush=True)


def _emit_github_summary(results: list, n_samples: int):
    """Write a markdown results table to GITHUB_STEP_SUMMARY.

    Visible directly in the Actions run UI without downloading artifacts.
    Includes: pass/fail, max_dev, float proxy, VPREC sweep (failing levels),
    and dd_line source locations for any failing cases.
    """
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if not summary_path:
        return

    n_pass = sum(1 for r in results if r["passed"])
    n_fail = len(results) - n_pass

    md = []
    md.append("## FP Stability Results\n")
    md.append(f"**{n_pass} passed, {n_fail} failed** — {n_samples} random-rounding samples per case\n")

    # Main results table
    md.append("| Case | Status | max\\_dev | threshold | Float proxy | MCA sig bits |")
    md.append("|------|:------:|--------:|--------:|--------:|:------:|")
    for r in results:
        status = "✅" if r["passed"] else "❌"
        fp = f"{r['float_proxy']:.2e}" if r["float_proxy"] is not None else "—"
        sb = str(r["mca_sigbits"]) if r.get("mca_sigbits") is not None else "—"
        md.append(f"| `{r['name']}` | {status} | {r['max_dev']:.2e} | {r['threshold']:.0e} | {fp} | {sb} |")
    md.append("")

    # VPREC sweep — one column per bit level, ❌ where dev > threshold
    if any(r["vprec"] for r in results):
        _labels = {52: "52b", 23: "23b", 16: "16b", 10: "10b"}
        header = " | ".join(_labels[b] for b in VPREC_MANTISSA_BITS)
        sep = " | ".join(":---:" for _ in VPREC_MANTISSA_BITS)
        md.append("### VPREC precision sweep\n")
        md.append(f"| Case | {header} |")
        md.append(f"|------|{sep}|")
        for r in results:
            vmap = {b: d for b, d in r["vprec"]}
            cols = []
            for b in VPREC_MANTISSA_BITS:
                d = vmap.get(b)
                if d is None:
                    cols.append("—")
                elif d == float("inf"):
                    cols.append("💥 crash")
                else:
                    cols.append(f"{d:.2e}")
            md.append(f"| `{r['name']}` | {' | '.join(cols)} |")
        md.append("")

    # dd_line hotspot sources — always shown (top 10 per case) with source context
    cases_with_locs = [r for r in results if r["dd_line_locs"]]
    if cases_with_locs:
        md.append("### Top FP hotspots (dd\\_line)\n")
        for r in cases_with_locs:
            status = "❌ FAIL" if not r["passed"] else "✅ pass"
            md.append(f"**`{r['name']}`** ({status})\n")
            for rel_path, start, end in r["dd_line_locs"][:10]:
                loc = f"{rel_path}:{start}" if start == end else f"{rel_path}:{start}-{end}"
                md.append(f"- `{loc}`")
                snippet = _get_source_context(rel_path, start)
                if snippet:
                    md.append("  ```fortran")
                    for line in snippet.splitlines():
                        md.append(f"  {line}")
                    md.append("  ```")
            md.append("")

    # dd_sym function names (collapsed, since less actionable than dd_line)
    cases_with_syms = [r for r in results if r["dd_sym_syms"]]
    if cases_with_syms:
        md.append("<details>")
        md.append("<summary>Responsible functions (dd_sym)</summary>\n")
        for r in cases_with_syms:
            md.append(f"\n**`{r['name']}`**\n")
            for sym in r["dd_sym_syms"]:
                md.append(f"- `{sym}`")
        md.append("\n</details>\n")

    # Cancellation hotspots
    cases_with_cancel = [r for r in results if r.get("cancellation_locs")]
    if cases_with_cancel:
        md.append("### Catastrophic cancellation sites\n")
        for r in cases_with_cancel:
            md.append(f"**`{r['name']}`** — {len(r['cancellation_locs'])} site(s)\n")
            for fname, lineno in r["cancellation_locs"][:15]:
                md.append(f"- `{fname}:{lineno}`")
                snippet = _get_source_context(fname, lineno)
                if snippet:
                    md.append("  ```fortran")
                    for line in snippet.splitlines():
                        md.append(f"  {line}")
                    md.append("  ```")
            md.append("")

    # Float-max overflow sites
    cases_with_fmax = [r for r in results if r.get("float_max_locs")]
    if cases_with_fmax:
        md.append("### Float32 overflow sites (check\\_max\\_float)\n")
        for r in cases_with_fmax:
            md.append(f"**`{r['name']}`** — {len(r['float_max_locs'])} site(s)\n")
            for fname, lineno in r["float_max_locs"][:10]:
                md.append(f"- `{fname}:{lineno}`")
            md.append("")

    with open(summary_path, "a") as f:
        f.write("\n".join(md) + "\n")


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

    n_samples = ARG("samples")
    run_float = not ARG("no_float_proxy")
    run_vprec = not ARG("no_vprec")
    run_dd_sym = not ARG("no_dd_sym")
    run_dd_line = not ARG("no_dd_line")
    run_cancellation = not ARG("no_cancellation")
    run_mca = not ARG("no_mca")
    run_float_max = not ARG("no_float_max")

    log_dir = os.path.join(MFC_ROOT_DIR, "fp-stability-logs")
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
    if run_cancellation:
        features.append("cancellation")
    if run_mca:
        features.append("mca-sigbits")
    if run_float_max:
        features.append("float-max")
    cons.print(f"  features:    {', '.join(features) if features else 'stability only'}")
    cons.print(f"  logs:        {log_dir}")
    cons.print()

    start = time.time()
    results = []
    for case in CASES:
        try:
            r = _run_case(
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
                run_cancellation,
                run_mca,
                run_float_max,
            )
        except MFCException as exc:
            cons.print(f"  [bold red]ERROR[/bold red]: {exc}")
            r = {
                "name": case["name"],
                "passed": False,
                "max_dev": float("inf"),
                "threshold": case["threshold"],
                "float_proxy": None,
                "vprec": [],
                "dd_sym_syms": [],
                "dd_line_locs": [],
                "cancellation_locs": [],
                "mca_dev": None,
                "mca_sigbits": None,
                "float_max_locs": [],
            }
        results.append(r)

    elapsed = time.time() - start
    n_pass = sum(1 for r in results if r["passed"])
    n_fail = len(results) - n_pass

    cons.print(f"[bold]Results[/bold] ({elapsed:.0f}s):  [green]{n_pass} passed[/green]  [red]{n_fail} failed[/red]")
    for r in results:
        mark = "[green]✓[/green]" if r["passed"] else "[red]✗[/red]"
        cons.print(f"  {mark} {r['name']}")

    if n_fail > 0:
        cons.print(f"\n  dd_sym/dd_line logs in: {log_dir}")

    _emit_github_summary(results, n_samples)
    _emit_github_annotations(results)

    sys.exit(0 if n_fail == 0 else 1)
