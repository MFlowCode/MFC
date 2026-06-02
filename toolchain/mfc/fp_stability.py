"""
Floating-point stability test suite using Verrou.

Features
--------
A. Stability suite (always)
   N random-rounding samples per case; PASS/FAIL on significant bits retained
   (scale-free: -log2(max_dev/scale) vs one global floor, no per-case threshold).

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
   Each reported line is then *confirmed* by a positive control: --gen-source
   captures the symbol-correct executed lines, those are filtered to the suspect
   set, and a float-mode run with --source restricted to just them must
   reproduce the instability.  If perturbing the suspect set does not reproduce
   it, the case's hotspots are reported as unconfirmed (downgraded from
   ::warning:: to ::notice::) — this is a single set-level verdict, not per line.
   Each line is then perturbed alone and ranked by the share of the single-
   precision deviation it reproduces.  NOTE: that share is a *sensitivity*
   measure — where reduced precision most moves the output — typically dominated
   by the time integrator / final accumulation, NOT by where cancellation
   originates.  Stage F is the cancellation-origin view; the two usually differ.
   Hotspots are cross-referenced against the stage-F cancellation sites and
   flagged as instance-ambiguous when the .fpp line sits inside a #:for/#:def
   expansion.

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
  ./mfc.sh fp-stability                       # built-in 1-D suite
  ./mfc.sh fp-stability my_case.py            # your own case (small/short, serial, CPU)
  ./mfc.sh fp-stability --no-vprec --no-dd-line
  ./mfc.sh fp-stability --sim-binary PATH --pre-binary PATH

A user case .py is run as a single serial CPU process under Verrou, so it must be
a small, short proxy (a feasibility guard rejects large grids / long runs); output
is forced to serial .dat I/O and the files to diff are auto-detected.
"""

import math
import os
import shutil
import sys
import tempfile
import time

from .common import MFC_ROOT_DIR, MFCException
from .fp_stability_metrics import (
    CANCEL_BIT_LEVELS,
    MIN_SIG_BITS,
    _autodetect_compare,
    _cancellation_severity,
    _mark_cancellation,
    _max_abs_np,
    _max_diff_np,
    _sig_bits,
)
from .fp_stability_report import (
    _emit_github_annotations,
    _emit_github_summary,
)
from .fp_stability_runners import (
    _find_binary,
    _find_verrou,
    _run_cancellation_check,
    _run_confirmation,
    _run_dd_line,
    _run_dd_sym,
    _run_float_max_check,
    _run_float_proxy,
    _run_mca_samples,
    _run_preprocess,
    _run_simulation_verrou,
    _run_vprec_sweep,
    _write_inp,
)
from .printer import cons
from .state import ARG


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
#   ill_cond  - known source of cancellation (empty string = none expected)
# Pass/fail is scale-free (>= MIN_SIG_BITS significant bits retained), so cases
# need no per-case deviation threshold regardless of field magnitude.
#   pre       - parameters for pre_process (generates initial conditions)
#   sim       - parameters for simulation
CASES = [
    {
        "name": "sod_standard",
        "description": "1-D standard Sod, p_L/p_R=10, ideal gas (well-conditioned baseline)",
        "compare": ["cons.1.00.000050.dat", "cons.3.00.000050.dat"],
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
        "ill_cond": "Pressure recovery: p=(E-pi_inf)/gamma loses ~4 digits (pi_inf/p_right~40,000)",
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
        "ill_cond": "low_Mach correction: velocity perturbation ~u/c cancels severely at M≈0",
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


def _blank_result(name: str) -> dict:
    """A result dict with every field at its empty/unmeasured default."""
    return {
        "name": name,
        "passed": False,
        "max_dev": float("inf"),
        "sig_bits": None,
        "float_proxy": None,
        "vprec": [],
        "dd_sym_syms": [],
        "dd_line_locs": [],
        "dd_line_confirmed": None,
        "dd_line_confirm_dev": None,
        "cancellation_locs": [],
        "cancellation_bits": {},
        "mca_dev": None,
        "mca_sigbits": None,
        "float_max_locs": [],
    }


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
    compare = case["compare"]

    cons.print(f"[bold]{name}[/bold]: {case['description']}")
    cons.indent()
    if case["ill_cond"]:
        cons.print(f"  ill-conditioning: {case['ill_cond']}")
    cons.print(f"  pass floor: >= {MIN_SIG_BITS} significant bits retained")

    work_dir = tempfile.mkdtemp(prefix=f"mfc-fps-{name}-")
    result = _blank_result(name)
    try:
        cons.print("  [dim]running pre_process...[/dim]")
        _write_inp(case["sim"], "simulation", work_dir)
        _run_preprocess(pp_bin, case["pre"], work_dir)

        ref_dir = os.path.join(work_dir, "ref")
        os.makedirs(ref_dir)
        cons.print("  [dim]reference run (rounding=nearest)...[/dim]")
        _run_simulation_verrou(verrou_bin, sim_bin, work_dir, ref_dir, rounding_mode="nearest")

        # For a user case with no fixed compare list, diff whatever the reference
        # run actually wrote (conserved vars at the final step).
        if not compare:
            compare = _autodetect_compare(os.listdir(ref_dir))
            case["compare"] = compare
            if not compare:
                raise MFCException("case produced no cons.*/prim.* output to compare (check t_step_save/t_step_stop and parallel_io)")
            cons.print(f"  [dim]comparing: {', '.join(compare)}[/dim]")

        # --- A: random-rounding stability samples ---
        # Pass/fail is scale-free: bits retained = -log2(max_dev / field-scale),
        # vs one global floor (no per-case hand-tuned absolute threshold).
        ref_scale = _max_abs_np(ref_dir, compare)
        max_dev = 0.0
        cons.print(f"  [dim]random-rounding runs (N={n_samples})...[/dim]")
        for i in range(n_samples):
            run_dir = os.path.join(work_dir, f"run_{i:02d}")
            os.makedirs(run_dir)
            _run_simulation_verrou(verrou_bin, sim_bin, work_dir, run_dir, rounding_mode="random")
            max_dev = max(max_dev, _max_diff_np(ref_dir, run_dir, compare))

        sig_bits = _sig_bits(max_dev, ref_scale)
        passed = sig_bits >= MIN_SIG_BITS
        result["passed"] = passed
        result["max_dev"] = max_dev
        result["sig_bits"] = sig_bits
        tag = "[bold green]PASS[/bold green]" if passed else "[bold red]FAIL[/bold red]"
        cons.print(f"  {tag}  {sig_bits:.1f} bits retained (floor {MIN_SIG_BITS})  max_dev={max_dev:.3e}")

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
                elif _sig_bits(dev, ref_scale) < MIN_SIG_BITS:
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
                result["dd_line_locs"] = _run_dd_line(
                    case,
                    verrou_bin,
                    sim_bin,
                    work_dir,
                    log_dir,
                    threshold=dd_threshold,
                )
                macro_n = sum(1 for loc in result["dd_line_locs"] if loc["macro"])
                if macro_n:
                    cons.print(f"  [dim]dd_line: {macro_n} hotspot(s) inside fypp-expanded code (instance-ambiguous)[/dim]")
            except Exception as exc:
                cons.print(f"  [bold yellow]dd_line error[/bold yellow]: {exc}")

        # --- E2: confirm dd_line hotspots and rank each by its individual share ---
        if dd_threshold > 0 and run_dd_line and result["dd_line_locs"]:
            cons.print("  [dim]confirming + ranking dd_line hotspots (per-line perturbation)...[/dim]")
            try:
                confirmed, cdev, ranked = _run_confirmation(case, verrou_bin, sim_bin, work_dir, ref_dir, result["dd_line_locs"], dd_threshold, float_proxy)
                result["dd_line_locs"] = ranked
                result["dd_line_confirmed"] = confirmed
                result["dd_line_confirm_dev"] = cdev
                if confirmed is True:
                    cons.print(f"  [bold green]dd_line confirmed[/bold green]: suspect-only dev={cdev:.3e} >= {dd_threshold:.1e}")
                elif confirmed is False:
                    cons.print(f"  [bold yellow]dd_line UNCONFIRMED[/bold yellow]: suspect-only dev={cdev:.3e} < {dd_threshold:.1e} (attribution suspect)")
                top = ranked[0] if ranked else None
                if top and top.get("share") is not None:
                    cons.print(f"  highest single-precision sensitivity: {top['path']}:{top['start']} ({top['share'] * 100:.0f}% of float-proxy)")
                    cons.print("  [dim](sensitivity = where reduced precision most moves the output, often the time")
                    cons.print("  [dim] integrator; not necessarily where cancellation originates — see cancellation sites)[/dim]")
            except Exception as exc:
                cons.print(f"  [bold yellow]dd_line confirmation error[/bold yellow]: {exc}")

        # --- F: cancellation detection ---
        if run_cancellation:
            cons.print("  [dim]cancellation detection...[/dim]")
            try:
                # sweep bit thresholds to get per-site severity (bits lost); each
                # run returns None if it failed (distinct from [] = ran, found none)
                level_sites = [(level, _run_cancellation_check(verrou_bin, sim_bin, work_dir, threshold=level)) for level in CANCEL_BIT_LEVELS]
                locs = next((s for lvl, s in level_sites if lvl == CANCEL_BIT_LEVELS[0]), None)
                if locs is None:
                    cons.print("  [bold yellow]cancellation: detection run failed (see logs); not reported[/bold yellow]")
                else:
                    bits = _cancellation_severity([(lvl, s) for lvl, s in level_sites if s is not None])
                    result["cancellation_locs"] = locs
                    result["cancellation_bits"] = bits
                    if locs:
                        worst = max(bits.values()) if bits else 0
                        cons.print(f"  cancellation: {len(locs)} site(s), worst loses ≥ {worst / math.log2(10):.0f} of ~16 digits")
                    else:
                        cons.print("  cancellation: none detected")
                    # cross-reference: label dd_line hotspots that sit on a cancellation site
                    if result["dd_line_locs"] and locs:
                        _mark_cancellation(result["dd_line_locs"], locs)
                        n_xref = sum(1 for loc in result["dd_line_locs"] if loc.get("cancellation"))
                        if n_xref:
                            cons.print(f"  {n_xref} hotspot(s) coincide with a catastrophic-cancellation site")
            except Exception as exc:
                cons.print(f"  [bold yellow]cancellation check error[/bold yellow]: {exc}")

        # --- G: MCA significant-bits estimate ---
        if run_mca:
            cons.print(f"  [dim]MCA significant-bits estimate (N={n_samples})...[/dim]")
            try:
                mca_dev, mca_sigbits, n_ok = _run_mca_samples(case, verrou_bin, sim_bin, work_dir, ref_dir, n_samples)
                if n_ok == 0:
                    cons.print(f"  [bold yellow]MCA: no samples completed (0/{n_samples}; see logs)[/bold yellow]")
                else:
                    result["mca_dev"] = mca_dev
                    result["mca_sigbits"] = mca_sigbits
                    bits_str = f"~{mca_sigbits} sig bits" if mca_sigbits is not None else "n/a"
                    cons.print(f"  MCA: dev={mca_dev:.3e}  ({bits_str})  [{n_ok}/{n_samples} samples]")
            except Exception as exc:
                cons.print(f"  [bold yellow]MCA error[/bold yellow]: {exc}")

        # --- H: float-max overflow detection ---
        if run_float_max:
            cons.print("  [dim]float-max overflow check...[/dim]")
            try:
                locs = _run_float_max_check(verrou_bin, sim_bin, work_dir)
                if locs is None:
                    cons.print("  [bold yellow]float-max: run failed (see logs); not reported[/bold yellow]")
                else:
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


# Verrou is ~30x slower and the suite runs the simulation many times, so a user
# case must be a small, short, single-process proxy. Work = cells x time steps;
# both a huge grid and a long run are rejected (built-in cases are ~1k cell-steps).
FP_CASE_MAX_CELLS = 100_000
FP_CASE_MAX_WORK = 200_000  # cells x t_step_stop


def _load_user_case(input_path: str) -> dict:
    """Build a single fp-stability case from a user case .py.

    The case is run as ONE serial CPU process under Verrou (so it must be small
    and short — a coarsened proxy of a production run, not the real thing); a grid
    too large to be feasible errors. The output files to compare are auto-detected
    from the reference run, so 'compare' is left empty here.
    """
    from .run import input as run_input  # lazy import: avoids a circular import

    params = run_input.load(input_path, None, {}, do_print=False).params
    # Force serial .dat I/O: the suite runs the no-MPI binary as one process and
    # diffs serial cons.*/prim.* files (not the parallel SILO/HDF5 path).
    params["parallel_io"] = "F"
    m, n, p = (int(params.get(k, 0) or 0) for k in ("m", "n", "p"))
    cells = (m + 1) * (n + 1) * (p + 1)
    t_stop = int(params.get("t_step_stop", 0) or 0)
    work = cells * max(t_stop, 1)
    if cells > FP_CASE_MAX_CELLS:
        raise MFCException(f"case has {cells:,} cells — too large for Verrou (~30x slowdown, run many times). " f"Use a coarsened proxy (<= {FP_CASE_MAX_CELLS:,} cells).")
    if work > FP_CASE_MAX_WORK:
        raise MFCException(
            f"case is ~{work:,} cell-steps ({cells:,} cells x {t_stop} time steps) — too slow under "
            f"Verrou (~30x, run many times). Reduce m/n/p or t_step_stop (target <= {FP_CASE_MAX_WORK:,} cell-steps)."
        )
    stem = os.path.splitext(os.path.basename(input_path))[0]
    if stem == "case":  # examples/<name>/case.py — the dir name is more telling
        stem = os.path.basename(os.path.dirname(os.path.abspath(input_path))) or stem
    return {
        "name": stem,
        "description": f"user case {input_path} ({cells} cells, run single-rank on CPU)",
        "compare": [],  # auto-detected from the reference run's output
        "ill_cond": "",
        "pre": params,
        "sim": params,
    }


def fp_stability():
    verrou_bin = ARG("verrou_binary") or _find_verrou()
    if not verrou_bin or not os.path.isfile(verrou_bin):
        cons.print("[bold yellow]SKIP[/bold yellow]: Verrou not found (it is a compiled Valgrind tool, not a pip package).")
        cons.print("  Install it (Linux; ~20 min source build) with:  [bold]bash toolchain/bootstrap/verrou.sh[/bold]")
        cons.print("  Or point at an existing build with --verrou-binary PATH or $VERROU_HOME.")
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

    cases_to_run = [_load_user_case(ARG("input"))] if ARG("input") else CASES

    log_dir = os.path.join(MFC_ROOT_DIR, "fp-stability-logs")
    os.makedirs(log_dir, exist_ok=True)

    cons.print()
    cons.print("[bold]MFC Floating-Point Stability Suite[/bold]")
    cons.print(f"  verrou:      {verrou_bin}")
    cons.print(f"  simulation:  {sim_bin}")
    cons.print(f"  pre_process: {pp_bin}")
    if ARG("input"):
        cons.print(f"  case:        {ARG('input')}  (single serial CPU run under Verrou)")
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
    for case in cases_to_run:
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
            r = _blank_result(case["name"])
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
