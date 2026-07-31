"""Convergence-rate tests for ./mfc.sh test.

A "convergence case" runs MFC at several resolutions (or several CFLs), reads
back the conserved variable, fits log(error) vs log(h) for the rate, and
checks that the rate matches the scheme's nominal order. Each case is one
TestCase with kind="convergence" and a ConvergenceSpec attached; test.py
routes it here via run_case().

Four runner flavours, picked by setting spec.runner in cases.py:

  run_h_sweep   vary N (fixed CFL). cell_shift > 0 makes T = K*h and compares
                to np.roll(q(0), +K) — cheap (Nt = O(1) in N) and reports the
                scheme spatial order p directly. cell_shift == 0 runs to one
                full period and compares to q(0).
  run_dt_sweep  fix N, vary CFL. Measures the time-stepper order.
  run_sod_l1    1D L1 self-convergence: compare each N against 2N after
                cell-averaging the fine grid.
  run_amp_sweep fix the grid, sweep an IC jump amplitude A (one step). Fits
                the order of the pressure response in A — an algebraic flux
                regression, not a spatial-convergence test — and asserts it
                two-sided: |fitted - expected_order| <= tol.

Helpers below: _read_field reads a Fortran-unformatted record per rank and
concatenates; _l2_norm and _fit_slope are obvious one-liners; _run_mfc shells
out to ./mfc.sh run and stashes p_all/ in a tempdir.
"""

import dataclasses
import json
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import typing

import numpy as np

from .. import common

CONS_TOL = 1e-10
MFC = ".\\mfc.bat" if os.name == "nt" else "./mfc.sh"


@dataclasses.dataclass
class ConvergenceSpec:
    """One convergence case. `runner` is the function that drives it."""

    runner: typing.Callable
    case_path: str
    expected_order: float  # scheme order p (or temporal q for dt sweeps)
    tol: float  # passes if measured >= expected_order - tol
    extra_args: typing.List[str] = dataclasses.field(default_factory=list)
    cons_vars: typing.List[typing.Tuple[str, int]] = dataclasses.field(default_factory=list)
    primary_idx: int = 1  # which q_cons_vf{idx} to use for the L2/L1 norm
    num_ranks: int = 1
    # Spatial-sweep cases:
    resolutions: typing.List[int] = dataclasses.field(default_factory=list)
    ndim: int = 1
    domain_len: float = 1.0
    cell_shift: int = 0  # 0 -> period mode; >0 -> compare to K-cell-shifted IC
    # Temporal-sweep cases:
    cfls: typing.List[float] = dataclasses.field(default_factory=list)
    fixed_N: int = 0
    # Amplitude-sweep cases:
    amps: typing.List[float] = dataclasses.field(default_factory=list)


# Low-level helpers.


def _read_field(run_dir: str, step: int, var_idx: int, num_ranks: int, expected_size: int) -> np.ndarray:
    """Read q_cons_vf<var_idx>.dat across all ranks, in rank order."""
    chunks = []
    for rank in range(num_ranks):
        path = os.path.join(run_dir, "p_all", f"p{rank}", str(step), f"q_cons_vf{var_idx}.dat")
        with open(path, "rb") as f:
            rec_len = struct.unpack("i", f.read(4))[0]
            chunks.append(np.frombuffer(f.read(rec_len), dtype=np.float64).copy())
            f.read(4)
    arr = np.concatenate(chunks)
    if arr.size != expected_size:
        raise common.MFCException(f"Expected {expected_size} cells, got {arr.size}")
    return arr


def _l2_norm(diff: np.ndarray, cell_vol: float) -> float:
    return float(np.sqrt(np.sum(diff**2) * cell_vol))


def _fit_slope(errors: typing.List[float], xs: typing.List[float]) -> float:
    """Least-squares slope of log(errors) vs log(xs)."""
    return float(np.polyfit(np.log(xs), np.log(errors), 1)[0])


def _pairwise_slope(err1: float, err2: float, x1: float, x2: float) -> float:
    return math.log(err2 / err1) / math.log(x2 / x1)


def _run_mfc(case_path: str, tmpdir: str, run_tag: str, args: typing.List[str], num_ranks: int) -> typing.Tuple[dict, str]:
    """Run case.py through ./mfc.sh and stash p_all/ in tmpdir/run_tag for reading."""
    cfg_run = subprocess.run(
        [sys.executable, case_path, "--mfc", "{}"] + args,
        capture_output=True,
        text=True,
        check=False,
    )
    if cfg_run.returncode != 0:
        raise common.MFCException(f"case.py failed:\n{cfg_run.stderr}")
    cfg = json.loads(cfg_run.stdout)

    sim = subprocess.run(
        [MFC, "run", case_path, "-t", "pre_process", "simulation", "-n", str(num_ranks), "--"] + args,
        capture_output=True,
        text=True,
        check=False,
    )
    if sim.returncode != 0:
        raise common.MFCException(f"./mfc.sh run failed for {run_tag}\n{sim.stdout[-3000:]}\n{sim.stderr}")

    case_dir = os.path.dirname(case_path)
    src = os.path.join(case_dir, "p_all")
    dst = os.path.join(tmpdir, run_tag, "p_all")
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    shutil.rmtree(src, ignore_errors=True)
    shutil.rmtree(os.path.join(case_dir, "D"), ignore_errors=True)
    return cfg, os.path.join(tmpdir, run_tag)


def _conservation_at_step(run_dir: str, Nt: int, cell_vol: float, cons_vars, num_ranks: int, cell_count: int) -> dict:
    """Per-variable |∫q(T) - ∫q(0)| / |∫q(0)|."""
    out = {}
    for name, idx in cons_vars:
        s0 = float(np.sum(_read_field(run_dir, 0, idx, num_ranks, cell_count))) * cell_vol
        sT = float(np.sum(_read_field(run_dir, Nt, idx, num_ranks, cell_count))) * cell_vol
        out[name] = abs(sT - s0) / (abs(s0) + 1e-300)
    return out


def _conservation_lines(history: typing.List[dict], cons_vars) -> typing.Tuple[bool, typing.List[str]]:
    lines = [f"\n  Conservation (need rel. error < {CONS_TOL:.0e}):"]
    passed = True
    for name, _ in cons_vars:
        max_err = max(c[name] for c in history)
        ok = max_err < CONS_TOL
        passed = passed and ok
        lines.append(f"    {name:<14}: max = {max_err:.2e}  {'OK' if ok else 'FAIL'}")
    return passed, lines


def _table_line(values: typing.List[str], widths: typing.List[int]) -> str:
    return "  " + "  ".join(f"{v:>{w}}" for v, w in zip(values, widths))


# Runners.


def run_h_sweep(spec: ConvergenceSpec) -> typing.Tuple[bool, str]:
    """Vary N at fixed CFL.

    cell_shift > 0: T = K*h with v=1 — compare q(T) to np.roll(q(0), +K).
        Cost is O(1) in N; raw rate is p+1, displayed as scheme order p.
    cell_shift == 0: T = full period — compare q(T) to q(0). Rate = p.
    """
    is_shift = spec.cell_shift > 0
    num_ranks = 1 if is_shift else spec.num_ranks  # cell-shift reshape needs single rank
    label = "spatial order" if is_shift else "rate"
    threshold = spec.expected_order - spec.tol

    errors, hs, nts, history = [], [], [], []
    with tempfile.TemporaryDirectory() as tmpdir:
        for N in spec.resolutions:
            h = spec.domain_len / N
            cell_vol = h**spec.ndim
            cell_count = N**spec.ndim
            args = ["-N", str(N)] + spec.extra_args
            if is_shift:
                args += ["--t-end", str(spec.cell_shift * h)]

            cfg, run_dir = _run_mfc(spec.case_path, tmpdir, f"N{N}", args, num_ranks)
            Nt = int(cfg["t_step_stop"])
            q0 = _read_field(run_dir, 0, spec.primary_idx, num_ranks, cell_count)
            qT = _read_field(run_dir, Nt, spec.primary_idx, num_ranks, cell_count)
            ref = q0
            if is_shift:
                # Wave moves at v=+1: q(T)[i] = q(0)[i-K], so np.roll shift=+K.
                shape = (N,) * spec.ndim
                ref = np.roll(
                    q0.reshape(shape, order="F"),
                    shift=spec.cell_shift,
                    axis=tuple(range(spec.ndim)),
                ).flatten(order="F")

            errors.append(_l2_norm(qT - ref, cell_vol))
            hs.append(h)
            nts.append(Nt)
            history.append(_conservation_at_step(run_dir, Nt, cell_vol, spec.cons_vars, num_ranks, cell_count))

    offset = 1 if is_shift else 0
    widths = [6, 6, 10, 14, 13]
    lines = [
        f"  (need {label} >= {threshold:.1f})",
        "",
        _table_line(["N", "Nt", "dx", "L2 error", label.replace(" ", "_")], widths),
        _table_line(["-" * w for w in widths], widths),
    ]
    for i, N in enumerate(spec.resolutions):
        if i == 0:
            r_str = "---"
        else:
            r_str = f"{_pairwise_slope(errors[i - 1], errors[i], hs[i - 1], hs[i]) - offset:.2f}"
        lines.append(_table_line([str(N), str(nts[i]), f"{hs[i]:.6f}", f"{errors[i]:.6e}", r_str], widths))

    rate_passed = True
    if len(errors) > 1:
        fitted = _fit_slope(errors, hs) - offset
        lines.append(f"\n  Fitted {label}: {fitted:.2f}  (need >= {threshold:.1f})")
        rate_passed = fitted >= threshold

    cons_passed, cons_lines = _conservation_lines(history, spec.cons_vars)
    lines += cons_lines
    return rate_passed and cons_passed, "\n".join(lines)


def run_dt_sweep(spec: ConvergenceSpec) -> typing.Tuple[bool, str]:
    """Fix N, vary CFL — measures temporal scheme order."""
    threshold = spec.expected_order - spec.tol
    h = 1.0 / spec.fixed_N

    errors, dts, nts, history = [], [], [], []
    with tempfile.TemporaryDirectory() as tmpdir:
        for cfl in spec.cfls:
            tag = f"cfl{cfl:.4f}".replace(".", "p")
            args = ["-N", str(spec.fixed_N), "--cfl", str(cfl)] + spec.extra_args
            cfg, run_dir = _run_mfc(spec.case_path, tmpdir, tag, args, spec.num_ranks)
            Nt = int(cfg["t_step_stop"])
            q0 = _read_field(run_dir, 0, spec.primary_idx, spec.num_ranks, spec.fixed_N)
            qT = _read_field(run_dir, Nt, spec.primary_idx, spec.num_ranks, spec.fixed_N)
            errors.append(_l2_norm(qT - q0, h))
            dts.append(float(cfg["dt"]))
            nts.append(Nt)
            history.append(_conservation_at_step(run_dir, Nt, h, spec.cons_vars, spec.num_ranks, spec.fixed_N))

    widths = [7, 12, 6, 14, 8]
    lines = [
        f"  N={spec.fixed_N}  (need rate >= {threshold:.1f})",
        "",
        _table_line(["CFL", "dt", "Nt", "L2 error", "rate"], widths),
        _table_line(["-" * w for w in widths], widths),
    ]
    for i, cfl in enumerate(spec.cfls):
        if i == 0:
            r_str = "---"
        else:
            r_str = f"{_pairwise_slope(errors[i - 1], errors[i], dts[i - 1], dts[i]):.2f}"
        lines.append(_table_line([f"{cfl:.3f}", f"{dts[i]:.6e}", str(nts[i]), f"{errors[i]:.6e}", r_str], widths))

    rate_passed = True
    if len(errors) > 1:
        fitted = _fit_slope(errors, dts)
        lines.append(f"\n  Fitted rate: {fitted:.2f}  (need >= {threshold:.1f})")
        rate_passed = fitted >= threshold

    cons_passed, cons_lines = _conservation_lines(history, spec.cons_vars)
    lines += cons_lines
    return rate_passed and cons_passed, "\n".join(lines)


def run_sod_l1(spec: ConvergenceSpec) -> typing.Tuple[bool, str]:
    """1D L1 self-convergence: compare each N against 2N (fine cell-averaged 2:1)."""
    threshold = spec.expected_order - spec.tol

    nts, run_dirs = [], []
    with tempfile.TemporaryDirectory() as tmpdir:
        for N in spec.resolutions:
            cfg, run_dir = _run_mfc(spec.case_path, tmpdir, f"N{N}", ["-N", str(N)] + spec.extra_args, spec.num_ranks)
            nts.append(int(cfg["t_step_stop"]))
            run_dirs.append(run_dir)

        errors, error_resolutions = [], []
        for i in range(len(spec.resolutions) - 1):
            N_c, N_f = spec.resolutions[i], spec.resolutions[i + 1]
            if N_f != 2 * N_c:
                continue
            rho_c = _read_field(run_dirs[i], nts[i], 1, spec.num_ranks, N_c)
            rho_f = _read_field(run_dirs[i + 1], nts[i + 1], 1, spec.num_ranks, N_f)
            avg = 0.5 * (rho_f[0::2] + rho_f[1::2])
            errors.append(float(np.sum(np.abs(rho_c - avg)) * (1.0 / N_c)))
            error_resolutions.append(N_c)

    hs = [1.0 / N for N in error_resolutions]
    widths = [6, 6, 14, 8]
    lines = [
        f"  (need L1 rate >= {threshold:.1f})",
        "",
        _table_line(["N", "Nt", "L1 self-err", "rate"], widths),
        _table_line(["-" * w for w in widths], widths),
    ]
    for i, N in enumerate(error_resolutions):
        if i == 0:
            r_str = "---"
        else:
            r_str = f"{_pairwise_slope(errors[i - 1], errors[i], hs[i - 1], hs[i]):.2f}"
        lines.append(_table_line([str(N), str(nts[i]), f"{errors[i]:.6e}", r_str], widths))

    if len(errors) >= 2:
        fitted = _fit_slope(errors, hs)
        lines.append(f"\n  Fitted rate: {fitted:.2f}  (need >= {threshold:.1f})")
        return fitted >= threshold, "\n".join(lines)
    if len(errors) == 1:
        rate = _pairwise_slope(errors[0], errors[0], hs[0], hs[0])
        lines.append(f"\n  Single pair rate (need >= {threshold:.1f})")
        return rate >= threshold, "\n".join(lines)
    lines.append("\n  ERROR: need >= 2 consecutive 2x-apart resolutions")
    return False, "\n".join(lines)


def _hypo_pressure_2d(run_dir: str, step: int, cell_count: int, cfg: dict) -> np.ndarray:
    """Pressure from conserved fields for the 2-fluid 2D hypoelastic layout.

    Layout: alpha_rho(1:2), rho*u/rho*v(3:4), E(5), alpha(6:7),
    rho*tau_xx/xy/yy(8:10). Inverts the stiffened-gas EOS with the hypoelastic
    strain energy E_e = sum(tau_c^2)/(4G) per volume, shear component doubled
    (mirrors s_compute_pressure in m_variables_conversion.fpp). The caller
    checks p(step 0) against the analytic IC so a layout or EOS drift fails
    loudly instead of silently mismeasuring.
    """
    f = {i: _read_field(run_dir, step, i, 1, cell_count) for i in (1, 2, 3, 4, 5, 8, 9, 10)}
    gamma_coef = float(cfg["fluid_pp(1)%gamma"])  # 1/(gamma - 1); identical fluids
    G = float(cfg["fluid_pp(1)%G"])
    rho = f[1] + f[2]
    dyn_p = 0.5 * (f[3] ** 2 + f[4] ** 2) / rho
    E_e = ((f[8] / rho) ** 2 + 2.0 * (f[9] / rho) ** 2 + (f[10] / rho) ** 2) / (4.0 * G)
    return (f[5] - dyn_p - E_e) / gamma_coef


def run_amp_sweep(spec: ConvergenceSpec) -> typing.Tuple[bool, str]:
    """Fix the grid, sweep the IC tau_xy jump amplitude A; fit response order.

    The tau_xy jump launches shear waves. Measure the baseline-subtracted
    (A = 0 run) pressure response after one step as Linf over the
    interface-adjacent cells. HLLD pairs the same tangential double-star state
    in its momentum and energy fluxes, cancelling the velocity-linear mismatch:
    order 2 in A. HLLC weights them differently: order 1 (the nonzero uniform
    tangential velocity in the case is what exposes it). Asserted two-sided so
    the HLLC control also guards the case's discriminating power.
    """
    errors, tags = [], []
    baseline = None
    with tempfile.TemporaryDirectory() as tmpdir:
        for amp in [0.0] + list(spec.amps):
            tag = f"amp{amp:.4e}".replace(".", "p").replace("-", "m").replace("+", "")
            cfg, run_dir = _run_mfc(spec.case_path, tmpdir, tag, ["--amp", str(amp)] + spec.extra_args, 1)
            Nt = int(cfg["t_step_stop"])
            mx, my = int(cfg["m"]) + 1, int(cfg["n"]) + 1
            p_start = _hypo_pressure_2d(run_dir, 0, mx * my, cfg)
            p_end = _hypo_pressure_2d(run_dir, Nt, mx * my, cfg)
            p_ic = float(cfg["patch_icpp(1)%pres"])
            if np.max(np.abs(p_start - p_ic)) > 1e-8 * abs(p_ic):
                raise common.MFCException(f"{tag}: step-0 pressure differs from IC p = {p_ic} — conserved-layout or EOS drift")
            if not (np.isfinite(p_end).all() and np.min(p_end) > 0.0):
                raise common.MFCException(f"{tag}: nonfinite or nonpositive pressure after {Nt} step(s)")
            resp = (p_end - p_start).reshape((mx, my), order="F")
            if amp == 0.0:
                baseline = resp
                continue
            # Interface-adjacent columns: the contact sits at x = 0.5.
            x_cc = (np.arange(mx) + 0.5) / mx
            mask = np.abs(x_cc - 0.5) <= 1.01 / mx
            errors.append(float(np.max(np.abs((resp - baseline)[mask, :]))))
            tags.append(amp)

    widths = [12, 14, 8]
    lines = [
        f"  (need |fitted order - {spec.expected_order:.1f}| <= {spec.tol:.1f})",
        "",
        _table_line(["amplitude", "interface Linf", "order"], widths),
        _table_line(["-" * w for w in widths], widths),
    ]
    for i, amp in enumerate(tags):
        r_str = "---" if i == 0 else f"{_pairwise_slope(errors[i - 1], errors[i], tags[i - 1], amp):.2f}"
        lines.append(_table_line([f"{amp:.3e}", f"{errors[i]:.6e}", r_str], widths))

    if min(errors) <= 0.0:
        lines.append("\n  ERROR: zero response — case no longer perturbs the solution")
        return False, "\n".join(lines)
    fitted = _fit_slope(errors, tags)
    passed = abs(fitted - spec.expected_order) <= spec.tol
    lines.append(f"\n  Fitted amplitude order: {fitted:.3f}  (need {spec.expected_order:.1f} +/- {spec.tol:.1f})")
    return passed, "\n".join(lines)


# Entry point used by test.py.


def run_case(spec: ConvergenceSpec) -> typing.Tuple[bool, str]:
    """Run one convergence case. Returns (passed, full output report)."""
    try:
        return spec.runner(spec)
    except Exception as exc:
        return False, f"  ERROR: {exc}\n"
