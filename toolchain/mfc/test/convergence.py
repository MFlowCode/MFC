"""Convergence-rate / temporal-order test runners.

Invoked from `_handle_case` in test.py when `case.kind == "convergence"`.
Each scheme is registered as its own TestCase in cases.py with a `spec` dict
that names a `runner` ("1d_advection", "2d_vortex", "temporal", "sod_l1") and
the parameters for that runner.
"""

import io
import json
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile
from contextlib import redirect_stdout

import numpy as np

MFC = ".\\mfc.bat" if os.name == "nt" else "./mfc.sh"
CONS_TOL = 1e-10


def _read_cons_var(run_dir, step, var_idx, num_ranks=1, expected_size=None):
    """Read q_cons_vf{var_idx} from all ranks; rank-order concatenation."""
    chunks = []
    for rank in range(num_ranks):
        path = os.path.join(run_dir, "p_all", f"p{rank}", str(step), f"q_cons_vf{var_idx}.dat")
        with open(path, "rb") as f:
            rec_len = struct.unpack("i", f.read(4))[0]
            data = np.frombuffer(f.read(rec_len), dtype=np.float64)
            f.read(4)
        chunks.append(data.copy())
    combined = np.concatenate(chunks)
    if expected_size is not None and combined.size != expected_size:
        raise ValueError(f"Expected {expected_size} values across {num_ranks} ranks, got {combined.size}")
    return combined


def _conservation_errors(run_dir, Nt, cell_vol, var_list, num_ranks, expected_size=None):
    """|Σq(T) - Σq(0)| / |Σq(0)| for each named variable."""
    errs = {}
    for name, idx in var_list:
        q0 = _read_cons_var(run_dir, 0, idx, num_ranks, expected_size)
        qT = _read_cons_var(run_dir, Nt, idx, num_ranks, expected_size)
        s0 = float(np.sum(q0)) * cell_vol
        sT = float(np.sum(qT)) * cell_vol
        errs[name] = abs(sT - s0) / (abs(s0) + 1e-300)
    return errs


def _l2_norm(diff, scale):
    return float(np.sqrt(np.sum(diff**2) * scale))


def _fit_rate(errors, h_values):
    log_h = np.log(np.array(h_values, dtype=float))
    log_err = np.log(np.array(errors, dtype=float))
    slope, _ = np.polyfit(log_h, log_err, 1)
    return float(slope)


def _pairwise_rates(errors, h_values):
    rates = [None]
    for i in range(1, len(errors)):
        log_h0 = math.log(h_values[i - 1])
        log_h1 = math.log(h_values[i])
        rates.append((math.log(errors[i]) - math.log(errors[i - 1])) / (log_h1 - log_h0))
    return rates


def _run_mfc_case(case_path, tmpdir, run_tag, case_args, num_ranks=1):
    """Run case.py once and copy p_all to tmpdir/run_tag. Returns (cfg_dict, run_dir)."""
    result = subprocess.run(
        [sys.executable, case_path, "--mfc", "{}"] + case_args,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"case.py failed:\n{result.stderr}")
    cfg = json.loads(result.stdout)

    cmd = [MFC, "run", case_path, "-t", "pre_process", "simulation", "-n", str(num_ranks), "--"] + case_args
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.getcwd(), check=False)
    if result.returncode != 0:
        print(result.stdout[-3000:])
        print(result.stderr)
        raise RuntimeError(f"./mfc.sh run failed for {run_tag}")

    case_dir = os.path.dirname(case_path)
    src = os.path.join(case_dir, "p_all")
    dst = os.path.join(tmpdir, run_tag, "p_all")
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    shutil.rmtree(src, ignore_errors=True)
    shutil.rmtree(os.path.join(case_dir, "D"), ignore_errors=True)

    return cfg, os.path.join(tmpdir, run_tag)


def _print_conservation_check(all_cons_errs, var_list, tol=CONS_TOL):
    print(f"\n  Conservation (need rel. error < {tol:.0e}):")
    passed = True
    for name, _ in var_list:
        max_err = max(ce[name] for ce in all_cons_errs)
        ok = max_err < tol
        print(f"    {name:<14}: max = {max_err:.2e}  {'OK' if ok else 'FAIL'}")
        if not ok:
            passed = False
    return passed


# Runners.


def _run_resolution_sweep(spec):
    """1D/2D resolution sweep: L2(q(T) - q(0)) vs h, expect rate >= expected_order - tol."""
    case_path = spec["case_path"]
    extra_args = list(spec.get("extra_args", []))
    expected_order = spec["expected_order"]
    tol = spec["tol"]
    resolutions = list(spec["resolutions"])
    num_ranks = spec.get("num_ranks", 1)
    domain_len = spec.get("domain_len", 1.0)
    ndim = spec.get("ndim", 1)
    cons_vars = spec["cons_vars"]
    primary_idx = spec.get("primary_idx", 1)

    if "min_N" in spec and spec["min_N"] is not None:
        resolutions = [N for N in resolutions if N >= spec["min_N"]]
    if "max_N" in spec and spec["max_N"] is not None:
        resolutions = [N for N in resolutions if N <= spec["max_N"]]

    print(f"  (need rate >= {expected_order - tol:.1f})")

    errors = []
    nts = []
    all_cons = []
    with tempfile.TemporaryDirectory() as tmpdir:
        for N in resolutions:
            dx = domain_len / N
            cell_vol = dx**ndim
            cell_count = N**ndim
            cfg, run_dir = _run_mfc_case(case_path, tmpdir, f"N{N}", ["-N", str(N)] + extra_args, num_ranks)
            Nt = int(cfg["t_step_stop"])
            nts.append(Nt)
            q0 = _read_cons_var(run_dir, 0, primary_idx, num_ranks, expected_size=cell_count)
            qT = _read_cons_var(run_dir, Nt, primary_idx, num_ranks, expected_size=cell_count)
            errors.append(_l2_norm(qT - q0, cell_vol))
            all_cons.append(_conservation_errors(run_dir, Nt, cell_vol, cons_vars, num_ranks, expected_size=cell_count))

    dxs = [domain_len / N for N in resolutions]
    rates = _pairwise_rates(errors, dxs)

    print(f"\n  {'N':>6}  {'Nt':>6}  {'dx':>10}  {'L2 error':>14}  {'rate':>8}")
    print(f"  {'-' * 6}  {'-' * 6}  {'-' * 10}  {'-' * 14}  {'-' * 8}")
    for i, N in enumerate(resolutions):
        r_str = f"{rates[i]:>8.2f}" if rates[i] is not None else f"{'---':>8}"
        print(f"  {N:>6}  {nts[i]:>6}  {dxs[i]:>10.6f}  {errors[i]:>14.6e}  {r_str}")

    if len(resolutions) > 1:
        overall = _fit_rate(errors, dxs)
        print(f"\n  Fitted rate: {overall:.2f}  (need >= {expected_order - tol:.1f})")
        rate_passed = overall >= expected_order - tol
    else:
        rate_passed = True

    cons_passed = _print_conservation_check(all_cons, cons_vars, CONS_TOL)
    return rate_passed and cons_passed


def _run_temporal(spec):
    """Fixed N, vary CFL: rate of L2 error vs dt should match temporal scheme order."""
    case_path = spec["case_path"]
    extra_args = list(spec.get("extra_args", []))
    expected_order = spec["expected_order"]
    tol = spec["tol"]
    cfls = list(spec["cfls"])
    num_ranks = spec.get("num_ranks", 1)
    N = spec["N"]
    cons_vars = spec["cons_vars"]
    primary_idx = spec.get("primary_idx", 1)

    print(f"  N={N}  (need rate >= {expected_order - tol:.1f})")

    dx = 1.0 / N
    errors = []
    dts = []
    nts = []
    all_cons = []
    with tempfile.TemporaryDirectory() as tmpdir:
        for cfl in cfls:
            tag = f"cfl{cfl:.4f}".replace(".", "p")
            args = ["-N", str(N), "--cfl", str(cfl)] + extra_args
            cfg, run_dir = _run_mfc_case(case_path, tmpdir, tag, args, num_ranks)
            Nt = int(cfg["t_step_stop"])
            dts.append(float(cfg["dt"]))
            nts.append(Nt)
            q0 = _read_cons_var(run_dir, 0, primary_idx, num_ranks, expected_size=N)
            qT = _read_cons_var(run_dir, Nt, primary_idx, num_ranks, expected_size=N)
            errors.append(_l2_norm(qT - q0, dx))
            all_cons.append(_conservation_errors(run_dir, Nt, dx, cons_vars, num_ranks, expected_size=N))

    rates = _pairwise_rates(errors, dts)

    print(f"\n  {'CFL':>7}  {'dt':>12}  {'Nt':>6}  {'L2 error':>14}  {'rate':>8}")
    print(f"  {'-' * 7}  {'-' * 12}  {'-' * 6}  {'-' * 14}  {'-' * 8}")
    for i, cfl in enumerate(cfls):
        r_str = f"{rates[i]:>8.2f}" if rates[i] is not None else f"{'---':>8}"
        print(f"  {cfl:>7.3f}  {dts[i]:>12.6e}  {nts[i]:>6}  {errors[i]:>14.6e}  {r_str}")

    if len(cfls) > 1:
        overall = _fit_rate(errors, dts)
        print(f"\n  Fitted rate: {overall:.2f}  (need >= {expected_order - tol:.1f})")
        rate_passed = overall >= expected_order - tol
    else:
        rate_passed = True

    cons_passed = _print_conservation_check(all_cons, cons_vars, CONS_TOL)
    return rate_passed and cons_passed


def _l1_self_error(coarse, fine, dx_coarse):
    """L1 diff between coarse solution and 2:1 cell-averaged fine solution."""
    assert len(fine) == 2 * len(coarse), f"Expected 2:1 ratio, got {len(fine)}:{len(coarse)}"
    fine_avg = (fine[0::2] + fine[1::2]) / 2.0
    return float(np.sum(np.abs(coarse - fine_avg)) * dx_coarse)


def _run_sod_l1(spec):
    """1D L1 self-convergence (consecutive 2x-apart resolutions)."""
    case_path = spec["case_path"]
    extra_args = list(spec.get("extra_args", []))
    expected_order = spec["expected_order"]
    tol = spec["tol"]
    resolutions = list(spec["resolutions"])
    num_ranks = spec.get("num_ranks", 1)

    if "min_N" in spec and spec["min_N"] is not None:
        resolutions = [N for N in resolutions if N >= spec["min_N"]]

    print(f"  (need L1 rate >= {expected_order - tol:.1f})")

    nts = []
    run_dirs = []
    with tempfile.TemporaryDirectory() as tmpdir:
        for N in resolutions:
            cfg, run_dir = _run_mfc_case(case_path, tmpdir, f"N{N}", ["-N", str(N)] + extra_args, num_ranks)
            nts.append(int(cfg["t_step_stop"]))
            run_dirs.append(run_dir)

        errors = []
        error_resolutions = []
        for i in range(len(resolutions) - 1):
            N_c, N_f = resolutions[i], resolutions[i + 1]
            if N_f != 2 * N_c:
                continue
            rho_c = _read_cons_var(run_dirs[i], nts[i], 1, num_ranks, expected_size=N_c)
            rho_f = _read_cons_var(run_dirs[i + 1], nts[i + 1], 1, num_ranks, expected_size=N_f)
            errors.append(_l1_self_error(rho_c, rho_f, 1.0 / N_c))
            error_resolutions.append(N_c)

    dxs = [1.0 / N for N in error_resolutions]
    rates = [None]
    for i in range(1, len(errors)):
        rates.append((math.log(errors[i]) - math.log(errors[i - 1])) / (math.log(dxs[i]) - math.log(dxs[i - 1])))

    print(f"\n  {'N':>6}  {'Nt':>6}  {'L1 self-err':>14}  {'rate':>8}")
    print(f"  {'-' * 6}  {'-' * 6}  {'-' * 14}  {'-' * 8}")
    for i, N in enumerate(error_resolutions):
        r_str = f"{rates[i]:>8.2f}" if rates[i] is not None else f"{'---':>8}"
        print(f"  {N:>6}  {nts[i]:>6}  {errors[i]:>14.6e}  {r_str}")

    if len(errors) >= 2:
        overall = _fit_rate(errors, dxs)
        print(f"\n  Fitted rate: {overall:.2f}  (need >= {expected_order - tol:.1f})")
        return overall >= expected_order - tol
    if len(errors) == 1:
        print(f"\n  Single pair rate: {rates[-1]:.2f}  (need >= {expected_order - tol:.1f})")
        return rates[-1] >= expected_order - tol
    print("\n  ERROR: need >= 2 consecutive 2x-apart resolutions to compute a rate")
    return False


_RUNNERS = {
    "1d_advection": _run_resolution_sweep,
    "2d_advection": _run_resolution_sweep,
    "temporal": _run_temporal,
    "sod_l1": _run_sod_l1,
}


def run_convergence_case(spec):
    """Dispatch on spec['runner']. Returns (passed: bool, captured_output: str)."""
    runner = _RUNNERS.get(spec["runner"])
    if runner is None:
        raise ValueError(f"unknown convergence runner '{spec['runner']}' (have: {sorted(_RUNNERS)})")

    buf = io.StringIO()
    try:
        with redirect_stdout(buf):
            passed = runner(spec)
    except Exception as exc:
        return False, buf.getvalue() + f"\nERROR: {exc}\n"
    return passed, buf.getvalue()
