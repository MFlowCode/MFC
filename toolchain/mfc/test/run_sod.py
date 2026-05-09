#!/usr/bin/env python3
"""
L1 self-convergence study for the 1D Sod shock tube across all MFC schemes.

Sod problem: rho_L=1, u_L=0, p_L=1; rho_R=0.125, u_R=0, p_R=0.1; T=0.2.
Contains a shock, contact discontinuity, and rarefaction fan.

By Godunov's theorem, any conservative monotone scheme converges at 1st order
in L1 for problems with shocks.  Higher-order schemes (WENO5, TENO7, ...) also
achieve L1 rate ~1 globally because the shock contributes an O(h) error that
dominates the smooth-region high-order accuracy.

Self-convergence method: run at N and 2N, cell-average the finer solution to
the coarse grid, compute L1(rho_N - avg(rho_{2N})).  No exact solution needed.
Ranks are read in rank order, which equals spatial order for 1D decomposition.

Usage:
    python toolchain/mfc/test/run_sod.py
    python toolchain/mfc/test/run_sod.py --resolutions 64 128 256 512 --schemes WENO5 TENO5
"""

import argparse
import math
import sys
import tempfile

import numpy as np
from _convergence_common import (
    fit_rate,
    print_summary,
    read_cons_var,
    run_mfc_case,
    run_with_traceback,
)

CASE = "examples/1D_sod_convergence/case.py"

# (label, extra_args, expected_order, tolerance, min_N)
# WENO1 contact smears over O(sqrt(h*T)) → fitted L1 rate ~0.6-0.7.
# SUPERBEE is over-compressive near contacts; min_N=128 skips the
# pre-asymptotic point at N=64 (~0.40) for a reliable fit.
SCHEMES = [
    ("WENO1", ["--order", "1"], 1, 0.5, None),
    ("WENO3", ["--order", "3"], 1, 0.3, None),
    ("WENO5", ["--order", "5"], 1, 0.3, None),
    ("WENO7", ["--order", "7"], 1, 0.3, None),
    ("MUSCL-minmod", ["--muscl", "--muscl-lim", "1"], 1, 0.3, None),
    ("MUSCL-MC", ["--muscl", "--muscl-lim", "2"], 1, 0.3, None),
    ("MUSCL-VanLeer", ["--muscl", "--muscl-lim", "4"], 1, 0.3, None),
    ("MUSCL-SUPERBEE", ["--muscl", "--muscl-lim", "5"], 1, 0.5, 128),
    ("TENO5", ["--order", "5", "--teno", "--teno-ct", "1e-6"], 1, 0.3, None),
    ("TENO7", ["--order", "7", "--teno", "--teno-ct", "1e-9"], 1, 0.3, None),
]


def l1_self_error(coarse, fine, dx_coarse):
    """L1 diff between coarse solution and cell-averaged fine solution."""
    assert len(fine) == 2 * len(coarse), f"Expected 2:1 ratio, got {len(fine)}:{len(coarse)}"
    fine_avg = (fine[0::2] + fine[1::2]) / 2.0
    return float(np.sum(np.abs(coarse - fine_avg)) * dx_coarse)


def test_scheme(label, extra_args, expected_order, tol, resolutions, min_N=None, num_ranks=1):
    if min_N is not None:
        resolutions = [N for N in resolutions if N >= min_N]
    print(f"\n{'=' * 60}\n  {label}  (need L1 rate >= {expected_order - tol:.1f})\n{'=' * 60}")

    nts = []
    run_dirs = []
    with tempfile.TemporaryDirectory() as tmpdir:
        for N in resolutions:
            cfg, run_dir = run_mfc_case(CASE, tmpdir, f"N{N}", ["-N", str(N)] + extra_args, num_ranks)
            nts.append(int(cfg["t_step_stop"]))
            run_dirs.append(run_dir)

        # Compute L1 self-errors: compare each N against 2N
        errors = []
        error_resolutions = []
        for i in range(len(resolutions) - 1):
            N_c, N_f = resolutions[i], resolutions[i + 1]
            if N_f != 2 * N_c:
                continue
            rho_c = read_cons_var(run_dirs[i], nts[i], 1, num_ranks, expected_size=N_c)
            rho_f = read_cons_var(run_dirs[i + 1], nts[i + 1], 1, num_ranks, expected_size=N_f)
            errors.append(l1_self_error(rho_c, rho_f, 1.0 / N_c))
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
        overall = fit_rate(errors, dxs)
        print(f"\n  Fitted rate: {overall:.2f}  (need >= {expected_order - tol:.1f})")
        passed = overall >= expected_order - tol
    elif len(errors) == 1:
        print(f"\n  Single pair rate: {rates[-1]:.2f}  (need >= {expected_order - tol:.1f})")
        passed = rates[-1] >= expected_order - tol
    else:
        print("\n  ERROR: need >= 2 consecutive 2x-apart resolutions to compute a rate")
        passed = False

    print(f"  {'PASS' if passed else 'FAIL'}")
    return passed


def main():
    parser = argparse.ArgumentParser(description="MFC Sod shock tube L1 convergence")
    parser.add_argument("--resolutions", type=int, nargs="+", default=[128, 256, 512, 1024])
    parser.add_argument("--schemes", nargs="+", default=[s[0] for s in SCHEMES])
    parser.add_argument("--num-ranks", type=int, default=1)
    args = parser.parse_args()

    results = {}
    for label, extra_args, expected_order, tol, min_N in SCHEMES:
        if label not in args.schemes:
            continue
        results[label] = run_with_traceback(label, test_scheme, label, extra_args, expected_order, tol, args.resolutions, min_N, args.num_ranks)

    sys.exit(0 if print_summary(results, label_width=18) else 1)


if __name__ == "__main__":
    main()
