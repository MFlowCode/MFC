#!/usr/bin/env python3
"""
Convergence-rate verification for MFC's 1D single-fluid Euler equations.

Single fluid with a density sine wave: rho = 1 + 0.2*sin(2*pi*x), u=1, p=1.
After exactly one period (T=1, u=1, L=1), the exact solution equals the IC.
L2(rho(T) - rho(0)) measures the accumulated scheme spatial truncation error.
No non-conservative alpha equation — clean benchmark for all schemes.

WENO5/TENO5 use CFL=0.02: RK3 temporal error O(dt^3) is then negligible
relative to the O(h^5) spatial error at N=128-512.

WENO7/TENO7 use CFL=0.005: at CFL=0.02 the RK3 temporal error (~3.4e-12 at
N=128) is comparable to the spatial error (~4.4e-12), giving a spurious rate
of ~3.7.  With CFL=0.005 the temporal error drops by (0.005/0.02)^3 = 1/64
to ~5.3e-14, well below spatial, and the measured rate approaches 7.
N is capped at 256 — the machine-precision floor is reached near N=512.

WENO3-JS degrades to 2nd order at smooth extrema (Henrick et al. 2005).
The expected rate for WENO3 here is therefore 2, not 3; the 2D isentropic
vortex test (run_convergence.py) verifies WENO3 rate 3.

MUSCL2 uses muscl_lim=0 (unlimited central-difference) by default.  TVD
limiters clip slopes to zero at smooth extrema and stall at 1st order on the
sine wave; the unlimited limiter preserves 2nd-order convergence everywhere.

Usage:
    python toolchain/mfc/test/run_convergence_1d.py [--resolutions 128 256 512 1024]
"""

import argparse
import sys
import tempfile

from _convergence_common import (
    CONS_TOL,
    conservation_errors,
    fit_rate,
    l2_norm,
    pairwise_rates,
    print_conservation_check,
    print_summary,
    read_cons_var,
    run_mfc_case,
    run_with_traceback,
)

CASE = "examples/1D_euler_convergence/case.py"

# (label, extra_args, expected_order, tolerance, min_N, max_N)
# Per-scheme resolution bounds let each scheme run over the range where its
# asymptotic order is cleanly visible. WENO7/TENO7 use CFL=0.005 to push the
# RK3 temporal floor below the spatial error; everyone else uses CFL=0.02.
SCHEMES = [
    ("WENO5", ["--order", "5", "--cfl", "0.02"], 5, 0.2, 128, 512),
    ("WENO3", ["--order", "3", "--cfl", "0.02"], 2, 0.2, 256, None),
    ("WENO1", ["--order", "1", "--cfl", "0.02"], 1, 0.05, 128, None),
    ("MUSCL2", ["--muscl", "--cfl", "0.02"], 2, 0.1, 128, None),
    ("TENO5", ["--order", "5", "--teno", "--teno-ct", "1e-6", "--cfl", "0.02"], 5, 0.2, 128, 512),
    ("WENO7", ["--order", "7", "--cfl", "0.005"], 7, 0.5, 64, 128),
    ("TENO7", ["--order", "7", "--teno", "--teno-ct", "1e-9", "--cfl", "0.005"], 7, 0.5, 64, 128),
]

# 1D single-fluid Euler (model_eqns=2, num_fluids=1): vf1=ρ, vf2=ρu, vf3=E
CONS_VARS = [("density", 1), ("x-momentum", 2), ("energy", 3)]


def test_scheme(label, extra_args, expected_order, tol, resolutions, min_N=None, max_N=None, num_ranks=1):
    if min_N is not None:
        resolutions = [N for N in resolutions if N >= min_N]
    if max_N is not None:
        resolutions = [N for N in resolutions if N <= max_N]
    print(f"\n{'=' * 60}\n  {label}  (need rate >= {expected_order - tol:.1f})\n{'=' * 60}")

    errors = []
    nts = []
    all_cons = []
    with tempfile.TemporaryDirectory() as tmpdir:
        for N in resolutions:
            dx = 1.0 / N
            cfg, run_dir = run_mfc_case(CASE, tmpdir, f"N{N}", ["-N", str(N)] + extra_args, num_ranks)
            Nt = int(cfg["t_step_stop"])
            nts.append(Nt)
            vf0 = read_cons_var(run_dir, 0, 1, num_ranks, expected_size=N)
            vfT = read_cons_var(run_dir, Nt, 1, num_ranks, expected_size=N)
            errors.append(l2_norm(vfT - vf0, dx))
            all_cons.append(conservation_errors(run_dir, Nt, dx, CONS_VARS, num_ranks, expected_size=N))

    dxs = [1.0 / N for N in resolutions]
    rates = pairwise_rates(errors, dxs)

    print(f"\n  {'N':>6}  {'Nt':>6}  {'dx':>10}  {'L2 error':>14}  {'rate':>8}")
    print(f"  {'-' * 6}  {'-' * 6}  {'-' * 10}  {'-' * 14}  {'-' * 8}")
    for i, N in enumerate(resolutions):
        r_str = f"{rates[i]:>8.2f}" if rates[i] is not None else f"{'---':>8}"
        print(f"  {N:>6}  {nts[i]:>6}  {dxs[i]:>10.6f}  {errors[i]:>14.6e}  {r_str}")

    if len(resolutions) > 1:
        overall = fit_rate(errors, dxs)
        print(f"\n  Fitted rate: {overall:.2f}  (need >= {expected_order - tol:.1f})")
        rate_passed = overall >= expected_order - tol
    else:
        rate_passed = True

    cons_passed = print_conservation_check(all_cons, CONS_VARS, CONS_TOL)
    passed = rate_passed and cons_passed
    print(f"  {'PASS' if passed else 'FAIL'}")
    return passed


def main():
    parser = argparse.ArgumentParser(description="MFC 1D advection convergence-rate verification")
    parser.add_argument("--resolutions", type=int, nargs="+", default=[64, 128, 256, 512, 1024])
    parser.add_argument("--schemes", nargs="+", default=[s[0] for s in SCHEMES])
    parser.add_argument("--muscl-lim", type=int, default=0, help="MUSCL limiter (0=unlimited 1=minmod ...)")
    parser.add_argument("--num-ranks", type=int, default=1)
    args = parser.parse_args()

    muscl_extra = ["--muscl-lim", str(args.muscl_lim)]
    results = {}
    for label, extra_args, expected_order, tol, min_N, max_N in SCHEMES:
        if label not in args.schemes:
            continue
        results[label] = run_with_traceback(label, test_scheme, label, extra_args + muscl_extra, expected_order, tol, args.resolutions, min_N, max_N, args.num_ranks)

    sys.exit(0 if print_summary(results, label_width=12) else 1)


if __name__ == "__main__":
    main()
