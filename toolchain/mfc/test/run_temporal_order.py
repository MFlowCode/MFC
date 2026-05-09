#!/usr/bin/env python3
"""
Time-integration order verification for MFC's RK1, RK2, and RK3 time steppers.

Uses the 1D single-fluid Euler advection problem (rho = 1 + 0.2*sin(2*pi*x),
u=1, p=1, L=1, T=1) with a fine spatial grid (N=512, WENO5) so the spatial
error (~4e-12) is negligible compared to the temporal error at the CFLs tested.

L2(rho(T) - rho(0)) measures total accumulated error.  By fixing N and varying
CFL (and hence dt), the spatial contribution is constant and the measured rate
reflects the time integration order.

CFL ranges are chosen to be within each stepper's stability region and keep
temporal errors well above the ~4e-12 spatial floor:
  RK1 (Euler, 1st order): CFL=[0.10, 0.05] — stable limit ~0.1 with WENO5+LF
    (nearly-imaginary eigenvalues constrain Euler more than TVD RK);
    error ~2.5e-4 and ~1.2e-4 (rate ≈ 1.0)
  RK2 (TVD Heun, 2nd order): CFL=[0.50, 0.25];
    error ~1.2e-6 and ~2.9e-7 (rate ≈ 2.0)
  RK3 (TVD Shu-Osher, 3rd order): CFL=[0.50, 0.25];
    error ~8.3e-10 and ~1.1e-10 (rate ≈ 3.0)

Usage:
    python toolchain/mfc/test/run_temporal_order.py
    python toolchain/mfc/test/run_temporal_order.py --schemes RK3/WENO5 --cfls 0.5 0.25 0.125
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
N_SPATIAL = 512  # fixed spatial resolution

# (label, extra_args, expected_order, tolerance, cfls)
# N=512, WENO5: spatial error ~4e-12, well below temporal error at CFL>=0.05.
# RK1 stable to CFL~0.1 with WENO5+LF; RK2/RK3 stable to CFL~1.
SCHEMES = [
    ("RK1/WENO5", ["--order", "5", "--time-stepper", "1"], 1, 0.1, [0.10, 0.05]),
    ("RK2/WENO5", ["--order", "5", "--time-stepper", "2"], 2, 0.2, [0.50, 0.25]),
    ("RK3/WENO5", ["--order", "5", "--time-stepper", "3"], 3, 0.3, [0.50, 0.25]),
]

# 1D single-fluid Euler (model_eqns=2, num_fluids=1): vf1=ρ, vf2=ρu, vf3=E
CONS_VARS = [("density", 1), ("x-momentum", 2), ("energy", 3)]


def test_scheme(label, extra_args, expected_order, tol, cfls, num_ranks=1):
    print(f"\n{'=' * 60}\n  {label}  N={N_SPATIAL}  (need rate >= {expected_order - tol:.1f})\n{'=' * 60}")

    dx = 1.0 / N_SPATIAL
    errors = []
    dts = []
    nts = []
    all_cons = []
    with tempfile.TemporaryDirectory() as tmpdir:
        for cfl in cfls:
            tag = f"cfl{cfl:.4f}".replace(".", "p")
            args = ["-N", str(N_SPATIAL), "--cfl", str(cfl)] + extra_args
            cfg, run_dir = run_mfc_case(CASE, tmpdir, tag, args, num_ranks)
            Nt = int(cfg["t_step_stop"])
            dts.append(float(cfg["dt"]))
            nts.append(Nt)
            vf0 = read_cons_var(run_dir, 0, 1, num_ranks, expected_size=N_SPATIAL)
            vfT = read_cons_var(run_dir, Nt, 1, num_ranks, expected_size=N_SPATIAL)
            errors.append(l2_norm(vfT - vf0, dx))
            all_cons.append(conservation_errors(run_dir, Nt, dx, CONS_VARS, num_ranks, expected_size=N_SPATIAL))

    rates = pairwise_rates(errors, dts)

    print(f"\n  {'CFL':>7}  {'dt':>12}  {'Nt':>6}  {'L2 error':>14}  {'rate':>8}")
    print(f"  {'-' * 7}  {'-' * 12}  {'-' * 6}  {'-' * 14}  {'-' * 8}")
    for i, cfl in enumerate(cfls):
        r_str = f"{rates[i]:>8.2f}" if rates[i] is not None else f"{'---':>8}"
        print(f"  {cfl:>7.3f}  {dts[i]:>12.6e}  {nts[i]:>6}  {errors[i]:>14.6e}  {r_str}")

    if len(cfls) > 1:
        overall = fit_rate(errors, dts)
        print(f"\n  Fitted rate: {overall:.2f}  (need >= {expected_order - tol:.1f})")
        rate_passed = overall >= expected_order - tol
    else:
        print("\n  (need >= 2 CFL values to compute rate)")
        rate_passed = True

    cons_passed = print_conservation_check(all_cons, CONS_VARS, CONS_TOL)
    passed = rate_passed and cons_passed
    print(f"  {'PASS' if passed else 'FAIL'}")
    return passed


def main():
    parser = argparse.ArgumentParser(description="MFC RK3 temporal order verification")
    parser.add_argument("--cfls", type=float, nargs="+", default=None, help="Override per-scheme CFLs")
    parser.add_argument("--schemes", nargs="+", default=[s[0] for s in SCHEMES])
    parser.add_argument("--num-ranks", type=int, default=1)
    args = parser.parse_args()

    results = {}
    for label, extra_args, expected_order, tol, default_cfls in SCHEMES:
        if label not in args.schemes:
            continue
        cfls = args.cfls if args.cfls is not None else default_cfls
        results[label] = run_with_traceback(label, test_scheme, label, extra_args, expected_order, tol, cfls, args.num_ranks)

    sys.exit(0 if print_summary(results, label_width=14) else 1)


if __name__ == "__main__":
    main()
