#!/usr/bin/env python3
"""
Convergence-rate verification for MFC's 2D isentropic vortex problem.

Uses hcid=283: 3-pt Gauss-Legendre cell averages of conserved variables as IC.
The vortex strength eps=0.01 (set in case.py) is chosen so that the dominant
error source is the WENO spatial truncation error O(eps^2 * h^p), not the
primitive-to-conserved covariance floor O(eps^3 * h^2).  For h > eps^(1/3)=0.22
(i.e., N < 46 per dimension), the p-th order scheme shows rate p.

L2(rho(T) - rho(0)) measures accumulated scheme error; the comparison to rho(0)
(the numerical IC) eliminates IC discretisation error, isolating the scheme error.

WENO7/TENO7 are NOT tested here.  For the isentropic vortex, the IC
primitive→conserved covariance error is O(eps^3 * h^2).  The WENO7 scheme
error is O(eps^2 * h^7).  Scheme error dominates only when h > eps^(1/5);
with eps=0.01 that requires h > 0.40, i.e., N < 25.  At N=64-128 the
covariance floor dominates and the measured rate is ~2, not 7.
WENO7/TENO7 7th-order convergence is verified by the 1D test (run_convergence_1d.py)
which uses a pure advection problem that avoids this nonlinear floor.

Usage:
    python toolchain/mfc/test/run_convergence.py [--resolutions 32 64 128]
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

CASE = "examples/2D_isentropicvortex_convergence/case.py"
DOMAIN_LEN = 10.0  # vortex domain [-5, 5]

# (label, extra_args, expected_order, tolerance, min_N, max_N)
# With eps=0.01 and N=32..128 the prim->cons covariance error O(eps^3 h^2) is
# well below the scheme's spatial error O(eps^2 h^p), so each scheme shows its
# nominal rate.  WENO5/TENO5: min_N=64 satisfies the 25-cell/rank stencil
# requirement under 4 MPI ranks (2x2). WENO3: pre-asymptotic at N=32..128
# (~2.0-2.2; approaches 3 at finer grids). WENO7/TENO7 are omitted (covariance
# floor dominates at testable N — see module docstring).
SCHEMES = [
    ("WENO5", ["--order", "5"], 5, 1.0, 64, None),
    ("WENO3", ["--order", "3"], 3, 1.2, 32, None),
    ("WENO1", ["--order", "1"], 1, 0.4, 32, None),
    ("MUSCL2", ["--muscl"], 2, 0.5, 32, None),
    ("TENO5", ["--order", "5", "--teno", "--teno-ct", "1e-6"], 5, 1.0, 64, None),
]

# 2D single-fluid Euler: vf1=ρ, vf2=ρu, vf3=ρv, vf4=E. Momentum is excluded:
# the vortex has zero net linear momentum, making the relative error
# ill-conditioned. Density and energy have large nonzero integrals.
CONS_VARS = [("density", 1), ("energy", 4)]


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
            dx = DOMAIN_LEN / N
            cfg, run_dir = run_mfc_case(CASE, tmpdir, f"N{N}", ["-N", str(N)] + extra_args, num_ranks)
            Nt = int(cfg["t_step_stop"])
            nts.append(Nt)
            rho0 = read_cons_var(run_dir, 0, 1, num_ranks, expected_size=N * N)
            rhoT = read_cons_var(run_dir, Nt, 1, num_ranks, expected_size=N * N)
            errors.append(l2_norm(rhoT - rho0, dx**2))
            all_cons.append(conservation_errors(run_dir, Nt, dx**2, CONS_VARS, num_ranks, expected_size=N * N))

    dxs = [DOMAIN_LEN / N for N in resolutions]
    rates = pairwise_rates(errors, dxs)

    print(f"\n  {'N':>6}  {'Nt':>5}  {'dx':>10}  {'L2 error':>14}  {'rate':>8}")
    print(f"  {'-' * 6}  {'-' * 5}  {'-' * 10}  {'-' * 14}  {'-' * 8}")
    for i, N in enumerate(resolutions):
        r_str = f"{rates[i]:>8.2f}" if rates[i] is not None else f"{'---':>8}"
        print(f"  {N:>6}  {nts[i]:>5}  {dxs[i]:>10.5f}  {errors[i]:>14.6e}  {r_str}")

    if len(resolutions) > 1:
        overall = fit_rate(errors, dxs)
        print(f"\n  Fitted rate: {overall:.2f}  (need >= {expected_order - tol:.1f})")
        rate_passed = overall >= expected_order - tol
    else:
        print("\n  (need >= 2 resolutions to compute rate)")
        rate_passed = True

    cons_passed = print_conservation_check(all_cons, CONS_VARS, CONS_TOL)
    passed = rate_passed and cons_passed
    print(f"  {'PASS' if passed else 'FAIL'}")
    return passed


def main():
    parser = argparse.ArgumentParser(description="MFC convergence-rate verification")
    parser.add_argument("--resolutions", type=int, nargs="+", default=[32, 64, 128])
    parser.add_argument("--schemes", nargs="+", default=[s[0] for s in SCHEMES])
    parser.add_argument("--num-ranks", type=int, default=1)
    args = parser.parse_args()

    results = {}
    for label, extra_args, expected_order, tol, min_N, max_N in SCHEMES:
        if label not in args.schemes:
            continue
        results[label] = run_with_traceback(label, test_scheme, label, extra_args, expected_order, tol, args.resolutions, min_N, max_N, args.num_ranks)

    sys.exit(0 if print_summary(results, label_width=12) else 1)


if __name__ == "__main__":
    main()
