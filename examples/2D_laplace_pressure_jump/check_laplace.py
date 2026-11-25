#!/usr/bin/env python3
"""
Laplace pressure jump of the static 2D drop (2D_laplace_pressure_jump) against sigma/R.

    check_laplace.py CASE_DIR [--case-args '-N 50'] [--avg 0.2] [--tol TOL]

CASE_DIR is a finished run (pre_process, simulation, post_process). The drop is a quarter circle centred on the corner
(symmetry planes at x = 0 and y = 0). Over the saves in the last --avg fraction of the run it averages
  dp     mean pressure where alpha_water > 1 - 1e-3 minus the mean where alpha_water < 1e-3 and r > 1.5 R
  R      equivalent radius of the water area, sqrt(4 A / pi), A = sum(alpha_water) dx dy
and prints dp against sigma/R (and against sigma/r0 with the initial radius), the largest velocity (spurious currents)
and the capillary number mu-free scale |u|_max sqrt(rho_l R / sigma). Exits 1 if |dp R / sigma - 1| exceeds --tol.
"""

import argparse
import json
import math
import subprocess
import sys

import numpy as np

from mfc.viz import assemble_silo, discover_timesteps


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--case-args", default="", help="Arguments given to case.py for the run.")
    ap.add_argument("--avg", type=float, default=0.2, help="Averaging window, fraction of the saves (default 0.2).")
    ap.add_argument("--tol", type=float, default=0.0, help="Fail if |dp R/sigma - 1| exceeds this (0: no check).")
    args = ap.parse_args()

    d = args.case_dir
    case = json.loads(subprocess.run([sys.executable, "case.py", *args.case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)
    sigma, r0 = case["sigma"], case["patch_icpp(2)%radius"]
    rho_l = case["patch_icpp(2)%alpha_rho(1)"] / case["patch_icpp(2)%alpha(1)"]
    steps = discover_timesteps(d, "silo")
    use = steps[-max(1, int(args.avg * len(steps))) :]
    dps, rs, umax = [], [], []
    for s in use:
        data = assemble_silo(d, s)
        x, y = (np.asarray(c).ravel() for c in (data.x_cc, data.y_cc))
        V = {k: np.asarray(v).squeeze() for k, v in data.variables.items()}
        X, Y = np.meshgrid(x, y, indexing="ij")
        a, p = V["alpha1"], V["pres"]
        R = math.sqrt(4.0 * a.sum() * (x[1] - x[0]) * (y[1] - y[0]) / math.pi)
        inside, outside = a > 1.0 - 1e-3, (a < 1e-3) & (np.hypot(X, Y) > 1.5 * R)
        dps.append(p[inside].mean() - p[outside].mean())
        rs.append(R)
        umax.append(float(np.sqrt(V["vel1"] ** 2 + V["vel2"] ** 2).max()))
    dp, R, u = float(np.mean(dps)), float(np.mean(rs)), float(np.max(umax))
    err = dp * R / sigma - 1.0
    print(f"N = {case['m'] + 1}, t = {use[-1] * case['dt']:.4f} s, {len(use)} saves averaged; R = {R:.5f} (r0 = {r0})")
    print(f"dp = {dp:.4f} Pa, sigma/R = {sigma / R:.4f} Pa ({err:+.3%}), sigma/r0 = {sigma / r0:.4f} Pa ({dp * r0 / sigma - 1:+.3%})")
    print(f"spread of dp over the window {np.ptp(dps) / dp:.2e}; max |u| = {u:.3e} m/s ({u * math.sqrt(rho_l * R / sigma):.3e} sqrt(sigma/(rho_l R)))")
    ok = args.tol <= 0 or abs(err) <= args.tol
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
