#!/usr/bin/env python3
"""
Check a run of 2D_ibm_slip_cylinder (run pre_process, simulation and post_process in CASE_DIR first).

    check_slip.py CASE_DIR [--case-args '--steps 0 --cells_per_R 16'] [--tol TOL]

leak mode: total gas mass over the fluid cells (centres outside the body), relative to the first save, at every save.
           Exits 1 if the last |change| exceeds --tol (when given).
vortex mode: rms over wall-distance bands (in cells) of (p - p_exact)/(0.5 rho0 (M c0)^2), (u_theta - u_theta,exact)/(M c0) and
           u_r/(M c0), at the last save, against the exact steady isentropic vortex.
"""

import argparse
import json
import subprocess
import sys

import numpy as np

from mfc.viz import assemble_silo, discover_timesteps

GAM, RHO0, P0 = 1.4, 1.2, 1.0e5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--case-args", default="", help="Arguments given to case.py for the run.")
    ap.add_argument("--tol", type=float, default=None, help="leak mode: largest accepted |mass change| at the last save.")
    args = ap.parse_args()

    d = args.case_dir
    cargs = args.case_args.split()
    case = json.loads(subprocess.run([sys.executable, "case.py", *cargs], cwd=d, capture_output=True, text=True, check=True).stdout)
    mode = cargs[cargs.index("--mode") + 1] if "--mode" in cargs else "leak"
    mach = float(cargs[cargs.index("--mach") + 1]) if "--mach" in cargs else (0.3 if mode == "vortex" else 0.5)
    R = case["patch_ib(1)%radius"]
    steps = discover_timesteps(d, "silo")
    data = [assemble_silo(d, s) for s in steps]
    x, y = data[0].x_cc, data[0].y_cc
    dx = x[1] - x[0]
    X, Y = np.meshgrid(x, y, indexing="ij")
    r = np.hypot(X, Y)

    if mode == "leak":
        fluid = r > R
        mass = [float(np.sum(a.variables["rho"][fluid])) for a in data]
        change = [m / mass[0] - 1 for m in mass[1:]]
        print(f"{x.size} x {y.size} cells, {R / dx:.0f} cells/R, saves {steps}: mass change " + " ".join(f"{c:+.3e}" for c in change))
        sys.exit(1 if args.tol is not None and abs(change[-1]) > args.tol else 0)

    c0 = (GAM * P0 / RHO0) ** 0.5
    ut = mach * c0 * R / np.maximum(r, R)
    tr = 1.0 - 0.5 * (GAM - 1.0) * (ut / c0) ** 2
    p_ex = P0 * tr ** (GAM / (GAM - 1.0))
    a = data[-1]
    ux, uy = a.variables["vel1"], a.variables["vel2"]
    rr = np.maximum(r, 1e-30)
    u_t, u_r = (-Y * ux + X * uy) / rr, (X * ux + Y * uy) / rr
    qd = 0.5 * RHO0 * (mach * c0) ** 2
    s = (r - R) / dx
    print(f"{R / dx:.0f} cells/R, M = {mach}, save {steps[-1]}: rms error by wall-distance band (cells)")
    print("band       p        u_theta   u_r")
    for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 6)):
        mk = (s > lo) & (s <= hi)
        e_p = np.sqrt(np.mean(((a.variables["pres"] - p_ex)[mk] / qd) ** 2))
        e_t = np.sqrt(np.mean(((u_t - ut)[mk] / (mach * c0)) ** 2))
        e_r = np.sqrt(np.mean((u_r[mk] / (mach * c0)) ** 2))
        print(f"{lo}-{hi}    {e_p:.3e}  {e_t:.3e}  {e_r:.3e}")


if __name__ == "__main__":
    main()
