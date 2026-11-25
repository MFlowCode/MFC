#!/usr/bin/env python3
"""
Divergence error and mass conservation of the Orszag-Tang vortex (2D_orszag_tang_hyper_cleaning).

    check_divb.py RUN_DIR [CONTROL_DIR] [--tol TOL]

RUN_DIR is a finished run (pre_process, simulation, post_process); CONTROL_DIR, if given, is the same deck run with
--no-clean. For each save the script prints the normalized divergence error

    eps = mean |div B| dx / mean |B|,   div B by second-order central differences on the periodic grid,

its maximum max |div B| dx / max |B|, and the relative change of the total mass. With a control run it also prints the
ratio of the two eps at each save and the L1 difference of the densities at the last common save relative to the mean
density (without cleaning this deck fails with NaN after step 900, so run the control with --steps 800).
The initial field is divergence free (B = (-sin 2 pi y, sin 4 pi x) / sqrt(4 pi)), so eps measures what the scheme creates.
Exit 1 if eps of RUN_DIR at the last save exceeds --tol.
"""

import argparse
import os
import sys

import numpy as np

from mfc.viz import assemble, assemble_silo, discover_format, discover_timesteps


def load(d, step):
    fmt = discover_format(d)
    a = assemble_silo(d, step) if fmt == "silo" else assemble(d, step, fmt)
    return a.x_cc, a.y_cc, {k: np.asarray(a.variables[k]) for k in ("rho", "Bx", "By")}


def series(d):
    out = []
    steps = discover_timesteps(d, discover_format(d))
    for s in steps:
        x, y, v = load(d, s)
        dx, dy = x[1] - x[0], y[1] - y[0]
        div = (np.roll(v["Bx"], -1, 0) - np.roll(v["Bx"], 1, 0)) / (2 * dx) + (np.roll(v["By"], -1, 1) - np.roll(v["By"], 1, 1)) / (2 * dy)
        bmag = np.hypot(v["Bx"], v["By"])
        out.append((s, np.mean(np.abs(div)) * dx / np.mean(bmag), np.abs(div).max() * dx / bmag.max(), v["rho"].sum(), v["rho"]))
    m0 = out[0][3]
    return [(s, e, emax, mass / m0 - 1, rho) for s, e, emax, mass, rho in out]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("control", nargs="?")
    ap.add_argument("--tol", type=float, default=1e-2, help="Max eps at the last save")
    args = ap.parse_args()
    a = series(os.path.abspath(args.run))
    b = series(os.path.abspath(args.control)) if args.control else None
    print("step      eps      max eps   mass change" + ("   control eps   control/run" if b else ""))
    for i, (s, e, emax, dm, _) in enumerate(a):
        row = f"{s:5d}  {e:.3e}  {emax:.3e}  {dm:+.2e}"
        if b and i < len(b):
            row += f"     {b[i][1]:.3e}     {b[i][1] / e if e > 0 else float('nan'):8.2f}"
        print(row)
    if b:
        k = min(len(a), len(b)) - 1  # last common save (the control may stop early)
        ra, rb = a[k][4], b[k][4]
        print(f"L1(rho_run - rho_control) / mean rho at step {a[k][0]}: {np.mean(np.abs(ra - rb)) / np.mean(rb):.3e}")
    if a[-1][1] > args.tol:
        print(f"FAIL: eps = {a[-1][1]:.3e} > {args.tol}")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
