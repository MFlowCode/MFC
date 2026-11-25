#!/usr/bin/env python3
"""
Interface advection (2D_advection): pressure-velocity equilibrium and interface error of a water disk in air.

    check_advection.py CASE_DIR [CASE_DIR ...] [--case-args ARGS ...] [--tol-p TOL]

Each CASE_DIR is a finished run (pre_process, simulation, post_process). The exact solution keeps p = p0 and
u = v = U and, after a whole number of flow-throughs, returns the volume fraction to its initial field. Prints, per run,
  dp, du    max |p/p0 - 1| and max |u - U|/U over all saves and cells (equilibrium across the interface)
  L1        sum |alpha - alpha_0| dA / (2 pi R) at the last save, in cells (an interface displacement)
  width     mixed cells (0.01 < alpha < 0.99) per unit interface length, in cells, at the first and the last save
and, for several runs, the observed order of L1. Exits 1 if dp of any run exceeds --tol-p.
"""

import argparse
import json
import math
import subprocess
import sys

import numpy as np

from mfc.viz import assemble_silo, discover_timesteps


def analyze(d, case_args):
    case = json.loads(subprocess.run([sys.executable, "case.py", *case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)
    p0, U, R = case["patch_icpp(1)%pres"], case["patch_icpp(1)%vel(1)"], case["patch_icpp(2)%radius"]
    steps = discover_timesteps(d, "silo")
    dp = du = 0.0
    for s in steps:
        V = {k: np.asarray(v).squeeze() for k, v in assemble_silo(d, s).variables.items()}
        dp = max(dp, float(np.abs(V["pres"] / p0 - 1.0).max()))
        du = max(du, float(max(np.abs(V["vel1"] - U).max(), np.abs(V["vel2"] - U).max())) / U)
        if s == steps[0]:
            a0 = V["alpha1"]
    a1 = V["alpha1"]
    dx = (case["x_domain%end"] - case["x_domain%beg"]) / (case["m"] + 1)
    length = 2.0 * math.pi * R / dx  # interface length in cells
    width = [float(((a > 0.01) & (a < 0.99)).sum()) / length for a in (a0, a1)]
    return dict(N=case["m"] + 1, t=steps[-1] * case["dt"], dp=dp, du=du, l1=float(np.abs(a1 - a0).sum()) / length, w0=width[0], w1=width[1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dirs", nargs="+")
    ap.add_argument("--case-args", nargs="*", default=None, help="case.py arguments of each run (default: none).")
    ap.add_argument("--tol-p", type=float, default=0.0, help="Fail if max |p/p0 - 1| exceeds this (0: no check).")
    args = ap.parse_args()
    cargs = args.case_args or [""] * len(args.case_dirs)
    res = sorted((analyze(d, c) for d, c in zip(args.case_dirs, cargs)), key=lambda r: r["N"])
    print(f"{'N':>5} {'t':>9} {'max dp/p0':>10} {'max du/U':>10} {'L1 [cells]':>11} {'width0':>7} {'width':>7}")
    for k, r in enumerate(res):
        order = f"  order {math.log(res[k - 1]['l1'] / res[k - 1]['N'] / (r['l1'] / r['N'])) / math.log(r['N'] / res[k - 1]['N']):.2f}" if k else ""
        print(f"{r['N']:5d} {r['t']:9.3e} {r['dp']:10.2e} {r['du']:10.2e} {r['l1']:11.4f} {r['w0']:7.2f} {r['w1']:7.2f}{order}")
    ok = args.tol_p <= 0 or max(r["dp"] for r in res) <= args.tol_p
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
