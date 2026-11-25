#!/usr/bin/env python3
"""
Circularly polarized Alfven wave (1D_mhd_smooth_alfven_wave) against its exact solution.

    check_alfven.py RUN_DIR [RUN_DIR ...] [--tol TOL]

The wave v_perp = B_perp = 0.1 (sin 2 pi x, cos 2 pi x), B_x = 1, rho = 1, p = 0.1 is an exact nonlinear solution of
ideal MHD that travels at the Alfven speed B_x / sqrt(rho) = 1 in -x (v_perp = +B_perp / sqrt(rho)); |B_perp|, rho and p
stay constant. At the last save t, each RUN_DIR is compared with the initial state shifted by -t. Prints the L1 errors
(mean |error|) of v_y, v_z, B_y, B_z, the max deviation of rho and p, and, for a sweep, the observed order of the
B_y error. Exit 1 if the finest run's B_y L1 error exceeds --tol.
"""

import argparse
import math
import os
import re
import sys

import numpy as np

from mfc.viz import assemble, assemble_silo, discover_format, discover_timesteps


def inp(d):
    out = {}
    for line in open(os.path.join(d, "simulation.inp")):
        m = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
        if m:
            out[m.group(1)] = m.group(2).strip("'\",")
    return out


def check(d):
    fmt = discover_format(d)
    step = discover_timesteps(d, fmt)[-1]
    a = assemble_silo(d, step) if fmt == "silo" else assemble(d, step, fmt)
    t = step * float(inp(d)["dt"].replace("d", "e"))
    x, v = a.x_cc, a.variables
    ph = 2 * math.pi * (x + t)
    exact = {"vel2": 0.1 * np.sin(ph), "vel3": 0.1 * np.cos(ph), "By": 0.1 * np.sin(ph), "Bz": 0.1 * np.cos(ph)}
    l1 = {k: float(np.mean(np.abs(v[k] - e))) for k, e in exact.items()}
    drho = float(np.abs(v["rho"] - 1).max())
    dp = float(np.abs(v["pres"] - 0.1).max())
    print(f"{d}: N = {x.size}, t = {t:.6f}, " + ", ".join(f"L1({k}) = {e:.3e}" for k, e in l1.items()) + f", max|rho-1| = {drho:.2e}, max|p-0.1| = {dp:.2e}")
    return x.size, l1["By"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--tol", type=float, default=1e-3, help="Max L1(B_y) of the finest run")
    args = ap.parse_args()
    res = sorted(check(os.path.abspath(d)) for d in args.runs)
    for (n1, e1), (n2, e2) in ((a, b) for a, b in zip(res, res[1:]) if a[0] != b[0]):
        print(f"observed order {n1} -> {n2}: {math.log(e1 / e2) / math.log(n2 / n1):.2f}")
    if res[-1][1] > args.tol:
        print(f"FAIL: L1(B_y) = {res[-1][1]:.3e} > {args.tol}")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
