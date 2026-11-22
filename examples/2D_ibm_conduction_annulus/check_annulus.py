#!/usr/bin/env python3
"""
Errors of a 2D_ibm_conduction_annulus run against the exact steady solution T(r) = T1 + (T2 - T1) ln(r / R) / ln(r2 / R).

    check_annulus.py CASE_DIR [--r2 2] [--tol TOL] [--dtype float32]

Reads pre_process.inp and simulation.inp and the last restart_data/lustre_<k>.dat (parallel_io = T, one rank). Prints, relative to
T2 - T1 and with the sign (MFC - exact):
  wall error   a quadratic least-squares fit T(s) through the gas cells 0.5 dx < s < 3.5 dx from the wall (s the distance into
               the gas), extrapolated to s = 0, minus the wall temperature; one value per wall
  gradient     the dT/dr of the same fit at the wall, relative to the exact (T2 - T1) / (r ln(r2 / R)) (not a reliable measure: README)
  gap L2, Linf relative L2 and maximum error of T over the gas cells 0.5 dx from either wall
Exits 1 if the larger |wall error| exceeds --tol.
"""

import argparse
import glob
import math
import os
import re
import sys

import numpy as np

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("case_dir")
parser.add_argument("--r2", type=float, default=2.0, help="Outer wall radius in units of the inner radius.")
parser.add_argument("--dtype", default="float64", help="Storage type of restart_data (float32 for single-precision builds).")
parser.add_argument("--tol", type=float, default=None, help="Exit 1 if the larger |wall error| / (T2 - T1) exceeds this.")
args = parser.parse_args()


def read_inp(path):
    out = {}
    for line in open(path):
        match = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
        if match:
            out[match.group(1)] = match.group(2).strip("'\"").replace("d", "e")
    return out


inp = {**read_inp(os.path.join(args.case_dir, "pre_process.inp")), **read_inp(os.path.join(args.case_dir, "simulation.inp"))}
m, n = int(inp["m"]) + 1, int(inp["n"]) + 1
x0, x1, y0, y1 = (float(inp[k]) for k in ("x_domain%beg", "x_domain%end", "y_domain%beg", "y_domain%end"))
dx = (x1 - x0) / m
R = float(inp["patch_ib(1)%radius"])
r2 = args.r2 * R
T1, T2 = float(inp["patch_ib(1)%Twall"]), float(inp["patch_ib(2)%Twall"])
cv = float(inp["fluid_pp(1)%cv"])

found = [re.fullmatch(r"lustre_(\d+)\.dat", os.path.basename(f)) for f in glob.glob(os.path.join(args.case_dir, "restart_data", "lustre_*.dat"))]
last = max(int(f.group(1)) for f in found if f)
q = np.fromfile(os.path.join(args.case_dir, "restart_data", f"lustre_{last}.dat"), dtype=args.dtype).reshape(5, n, m).astype(float)
rho, mx, my, E = q[0], q[1], q[2], q[3]  # one fluid: alpha_rho, momentum, energy (indexed [y, x])
T = (E - 0.5 * (mx**2 + my**2) / rho) / (rho * cv)
X, Y = np.meshgrid(x0 + (np.arange(m) + 0.5) * dx, y0 + (np.arange(n) + 0.5) * dx)
r = np.hypot(X, Y)

exact = lambda rr: T1 + (T2 - T1) * np.log(rr / R) / np.log(r2 / R)
dT = T2 - T1
wall_err, grad_err = [], []
for rw, Tw, sign in ((R, T1, 1.0), (r2, T2, -1.0)):
    s = sign * (r - rw)
    f = (s > 0.5 * dx) & (s < 3.5 * dx)
    A = np.stack([np.ones(f.sum()), s[f], s[f] ** 2], 1)
    a, b, _ = np.linalg.lstsq(A, T[f], rcond=None)[0]
    wall_err.append((a - Tw) / dT)
    grad_err.append(sign * b / (dT / (rw * math.log(r2 / R))) - 1.0)

gap = (r > R + 0.5 * dx) & (r < r2 - 0.5 * dx)
err = (T - exact(np.clip(r, R, r2)))[gap] / dT
cells_per_R = round(R / dx)
print(
    f"cells/R {cells_per_R}  save {last}  wall error inner {wall_err[0]:+.3e} outer {wall_err[1]:+.3e}  gradient inner {grad_err[0]:+.3e} outer {grad_err[1]:+.3e}  gap L2 {math.sqrt(np.mean(err**2)):.3e}  Linf {np.abs(err).max():.3e}"
)
sys.exit(1 if args.tol is not None and max(abs(w) for w in wall_err) > args.tol else 0)
