#!/usr/bin/env python3
"""
Compare a 2D_ibm_plane_poiseuille run with u = g y (H - y) / (2 nu).

    check_poiseuille.py CASE_DIR [--tol TOL]

Reads the serial text output of the last save (CASE_DIR/D/cons.*.dat), CASE_DIR/D/ib_forces.dat and the parameters from
CASE_DIR/pre_process.inp and simulation.inp (H = 1; the walls are at y = 0 and y = H; the wall offset does not enter the reference).
Quantities, all relative to the exact value: L2 = rms of u - exact over the fluid cells, divided by u_max = g H^2 / (8 nu); Q = flow
rate (sum of u dy over the fluid cells) against g H^3 / (12 nu); F = wall force from the IB face fluxes (mean over the last 10 % of
the history) against g H Lx / 2 per wall (in steady state a body-force balance, so it does not separate the wall orders); u_wall =
velocity extrapolated to the wall by a quadratic fit through the three cells nearest it, over u_max (exact 0); the change of u over
the last save interval. Exits 1 if L2 exceeds --tol.
"""

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "2D_ibm_taylor_couette"))
from check_couette import last_tenth, read_ib_history, read_inp, saves, velocity_field  # noqa: E402

H = 1.0


def analyze(case_dir):
    p = read_inp(case_dir)
    dx = (p["y_domain%end"] - p["y_domain%beg"]) / (p["n"] + 1)
    lx = p["x_domain%end"] - p["x_domain%beg"]
    nu, g = 1.0 / p["fluid_pp(1)%Re(1)"], p["g_x"]
    u_max = g * H**2 / (8.0 * nu)
    st = saves(case_dir)

    def row_means(step):
        x, y, rho, u, v = velocity_field(case_dir, step, dx)
        rows = np.unique(np.round(y, 9))
        return rows, np.array([u[np.abs(y - r) < 1e-9].mean() for r in rows]), y, u

    rows, um, y, u = row_means(st[-1])
    fluid = (rows > 0) & (rows < H)
    exact = g * rows * (H - rows) / (2.0 * nu)
    out = {"cells": round(H / dx, 6), "step": st[-1]}
    out["L2"] = math.sqrt(np.mean(((um - exact)[fluid] / u_max) ** 2))
    out["Q_err"] = um[fluid].sum() * dx / (g * H**3 / (12.0 * nu)) - 1.0
    yf, uf = rows[fluid], um[fluid]
    for name, d, v in (("bottom", yf, uf), ("top", H - yf, uf)):
        o = np.argsort(d)[:3]
        out[f"u_wall_{name}"] = np.polynomial.polynomial.polyfit(d[o], v[o], 2)[0] / u_max
    if len(st) > 1:
        out["change_last_save"] = math.sqrt(np.mean(((um - row_means(st[-2])[1])[fluid] / u_max) ** 2))
    h = read_ib_history(case_dir)
    for i, name in ((1, "bottom"), (2, "top")):
        out[f"F_{name}_err"] = last_tenth(h[i])[:, 1].mean() / (g * H * lx / 2.0) - 1.0
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--tol", type=float, default=float("inf"), help="Exit 1 if L2 / u_max exceeds this.")
    a = ap.parse_args()
    o = analyze(a.case_dir)
    print(f"cells across the gap {o['cells']:g}, step {o['step']}")
    print(f"  L2 / u_max              {o['L2']:.3e}")
    print(f"  flow rate error         {o['Q_err']:+.3e}")
    print(f"  wall velocity / u_max   bottom {o['u_wall_bottom']:+.3e}  top {o['u_wall_top']:+.3e}  (quadratic fit)")
    print(f"  wall force error        bottom {o['F_bottom_err']:+.3e}  top {o['F_top_err']:+.3e}")
    if "change_last_save" in o:
        print(f"  change over the last save interval {o['change_last_save']:.2e} (steady state when small)")
    sys.exit(0 if o["L2"] <= a.tol else 1)
