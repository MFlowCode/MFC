#!/usr/bin/env python3
"""
Compare a 2D_ibm_taylor_couette run with the exact steady solution u_theta = A r + B / r.

    check_couette.py CASE_DIR [--tol TOL]

Reads the serial text output of the last save (CASE_DIR/D/cons.*.dat), CASE_DIR/D/ib_forces.dat and the case parameters from
CASE_DIR/pre_process.inp and simulation.inp. Exact (incompressible, constant mu, inner radius R = 1, outer radius r2 = 2 R):
A = -Omega R^2 / (r2^2 - R^2), B = Omega R^2 r2^2 / (r2^2 - R^2), wall shear mu S with S = r d(u_theta / r)/dr = -2 B / r^2, torque
on the inner body per unit depth T = -4 pi mu B (the outer wall receives +4 pi mu B).

Quantities: L2 = rms of (u_theta - exact) / (Omega R) over the fluid cells (centre more than half a cell from both walls); its
azimuthal-mean part (error averaged in bins of one cell width in r) and the rest; the torque on each wall from the face fluxes
(mean over the last 10 % of the force history) relative to the exact value; the wall angular velocity from a cubic fit to the
cells within 0.5 to 4 cells of the inner wall, relative to Omega; the inner-wall shear from local cubic fits at 96 wall points
(rms relative error, the same fit applied to the exact profile, and their difference, "net"); the change of u_theta over the
last save interval. Exits 1 if L2 exceeds --tol.
"""

import argparse
import glob
import math
import os
import re
import sys

import numpy as np

R, R2 = 1.0, 2.0


def read_inp(case_dir):
    """Numeric entries of pre_process.inp and simulation.inp as a dict (the domain extents are only in the first)."""
    out = {}
    for name in ("pre_process.inp", "simulation.inp"):
        for line in open(os.path.join(case_dir, name)):
            m = re.match(r"\s*([\w%().,]+)\s*=\s*([-+0-9.eEdD]+)\s*$", line)
            if m:
                out[m.group(1)] = float(m.group(2).replace("d", "e").replace("D", "e"))
    return out


def saves(case_dir):
    """Time-step numbers of the saved fields, in order."""
    return sorted({int(re.search(r"cons\.1\.00\.(\d+)\.dat", f).group(1)) for f in glob.glob(os.path.join(case_dir, "D", "cons.1.00.*.dat"))})


def read_cons(case_dir, step, nvar):
    """Conservative variables 1..nvar as columns (x, y, q) of the text output; x and y are right cell faces (see cell_centres)."""
    return [np.loadtxt(glob.glob(os.path.join(case_dir, "D", f"cons.{k}.00.{step:06d}.dat"))[0]) for k in range(1, nvar + 1)]


def read_ib_history(case_dir):
    """D/ib_forces.dat as {ib id: array of (t, Fx, Fy, Fz, Tx, Ty, Tz, ...)}."""
    rows = {}
    for line in open(os.path.join(case_dir, "D", "ib_forces.dat")):
        v = line.split()
        if len(v) >= 8 and "\x00" not in line:
            rows.setdefault(int(v[0]), []).append([float(t) for t in v[1:]])
    return {k: np.array(a) for k, a in rows.items()}


def last_tenth(h):
    """Rows of a force history in the last 10 % of its time span."""
    return h[h[:, 0] > h[-1, 0] - 0.1 * (h[-1, 0] - h[0, 0])]


def velocity_field(case_dir, step, dx):
    """x, y of the cell centres (the text output labels each cell with its right face), rho, u, v."""
    rho, mx, my = read_cons(case_dir, step, 3)
    return rho[:, 0] - 0.5 * dx, rho[:, 1] - 0.5 * dx, rho[:, 2], mx[:, 2] / rho[:, 2], my[:, 2] / rho[:, 2]


BASIS = [(i, j) for i in range(4) for j in range(4) if i + j <= 3]  # cubic in the wall-normal and tangential distances


def inner_wall_shear(x, y, ut, ut_exact, dx, s_exact, n_points=96):
    """Shear S = r d(u_theta/r)/dr at the inner wall from local cubic least-squares fits to the fluid cells within 3.5 to 5 cells of
    each of n_points wall points. Returns the rms over the points of (S - S_exact)/|S_exact| for the run and for the exact
    profile (the fit's own error), and of their difference (the solution's own error; the fit is linear)."""
    r, out = np.hypot(x, y), {"run": [], "exact": []}
    fluid = (r > R) & (r < R2)
    for a in (np.arange(n_points) + 0.5) * 2.0 * np.pi / n_points:
        nx, ny = np.cos(a), np.sin(a)
        xw, yw = R * nx, R * ny
        for rho in (3.5, 4.0, 4.5, 5.0):  # grow the fit radius until the fit is well conditioned
            m = fluid & (np.hypot(x - xw, y - yw) <= rho * dx)
            if m.sum() < len(BASIS) + 4:
                continue
            s_, t_ = ((x[m] - xw) * nx + (y[m] - yw) * ny) / dx, (-(x[m] - xw) * ny + (y[m] - yw) * nx) / dx
            design = np.stack([s_**i * t_**j for i, j in BASIS], axis=1)
            if np.linalg.cond(design) < 1e5:
                break
        else:
            continue
        for key, field in (("run", ut), ("exact", ut_exact)):
            c = np.linalg.lstsq(design, field[m], rcond=None)[0]
            out[key].append((c[BASIS.index((1, 0))] / dx - c[0] / R - s_exact) / abs(s_exact))
    run, exact = np.array(out["run"]), np.array(out["exact"])
    rms = lambda v: math.sqrt(np.mean(v**2))
    return rms(run), rms(exact), rms(run - exact)


def analyze(case_dir):
    p = read_inp(case_dir)
    dx = (p["x_domain%end"] - p["x_domain%beg"]) / (p["m"] + 1)
    om, nu = p["patch_ib(1)%angular_vel(3)"], 1.0 / p["fluid_pp(1)%Re(1)"]
    den = R2**2 - R**2
    A, B = -om * R**2 / den, om * R**2 * R2**2 / den
    t_exact = -4.0 * math.pi * nu * B
    st = saves(case_dir)
    x, y, rho, u, v = velocity_field(case_dir, st[-1], dx)
    r = np.hypot(x, y)
    ut = (-u * y + v * x) / np.maximum(r, 1e-30)
    fluid = (r > R + 0.5 * dx) & (r < R2 - 0.5 * dx)
    err = (ut - (A * r + B / np.maximum(r, 1e-30))) / (om * R)
    out = {"cells_per_R": round(R / dx, 6), "step": st[-1], "L2": math.sqrt(np.mean(err[fluid] ** 2))}
    # azimuthal mean of the error in bins of one cell width, and what is left
    nb = round((R2 - R) / dx)
    idx = np.digitize(r, np.linspace(R, R2, nb + 1)) - 1
    ok = fluid & (idx >= 0) & (idx < nb)
    prof = np.array([err[ok & (idx == b)].mean() if (ok & (idx == b)).any() else np.nan for b in range(nb)])
    out["L2_mean"] = math.sqrt(np.nanmean(prof**2))
    out["L2_rest"] = math.sqrt(np.mean((err[ok] - prof[idx[ok]]) ** 2))
    # inner wall: angular velocity from a cubic fit to the azimuthal-mean cells within 0.5 to 4 cells of it, and the local-fit shear
    w = ut / np.maximum(r, 1e-30)
    near = (r - R > 0.5 * dx) & (r - R < 4.0 * dx)
    out["w_wall_err"] = (np.polynomial.polynomial.polyfit(r[near] - R, w[near], 3)[0] - om) / om
    out["shear_rms"], out["shear_fit_rms"], out["shear_net_rms"] = inner_wall_shear(x, y, ut, A * r + B / np.maximum(r, 1e-30), dx, -2.0 * B / R**2)
    if len(st) > 1:
        _, _, _, u0, v0 = velocity_field(case_dir, st[-2], dx)
        ut0 = (-u0 * y + v0 * x) / np.maximum(r, 1e-30)
        out["change_last_save"] = math.sqrt(np.mean(((ut - ut0)[fluid] / (om * R)) ** 2))
    h = read_ib_history(case_dir)
    out["T_inner_err"] = last_tenth(h[1])[:, 6].mean() / t_exact - 1.0
    out["T_outer_err"] = last_tenth(h[2])[:, 6].mean() / -t_exact - 1.0
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--tol", type=float, default=float("inf"), help="Exit 1 if L2 / (Omega R) exceeds this.")
    a = ap.parse_args()
    o = analyze(a.case_dir)
    print(f"cells per R {o['cells_per_R']:g}, step {o['step']}")
    print(f"  L2 / (Omega R)          {o['L2']:.3e}  (azimuthal mean {o['L2_mean']:.3e}, rest {o['L2_rest']:.3e})")
    print(f"  torque error            inner {o['T_inner_err']:+.3e}  outer {o['T_outer_err']:+.3e}")
    print(f"  inner wall speed error  {o['w_wall_err']:+.3e}  (cubic fit)")
    print(f"  inner wall shear error  rms {o['shear_rms']:.3e} (local cubic fits at 96 wall points), the fit's own error {o['shear_fit_rms']:.3e}, net {o['shear_net_rms']:.3e}")
    if "change_last_save" in o:
        print(f"  change over the last save interval {o['change_last_save']:.2e} (steady state when small)")
    sys.exit(0 if o["L2"] <= a.tol else 1)
