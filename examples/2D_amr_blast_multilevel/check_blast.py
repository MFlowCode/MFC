#!/usr/bin/env python3
"""
Density of the circular blast (case.py, with or without --no-amr) against a fine 1D radial reference.

    check_blast.py RUN_DIR [RUN_DIR ...] [--ncell 16000] [--tol TOL]

The reference solves the radially symmetric Euler equations, d(r U)/dt + d(r F)/dr = (0, p, 0), with a second-order finite-volume
scheme (MUSCL with the minmod limiter, HLL fluxes, SSP-RK2) on --ncell cells over 0 <= r <= 0.75, from the run's disc (radius,
states and time read from its input files). Each run's last save (restart_data/lustre_<step>.dat or p_all/; level 0 for AMR runs) is compared
with the reference averaged over each cell (4 x 4 samples): the L1 density error over the box and the radius of the outgoing shock,
the outermost radius where the azimuthal-mean density is above (rho_ambient + rho_peak) / 2, measured the same way on the
cell-averaged reference. Exit 1 if the L1 error of the last
run exceeds --tol.
"""

import argparse
import functools
import glob
import os
import re
import sys

import numpy as np


def read_inp(d):
    out = {}
    for name in ("pre_process.inp", "simulation.inp"):
        for line in open(os.path.join(d, name)):
            match = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
            if match:
                out[match.group(1)] = match.group(2).strip("'\",")
    return out


def saves(d, nvar, shape):
    """{step: conservative fields 1..nvar} of a one-rank or parallel_io run (restart_data/lustre_<step>.dat or p_all/p0/<step>/)."""
    n = int(np.prod(shape))
    out = {}
    for f in glob.glob(os.path.join(d, "restart_data", "lustre_*.dat")):
        if re.search(r"lustre_\d+\.dat$", f):
            out[int(re.findall(r"\d+", os.path.basename(f))[0])] = np.fromfile(f, dtype=np.float64)[: nvar * n].reshape(nvar, *shape)
    if out:
        return out
    if len(glob.glob(os.path.join(d, "p_all", "p*"))) != 1:
        sys.exit(f"{d}: run on one rank or with parallel_io = T")
    for s in glob.glob(os.path.join(d, "p_all", "p0", "[0-9]*")):
        q = []
        for i in range(1, nvar + 1):
            with open(os.path.join(s, f"q_cons_vf{i}.dat"), "rb") as f:
                f.read(4)
                q.append(np.fromfile(f, dtype=np.float64, count=n).reshape(shape))
        out[int(os.path.basename(s))] = np.array(q)
    return out


@functools.lru_cache
def radial_reference(r0, inner, outer, g, t_end, ncell, r_max=0.75, cfl=0.4):
    """(r centres, rho, u, p) of the cylindrical Riemann problem at t_end."""
    edges = np.linspace(0.0, r_max, ncell + 1)
    rc, dr = 0.5 * (edges[:-1] + edges[1:]), r_max / ncell
    vol = 0.5 * (edges[1:] ** 2 - edges[:-1] ** 2)
    w = np.where(rc < r0, 1.0, 0.0)
    rho = w * inner[0] + (1 - w) * outer[0]
    p = w * inner[2] + (1 - w) * outer[2]
    q = np.array([rho, np.zeros(ncell), p / (g - 1.0)])

    def prim(q):
        rho = q[0]
        u = q[1] / rho
        return rho, u, (g - 1.0) * (q[2] - 0.5 * rho * u**2)

    def flux(rho, u, p):
        return np.array([rho * u, rho * u**2 + p, (p / (g - 1.0) + 0.5 * rho * u**2 + p) * u])

    def rhs(q):
        w = np.array(prim(q))
        wg = np.concatenate([w[:, 1::-1] * np.array([[1], [-1], [1]]), w, w[:, -1:], w[:, -1:]], axis=1)  # reflect at r = 0
        dl, dr_ = wg[:, 1:-1] - wg[:, :-2], wg[:, 2:] - wg[:, 1:-1]
        slope = np.where(dl * dr_ > 0, np.sign(dl) * np.minimum(abs(dl), abs(dr_)), 0.0)
        wl = wg[:, 1:-2] + 0.5 * slope[:, :-1]  # left state at faces 0 .. ncell
        wr = wg[:, 2:-1] - 0.5 * slope[:, 1:]
        fl, fr = flux(*wl), flux(*wr)
        ql = np.array([wl[0], wl[0] * wl[1], wl[2] / (g - 1.0) + 0.5 * wl[0] * wl[1] ** 2])
        qr = np.array([wr[0], wr[0] * wr[1], wr[2] / (g - 1.0) + 0.5 * wr[0] * wr[1] ** 2])
        cl, cr = np.sqrt(g * wl[2] / wl[0]), np.sqrt(g * wr[2] / wr[0])
        sl, sr = np.minimum(wl[1] - cl, wr[1] - cr), np.maximum(wl[1] + cl, wr[1] + cr)
        f = np.where(sl >= 0, fl, np.where(sr <= 0, fr, (sr * fl - sl * fr + sl * sr * (qr - ql)) / (sr - sl)))
        f = f * edges  # face area r
        src = np.zeros_like(q)
        src[1] = w[2] * dr  # integral of p dA/dr over the cell, A = r
        return (-(f[:, 1:] - f[:, :-1]) + src) / vol, np.max(np.abs(w[1]) + np.sqrt(g * w[2] / w[0]))

    t = 0.0
    while t < t_end:
        k1, smax = rhs(q)
        dt = min(cfl * dr / smax, t_end - t)
        q1 = q + dt * k1
        q = 0.5 * (q + q1 + dt * rhs(q1)[0])
        t += dt
    return (rc, *prim(q))


def analyze(d, ncell):
    inp = read_inp(d)
    nx, ny = int(inp["m"]) + 1, int(inp["n"]) + 1
    x0, x1, y0, y1 = (float(inp[k]) for k in ("x_domain%beg", "x_domain%end", "y_domain%beg", "y_domain%end"))
    g = 1.0 + 1.0 / float(inp["fluid_pp(1)%gamma"])
    state = lambda k: tuple(float(inp[f"patch_icpp({k})%{v}"]) for v in ("alpha_rho(1)", "vel(1)", "pres"))  # noqa: E731
    xc, yc, r0 = (float(inp[f"patch_icpp(2)%{v}"]) for v in ("x_centroid", "y_centroid", "radius"))
    data = saves(d, 1, (ny, nx))
    step = max(data)
    rho = data[step][0]
    t = step * float(inp["dt"])
    rr, rho_ref, _, _ = radial_reference(r0, state(2), state(1), g, t, ncell)
    hx, hy = (x1 - x0) / nx, (y1 - y0) / ny
    sub = (np.arange(4) + 0.5) / 4
    xs = x0 + (np.arange(nx)[:, None] + sub[None, :]).ravel() * hx - xc
    ys = y0 + (np.arange(ny)[:, None] + sub[None, :]).ravel() * hy - yc
    r = np.hypot(xs[None, :], ys[:, None])
    ex = np.interp(r, rr, rho_ref).reshape(ny, 4, nx, 4).mean(axis=(1, 3))
    l1 = np.abs(rho - ex).sum() * hx * hy

    def shock_radius(field):
        rcell = np.hypot(x0 + (np.arange(nx) + 0.5) * hx - xc, (y0 + (np.arange(ny) + 0.5) * hy - yc)[:, None])
        bins = np.arange(0.0, 0.5, max(hx, hy))
        idx = np.digitize(rcell.ravel(), bins)
        prof = np.array([field.ravel()[idx == i].mean() if np.any(idx == i) else np.nan for i in range(1, bins.size)])
        level = 0.5 * (state(1)[0] + np.nanmax(prof))
        above = np.where(prof > level)[0]
        i = above[-1]
        return bins[i] + (bins[i + 1] - bins[i]) * (prof[i] - level) / (prof[i] - prof[i + 1]) + 0.5 * (bins[1] - bins[0])

    return dict(name=os.path.basename(os.path.abspath(d)), n=nx, t=t, l1=l1, rs=shock_radius(rho), rs_ref=shock_radius(ex), amr=inp.get("amr", "F"))


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("runs", nargs="+", help="Run directories.")
parser.add_argument("--ncell", type=int, default=16000, help="Reference cells over 0 <= r <= 0.75 (default 16000).")
parser.add_argument("--tol", type=float, default=None, help="Exit 1 if the L1 density error of the last run exceeds this.")
args = parser.parse_args()

res = [analyze(d, args.ncell) for d in args.runs]
print(f"{'run':>24} {'N':>5} {'amr':>3} {'t':>8} {'L1 rho':>10} {'r_shock':>8} {'reference':>9}")
for r in res:
    print(f"{r['name']:>24} {r['n']:>5} {r['amr']:>3} {r['t']:>8.5f} {r['l1']:>10.3e} {r['rs']:>8.4f} {r['rs_ref']:>9.4f}")
if args.tol is not None:
    ok = res[-1]["l1"] <= args.tol
    print(f"{'PASS' if ok else 'FAIL'}: L1 density error {res[-1]['l1']:.3e} (tolerance {args.tol:g})")
    sys.exit(0 if ok else 1)
