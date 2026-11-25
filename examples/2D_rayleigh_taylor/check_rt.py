#!/usr/bin/env python3
"""
Early growth of the viscous Rayleigh-Taylor instability (2D_rayleigh_taylor) against linear theory.

    check_rt.py CASE_DIR [--case-args '-N 200 --tend 0.6'] [--ka-max 0.25] [--tol TOL]

CASE_DIR is a finished run (pre_process, simulation, post_process). The initial interface (hardcoded IC 204) is
y = h - a0 cos(k x), h = 0.7, lambda = 0.2, k a0 = 0.05, with the tanh profile alpha = (1 + tanh((y - y_i)/delta))/2,
delta = 2.5e-3, between rho = 3 (above) and 1. At each save the interface height eta(x) is the light-fluid volume per unit width in each
column and the amplitude a(t) is its cos(k x) Fourier coefficient. Over the saves with k a <= --ka-max it fits
  a(t) = A exp(n t) + B exp(-m t)
and prints n, and the local slope d ln a/dt at the end of the window. The reference is the growth rate of the linearized
incompressible Navier-Stokes equations about the same density and viscosity profiles (mu = alpha mu_H + (1 - alpha) mu_L),
a quadratic eigenvalue problem in the vertical velocity w(y), solved by finite differences; the sharp-interface inviscid
value sqrt(A g k) is printed for scale. Exits 1 if |n / n_linear - 1| exceeds --tol.
"""

import argparse
import json
import math
import subprocess
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.optimize import curve_fit

from mfc.viz import assemble_silo, discover_timesteps

LAM, H_INT, KA0, DELTA, RHO_H, RHO_L = 0.2, 0.7, 0.05, 2.5e-3, 3.0, 1.0  # hardcoded IC 204 (2dHardcodedIC.fpp)


def linear_growth(k, g, mu_h, mu_l, half=0.3, n_pts=3200):
    """Largest real n of n^2 M w + n C w + G w = 0, with w = w' = 0 at y = +-half about the interface, where
    M w = k^2 rho w - (rho w')',  C w = (mu (w'' + k^2 w))'' + k^2 mu (w'' + k^2 w) - 4 k^2 (mu w')',  G w = -g k^2 rho' w."""
    h = 2.0 * half / n_pts
    y = -half + h * np.arange(-1, n_pts + 2)  # nodes -1 .. n_pts + 1; the end nodes are ghosts
    alpha = lambda s: 0.5 * (1.0 + np.tanh(s / DELTA))  # noqa: E731
    rho = lambda s: alpha(s) * RHO_H + (1.0 - alpha(s)) * RHO_L  # noqa: E731
    mu = lambda s: alpha(s) * mu_h + (1.0 - alpha(s)) * mu_l  # noqa: E731
    ni, ne = n_pts - 1, n_pts + 3
    ext = sp.lil_matrix((ne, ni))  # interior unknowns -> all nodes: w = 0 on the ends, even ghosts (w' = 0)
    for j in range(ni):
        ext[j + 2, j] = 1.0
    ext[0, 0] = ext[ne - 1, ni - 1] = 1.0
    ext = ext.tocsr()

    def div_flux(c):  # (c w')' with c at the half nodes
        ch = c(0.5 * (y[:-1] + y[1:]))
        main = np.zeros(ne)
        main[:-1] -= ch
        main[1:] -= ch
        return sp.diags([ch, main, ch], [-1, 0, 1], format="csr") / h**2

    d2 = sp.diags([np.ones(ne - 1), -2.0 * np.ones(ne), np.ones(ne - 1)], [-1, 0, 1], format="csr") / h**2
    eye = sp.identity(ne, format="csr")
    q = sp.diags(mu(y), format="csr") @ (d2 + k**2 * eye)
    rows = slice(2, n_pts + 1)
    M = (k**2 * sp.diags(rho(y), format="csr") - div_flux(rho))[rows] @ ext
    C = (d2 @ q + k**2 * q - 4.0 * k**2 * div_flux(mu))[rows] @ ext
    G = (-g * k**2 * sp.diags((RHO_H - RHO_L) * 0.5 / DELTA * (1.0 - np.tanh(y / DELTA) ** 2), format="csr"))[rows] @ ext
    Z, I = sp.csr_matrix((ni, ni)), sp.identity(ni, format="csr")
    A, B = sp.bmat([[Z, I], [-G, -C]], format="csc"), sp.bmat([[I, Z], [Z, M]], format="csc")
    shift = math.sqrt((RHO_H - RHO_L) / (RHO_H + RHO_L) * g * k)
    lu = sla.splu((A - shift * B).tocsc())
    op = sla.LinearOperator(A.shape, matvec=lambda v: lu.solve(B @ v), dtype=float)
    n = shift + 1.0 / sla.eigs(op, k=6, which="LM", return_eigenvectors=False)
    return float(n[np.abs(n.imag) < 1e-8 * np.abs(n)].real.max())


def amplitude(d, step, k):
    data = assemble_silo(d, step)
    x, y = (np.asarray(c).ravel() for c in (data.x_cc, data.y_cc))
    a = np.asarray(data.variables["alpha1"]).squeeze()  # heavy-fluid volume fraction, (nx, ny)
    eta = ((1.0 - a) * (y[1] - y[0])).sum(axis=1)  # light-fluid column height: the interface height, exact to sub-cell displacements
    return -2.0 * float(np.mean(eta * np.cos(k * x)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--case-args", default="", help="Arguments given to case.py for the run.")
    ap.add_argument("--ka-max", type=float, default=0.25, help="Fit window: saves with k a <= this (default 0.25).")
    ap.add_argument("--tol", type=float, default=0.0, help="Fail if |n / n_linear - 1| exceeds this (0: no check).")
    args = ap.parse_args()

    d = args.case_dir
    case = json.loads(subprocess.run([sys.executable, "case.py", *args.case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)
    k, g = 2.0 * math.pi / LAM, -case["g_y"]
    mu_h, mu_l = 1.0 / case["fluid_pp(1)%Re(1)"], 1.0 / case["fluid_pp(2)%Re(1)"]
    n_lin = linear_growth(k, g, mu_h, mu_l)
    n_inv = math.sqrt((RHO_H - RHO_L) / (RHO_H + RHO_L) * g * k)

    steps = discover_timesteps(d, "silo")
    t = np.array(steps) * case["dt"]
    amp = np.array([amplitude(d, s, k) for s in steps])
    over = k * amp > args.ka_max
    end = int(np.argmax(over)) if over.any() else amp.size  # saves before the first one above the cap
    tw, aw = t[:end], amp[:end]
    if tw.size < 6:
        sys.exit(f"only {tw.size} saves with k a <= {args.ka_max}: save more often or run longer")
    a0 = KA0 / k
    popt, _ = curve_fit(lambda s, A, n, B, m: A * np.exp(n * s) + B * np.exp(-m * s), tw, aw, p0=(0.5 * a0, n_inv, 0.5 * a0, n_inv), maxfev=20000)
    n_fit = popt[1]
    slope = float(np.polyfit(tw[-3:], np.log(aw[-3:]), 1)[0])
    print(f"N = {case['m'] + 1} cells per wavelength, {tw.size} saves with k a <= {args.ka_max} (t <= {tw[-1]:.4f} s), k a(0) = {k * amp[0]:.4f}")
    print(f"linear theory: n = {n_lin:.4f} 1/s (viscous, diffuse interface); sqrt(A g k) = {n_inv:.4f} 1/s (inviscid, sharp)")
    print(f"MFC: n = {n_fit:.4f} 1/s ({n_fit / n_lin - 1:+.2%}) from a = A e^(n t) + B e^(-m t); d ln a/dt at the window end {slope:.4f} ({slope / n_lin - 1:+.2%})")
    ok = args.tol <= 0 or abs(n_fit / n_lin - 1.0) <= args.tol
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
