#!/usr/bin/env python3
"""
Flame + vortex-ring initial condition for 3D_flame_vortex_ring: overwrites RUN/restart_data/lustre_0.dat written by pre_process.

    ic3d.py RUN [--gamma G]

Flame: the Cantera profile along x in the lab frame (reactants at rest), T = 1200 K at x = 1.8 mm, p = 1 atm.
Ring: Gaussian-core vorticity w = G/(pi delta^2) exp(-s^2/delta^2) around the circle (x = 0.9 mm, radius 0.4 mm) about the y-z centre
line; the velocity comes from its Stokes stream function (thin-ring kernel with complete elliptic integrals, summed over the core),
tabulated on a fine (x, r) mesh and interpolated to the cells. The core's pressure deficit (about 0.2 % of p) is left out and relaxes
acoustically, identically on every grid. Layout of lustre_0.dat: float64 (var, z, y, x), vars rho, rho u, rho v, rho w, E, alpha,
rho Y_k, with E = rho e(T, Y) + rho |u|^2 / 2 (Cantera, the same NASA7 data as MFC's thermochemistry).
"""

import argparse
import importlib.util
import os
import re

import cantera as ct
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.special import ellipe, ellipk

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("fv_common", os.path.join(HERE, "..", "2D_flame_vortex_pair", "fv_common.py"))
fv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fv)

L, P0, MECH = fv.L, fv.P0, fv.MECH
XV, RR, DELTA, GAMMA = 0.9e-3, 0.4e-3, 0.15e-3, 0.022  # ring centre x, radius, Gaussian core radius, circulation [m^2/s]


def inp(run, key):
    return re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", open(os.path.join(run, "pre_process.inp")).read(), re.M).group(1)


def ring_table(gamma, h=5e-6):
    """u_x, u_r on an (x, r) mesh of spacing h covering the domain (r up to the cube corner)."""
    xs = np.arange(-2, int(L / h) + 3) * h
    rs = (np.arange(int(0.75 * L / h) + 3) + 0.5) * h
    ns, nphi = 48, 96  # midpoint rule on a polar patch around (XV, RR), out to 3.5 delta
    s = (np.arange(ns) + 0.5) * 3.5 * DELTA / ns
    ph = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
    S, PH = np.meshgrid(s, ph, indexing="ij")
    x0, r0 = (XV + S * np.cos(PH)).ravel(), (RR + S * np.sin(PH)).ravel()
    w0 = (gamma / (np.pi * DELTA**2) * np.exp(-(S**2) / DELTA**2) * S * (3.5 * DELTA / ns) * (2 * np.pi / nphi)).ravel()
    keep = r0 > 0  # the Gaussian tail across the axis (< 1e-3 of the circulation) is dropped
    x0, r0, w0 = x0[keep], r0[keep], w0[keep]
    X, R = np.meshgrid(xs, rs, indexing="ij")
    psi = np.zeros_like(X)
    for i in range(0, len(x0), 64):
        a, b, c = x0[i : i + 64, None, None], r0[i : i + 64, None, None], w0[i : i + 64, None, None]
        m = np.clip(4 * R * b / ((X - a) ** 2 + (R + b) ** 2), 0, 1 - 1e-15)
        k = np.sqrt(m)
        psi += (c * np.sqrt(R * b) / (2 * np.pi) * ((2 / k - k) * ellipk(m) - 2 / k * ellipe(m))).sum(0)
    # psi is even in r (psi ~ r^2 at the axis): a mirrored row keeps d/dr central at the first point
    ux = np.gradient(np.c_[psi[:, :1], psi], np.r_[-rs[:1], rs], axis=1)[:, 1:] / R
    ur = -np.gradient(psi, xs, axis=0) / R
    return xs, rs, ux, ur


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--gamma", type=float, default=GAMMA)
    a = ap.parse_args()
    m, n, p = (int(inp(a.run, k)) + 1 for k in "mnp")
    dx = L / m
    xc, yc, zc = (np.arange(m) + 0.5) * dx, (np.arange(n) + 0.5) * dx, (np.arange(p) + 0.5) * dx
    gas = ct.Solution(MECH)
    ns = gas.n_species
    T1, u1, Y1 = fv.flame_on(xc, fv.profile())
    sa = ct.SolutionArray(gas, m)
    sa.TPY = T1, P0, Y1.T
    rho1, e1 = sa.density, sa.int_energy_mass
    xs, rs, uxt, urt = ring_table(a.gamma) if a.gamma else (np.array([0, L]), np.array([0, L]), np.zeros((2, 2)), np.zeros((2, 2)))
    Ix = RegularGridInterpolator((xs, rs), uxt, bounds_error=False, fill_value=None)
    Ir = RegularGridInterpolator((xs, rs), urt, bounds_error=False, fill_value=None)
    yy, xx = np.meshgrid(yc - 0.5 * n * dx, xc, indexing="ij")  # one z-slab at a time: (y, x)
    out = os.path.join(a.run, "restart_data", "lustre_0.dat")
    q = np.memmap(out, "<f8", "r+").reshape(-1, p, n, m)
    assert q.shape[0] == 6 + ns, q.shape
    for kz, z in enumerate(zc - 0.5 * p * dx):
        r = np.hypot(yy, z)
        pts = np.c_[xx.ravel(), r.ravel()]
        ux, ur = Ix(pts).reshape(r.shape), Ir(pts).reshape(r.shape)
        rs_ = np.maximum(r, 1e-30)
        uy, uz = ur * yy / rs_, ur * z / rs_
        ux = ux + u1[None, :]
        rho = np.broadcast_to(rho1, r.shape)
        q[0, kz] = rho
        q[1, kz], q[2, kz], q[3, kz] = rho * ux, rho * uy, rho * uz
        q[4, kz] = rho * e1 + 0.5 * rho * (ux**2 + uy**2 + uz**2)
        q[5, kz] = 1.0
        for k in range(ns):
            q[6 + k, kz] = rho * Y1[k]
    q.flush()
    print(f"{out}: {m}x{n}x{p}, ring {'off' if not a.gamma else f'gamma {a.gamma}'}")


if __name__ == "__main__":
    main()
