#!/usr/bin/env python3
"""
Body-fitted reference for the sphere bow-shock stand-off (an independent solver, standalone NumPy): inviscid axisymmetric
Euler on a polar mesh in the meridional plane, independent of the immersed boundary.

Mesh: rho in [R, R_out] (geometric stretching from the wall), polar angle phi from the upstream axis (phi = 0, the stagnation
line) to phi_max past the shoulder. Finite volumes weighted by the cylindrical radius y (exact for axisymmetry): cell
d/dt (U y_c A) = - sum_faces F_n y_f L_f + [0, 0, p, 0] A, with Cartesian (x, y) velocity components, MUSCL (minmod) on
primitive variables along each grid direction, HLLC on the rotated face normal, SSP-RK3, global CFL step. Boundaries: wall
(reflective) at rho = R, free stream at rho = R_out (upstream of the shock), symmetry on the axis (phi = 0), zero gradient at
phi_max (supersonic outflow). Stand-off on the stagnation line (first phi row) with the same definition as the IB runs.

    python3 bodyfitted_sphere.py [--mach 2] [--nr 160] [--nphi 160] [--rout 2.2] [--tend 40]
"""

import argparse
import json

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--mach", type=float, default=2.0)
ap.add_argument("--nr", type=int, default=160)
ap.add_argument("--nphi", type=int, default=160)
ap.add_argument("--rout", type=float, default=2.2, help="Outer radius / R (must enclose the shock on the axis).")
ap.add_argument("--phimax", type=float, default=110.0, help="Degrees from the stagnation line.")
ap.add_argument("--stretch", type=float, default=1.0, help="Radial cell growth ratio.")
ap.add_argument("--tend", type=float, default=40.0, help="In R / u_inf.")
ap.add_argument("--cfl", type=float, default=0.4)
a = ap.parse_args()

G, R, M = 1.4, 1.0, a.mach
uinf = M  # c_inf = 1
T = lambda x: np.asarray(x, dtype=float)  # noqa: E731

# Mesh (node coordinates)
h = a.stretch ** np.arange(a.nr)
rn = R + np.concatenate([[0.0], np.cumsum(h)]) * (a.rout - R) / h.sum()
phn = np.linspace(0.0, np.radians(a.phimax), a.nphi + 1)
RN, PHN = np.meshgrid(rn, phn, indexing="ij")
# upstream axis is -x: x = -rho cos(phi), y = rho sin(phi) (sphere centred at the origin)
XN, YN = -RN * np.cos(PHN), RN * np.sin(PHN)


def quad_area_centroid(x, y):
    """Area and centroid of each quad cell (counter-clockwise nodes)."""
    x0, x1, x2, x3 = x[:-1, :-1], x[1:, :-1], x[1:, 1:], x[:-1, 1:]
    y0, y1, y2, y3 = y[:-1, :-1], y[1:, :-1], y[1:, 1:], y[:-1, 1:]
    A = 0.5 * ((x0 * y1 - x1 * y0) + (x1 * y2 - x2 * y1) + (x2 * y3 - x3 * y2) + (x3 * y0 - x0 * y3))
    cx = ((x0 + x1) * (x0 * y1 - x1 * y0) + (x1 + x2) * (x1 * y2 - x2 * y1) + (x2 + x3) * (x2 * y3 - x3 * y2) + (x3 + x0) * (x3 * y0 - x0 * y3)) / (6 * A)
    cy = ((y0 + y1) * (x0 * y1 - x1 * y0) + (y1 + y2) * (x1 * y2 - x2 * y1) + (y2 + y3) * (x2 * y3 - x3 * y2) + (y3 + y0) * (x3 * y0 - x0 * y3)) / (6 * A)
    return np.abs(A), cx, cy


A, XC, YC = quad_area_centroid(XN, YN)
# radial faces (constant rho, between phi nodes): i face between cells i-1 and i, i = 0..nr
fr_dx, fr_dy = XN[:, 1:] - XN[:, :-1], YN[:, 1:] - YN[:, :-1]
fr_L = np.hypot(fr_dx, fr_dy)
fr_nx, fr_ny = fr_dy / fr_L, -fr_dx / fr_L  # normal; orientation fixed below to point to +rho
s = np.sign(fr_nx * (-np.cos(0.5 * (PHN[:, 1:] + PHN[:, :-1]))) + fr_ny * np.sin(0.5 * (PHN[:, 1:] + PHN[:, :-1])))
fr_nx, fr_ny = fr_nx * s, fr_ny * s
fr_y = 0.5 * (YN[:, 1:] + YN[:, :-1])
# angular faces (constant phi, between rho nodes): j face between cells j-1 and j, j = 0..nphi
fp_dx, fp_dy = XN[1:, :] - XN[:-1, :], YN[1:, :] - YN[:-1, :]
fp_L = np.hypot(fp_dx, fp_dy)
fp_nx, fp_ny = -fp_dy / fp_L, fp_dx / fp_L
pm = PHN[0, :][None, :]
s = np.sign(fp_nx * np.sin(pm) + fp_ny * np.cos(pm))  # point to +phi: direction d/dphi of (-cos, sin) = (sin, cos)
fp_nx, fp_ny = fp_nx * s, fp_ny * s
fp_y = 0.5 * (YN[1:, :] + YN[:-1, :])

A_, YC_ = T(A), T(YC)
frN = (T(fr_nx), T(fr_ny), T(fr_L * fr_y))
fpN = (T(fp_nx), T(fp_ny), T(fp_L * fp_y))
nr, nphi = a.nr, a.nphi


def prim(U):
    r = U[0]
    u = U[1] / r
    v = U[2] / r
    return np.stack([r, u, v, (G - 1) * (U[3] - 0.5 * r * (u * u + v * v))])


def cons(W):
    r, u, v, p = W
    return np.stack([r, r * u, r * v, p / (G - 1) + 0.5 * r * (u * u + v * v)])


def minmod(a_, b_):
    return np.where(a_ * b_ > 0, np.sign(a_) * np.minimum(np.abs(a_), np.abs(b_)), np.zeros_like(a_))


def hllc_n(WL, WR, nx, ny):
    """HLLC flux through a face of unit normal (nx, ny), returned in Cartesian components."""
    rL, uL, vL, pL = WL
    rR, uR, vR, pR = WR
    unL, unR = uL * nx + vL * ny, uR * nx + vR * ny
    EL = pL / (G - 1) + 0.5 * rL * (uL * uL + vL * vL)
    ER = pR / (G - 1) + 0.5 * rR * (uR * uR + vR * vR)
    cL, cR = np.sqrt(G * np.abs(pL) / rL), np.sqrt(G * np.abs(pR) / rR)
    SL, SR = np.minimum(unL - cL, unR - cR), np.maximum(unL + cL, unR + cR)
    SM = (pR - pL + rL * unL * (SL - unL) - rR * unR * (SR - unR)) / (rL * (SL - unL) - rR * (SR - unR))

    def F(r, u, v, p, E, un):
        return np.stack([r * un, r * u * un + p * nx, r * v * un + p * ny, un * (E + p)])

    def Fstar(r, u, v, p, E, un, S):
        f = r * (S - un) / (S - SM)
        Us = np.stack([f, f * (u + (SM - un) * nx), f * (v + (SM - un) * ny), f * (E / r + (SM - un) * (SM + p / (r * (S - un))))])
        Uk = np.stack([r, r * u, r * v, E])
        return F(r, u, v, p, E, un) + S * (Us - Uk)

    FL, FR = F(rL, uL, vL, pL, EL, unL), F(rR, uR, vR, pR, ER, unR)
    FLs, FRs = Fstar(rL, uL, vL, pL, EL, unL, SL), Fstar(rR, uR, vR, pR, ER, unR, SR)
    return np.where(SL >= 0, FL, np.where(SM >= 0, FLs, np.where(SR > 0, FRs, FR)))


Winf = T([1.0, uinf, 0.0, 1 / G])


def ghosts(W):
    """Pad W [4, nr, nphi] by 2 cells on every side: wall (rho lo), free stream (rho hi), axis (phi lo), outflow (phi hi)."""

    def reflect(Wg, nx, ny):  # mirror the normal velocity about a face of normal (nx, ny)
        un = Wg[1] * nx + Wg[2] * ny
        return np.stack([Wg[0], Wg[1] - 2 * un * nx, Wg[2] - 2 * un * ny, Wg[3]])

    # rho direction
    lo = np.flip(W[:, :2, :], 1)
    nxw, nyw = frN[0][0][None, :], frN[1][0][None, :]
    lo = reflect(lo, nxw, nyw)
    hi = np.broadcast_to(Winf[:, None, None], (4, 2, W.shape[2]))
    W = np.concatenate([lo, W, hi], 1)
    # phi direction (after padding rho, so corners get values)
    alo = np.flip(W[:, :, :2], 2)
    alo = np.stack([alo[0], alo[1], -alo[2], alo[3]])  # axis: y-velocity odd
    ahi = np.broadcast_to(W[:, :, -1:], (4, W.shape[1], 2))
    return np.concatenate([alo, W, ahi], 2)


def rhs(U):
    W = prim(U)
    Wp = ghosts(W)  # [4, nr+4, nphi+4]
    # rho-direction faces i = 0..nr (between padded cells i+1 and i+2)
    Wc = Wp[:, :, 2:-2]
    dm, dp_ = Wc[:, 1:-1] - Wc[:, :-2], Wc[:, 2:] - Wc[:, 1:-1]
    sl = minmod(dm, dp_)  # slopes of padded cells 1..nr+2
    WL = Wc[:, 1:-2] + 0.5 * sl[:, :-1]  # left state at faces 0..nr: cell i+1 (padded) right face
    WR = Wc[:, 2:-1] - 0.5 * sl[:, 1:]
    Fr = hllc_n(WL, WR, frN[0], frN[1]) * frN[2]
    # wall face: exact pressure flux only
    pw = 0.5 * (WL[3][0] + WR[3][0])
    Fr[:, 0] = np.stack([np.zeros_like(pw), pw * frN[0][0], pw * frN[1][0], np.zeros_like(pw)]) * frN[2][0]
    # phi-direction faces j = 0..nphi
    Wc = Wp[:, 2:-2, :]
    dm, dp_ = Wc[:, :, 1:-1] - Wc[:, :, :-2], Wc[:, :, 2:] - Wc[:, :, 1:-1]
    sl = minmod(dm, dp_)
    WL = Wc[:, :, 1:-2] + 0.5 * sl[:, :, :-1]
    WR = Wc[:, :, 2:-1] - 0.5 * sl[:, :, 1:]
    Fp = hllc_n(WL, WR, fpN[0], fpN[1]) * fpN[2]
    div = (Fr[:, 1:, :] - Fr[:, :-1, :]) + (Fp[:, :, 1:] - Fp[:, :, :-1])
    S = np.zeros_like(U)
    S[2] = W[3] * A_
    return (S - div) / (YC_ * A_)


U = cons(np.broadcast_to(Winf[:, None, None], (4, nr, nphi)).copy())
dmin = float(min(np.min(fr_L), np.min(fp_L[:, 1:])))
t, t_end = 0.0, a.tend * R / uinf
p2 = (1 / G) * (1 + 2 * G / (G + 1) * (M**2 - 1))
pmid = 0.5 * (1 / G + p2)
rc = 0.5 * (rn[1:] + rn[:-1])


def standoff(U):
    p = prim(U)[3][:, 0]  # first phi row (next to the axis), from the wall outwards
    for i in range(len(p) - 1):
        if p[i] >= pmid > p[i + 1]:
            return float(rc[i] + (pmid - p[i]) / (p[i + 1] - p[i]) * (rc[i + 1] - rc[i]) - R)
    return float("nan")


hist, nxt = [], 0.0
while t < t_end - 1e-12:
    W = prim(U)
    c = np.sqrt(G * np.abs(W[3]) / W[0])
    smax = float((np.hypot(W[1], W[2]) + c).max())
    dt = min(a.cfl * dmin / smax, t_end - t)
    k1 = rhs(U)
    U1 = U + dt * k1
    U2 = 0.75 * U + 0.25 * (U1 + dt * rhs(U1))
    U = U / 3 + 2 / 3 * (U2 + dt * rhs(U2))
    t += dt
    if t >= nxt:
        hist.append((round(t * uinf / R, 1), round(standoff(U), 5)))
        nxt += 2.0 * R / uinf
d = standoff(U)
billig = 0.143 * np.exp(3.24 / M**2)
print(json.dumps({"mach": M, "nr": nr, "nphi": nphi, "rout": a.rout, "standoff": d, "billig": billig, "err_vs_billig": d / billig - 1, "history": hist[-5:]}))
