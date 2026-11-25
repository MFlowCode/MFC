#!/usr/bin/env python3
"""1D radial compressible Taylor-Couette reference for examples/2D_ibm_taylor_couette (numpy, shares no code with MFC).

Axisymmetric compressible Navier-Stokes in r on [R, r2] = [1, 2]: rho, rho u_r, rho u_theta, E; ideal gas gamma 1.4; constant mu, no
bulk viscosity and no heat conduction (MFC's model without chemistry); no-slip walls, inner wall turning at Omega from t = 0, outer
wall fixed; gas at rest at t = 0 with rho = 1, c = 1. Finite volumes, central second-order fluxes with a weak fourth-order filter on
rho, rho u_r and E, SSP-RK3. Prints the inner-wall shear S = r d(u_theta/r)/dr and the u_theta error against the incompressible
profile A r + B/r at the end time (the time of MFC's --ntau run).

    python3 ref1d.py --mach 0.05 --re 10 --ntau 12 --N 200

The shear printed is a one-sided cubic through the wall value and three cells, so it carries an O(1/N) error of its own (3.8e-3,
1.9e-3, 9.0e-4 at N = 100, 200, 400); the compressible change is the difference between two Mach numbers at the same N.
"""

import argparse
import math

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--mach", type=float, default=0.05)
ap.add_argument("--re", type=float, default=10.0)
ap.add_argument("--ntau", type=float, default=12.0)
ap.add_argument("--N", type=int, default=200)
ap.add_argument("--cfl", type=float, default=0.4)
ap.add_argument("--eps4", type=float, default=0.02)
ap.add_argument("--nprint", type=int, default=12)
a = ap.parse_args()

gam, R, r2 = 1.4, 1.0, 2.0
om = a.mach / R
nu = a.mach * R / a.re
mu = nu  # rho0 = 1
tau = ((r2 - R) / math.pi) ** 2 / nu
tend = a.ntau * tau
N = a.N
dr = (r2 - R) / N
rf = R + dr * np.arange(N + 1)
rc = 0.5 * (rf[1:] + rf[:-1])
vol = 0.5 * (rf[1:] ** 2 - rf[:-1] ** 2)  # per radian
den = r2**2 - R**2
A, B = -om * R**2 / den, om * R**2 * r2**2 / den
v_ex = A * rc + B / rc
S_ex = -2.0 * B / R**2

U = np.zeros((4, N))
U[0] = 1.0
U[3] = (1.0 / gam) / (gam - 1.0)


def prim(U):
    rho = U[0]
    u = U[1] / rho
    v = U[2] / rho
    p = (gam - 1.0) * (U[3] - 0.5 * rho * (u * u + v * v))
    return rho, u, v, p


def ghost(q, lo, hi):
    """Extend a cell array with one ghost cell each side; lo/hi are (kind, wall value): 'odd' mirror about the wall value or 'even'."""
    out = np.empty(N + 2)
    out[1:-1] = q
    for idx, nb, (kind, w) in ((0, 1, lo), (-1, -2, hi)):
        out[idx] = 2.0 * w - out[nb] if kind == "odd" else out[nb]
    return out


def rhs(U):
    rho, u, v, p = prim(U)
    E = U[3]
    rg = ghost(rho, ("even", 0), ("even", 0))
    ug = ghost(u, ("odd", 0.0), ("odd", 0.0))
    vg = ghost(v, ("odd", om * R), ("odd", 0.0))
    pg = ghost(p, ("even", 0), ("even", 0))
    Eg = ghost(E, ("even", 0), ("even", 0))
    # face values (central)
    rF = 0.5 * (rg[1:] + rg[:-1])
    uF = 0.5 * (ug[1:] + ug[:-1])
    vF = 0.5 * (vg[1:] + vg[:-1])
    pF = 0.5 * (pg[1:] + pg[:-1])
    EF = 0.5 * (Eg[1:] + Eg[:-1])
    uF[0] = uF[-1] = 0.0
    vF[0], vF[-1] = om * R, 0.0
    # face gradients
    du = (ug[1:] - ug[:-1]) / dr
    rcg = np.r_[rc[0] - dr, rc, rc[-1] + dr]
    dw = (vg[1:] / rcg[1:] - vg[:-1] / rcg[:-1]) / dr  # d(v/r)/dr at faces
    # divergence at faces: du/dr + u/r
    divF = du + uF / rf
    trr = mu * (2.0 * du - (2.0 / 3.0) * divF)
    trt = mu * rf * dw
    # cell-centred quantities for the hoop terms
    divc = (rf[1:] * uF[1:] - rf[:-1] * uF[:-1]) / vol
    ttt = mu * (2.0 * u / rc - (2.0 / 3.0) * divc)
    # fluxes times face radius (or r^2 for the angular momentum)
    Fm = rf * rF * uF
    Fu = rf * (rF * uF * uF - trr)
    Fv = rf**2 * (rF * uF * vF - trt)
    Fe = rf * ((EF + pF) * uF - uF * trr - vF * trt)
    out = np.empty_like(U)
    out[0] = -(Fm[1:] - Fm[:-1]) / vol
    # pressure gradient in non-conservative form: -dp/dr from faces, plus hoop (rho v^2 - t_tt)/r
    out[1] = -(Fu[1:] - Fu[:-1]) / vol - (pF[1:] - pF[:-1]) / dr + (rho * v * v - ttt) / rc
    out[2] = -(Fv[1:] - Fv[:-1]) / (rc * vol)
    out[3] = -(Fe[1:] - Fe[:-1]) / vol
    # fourth-order filter on rho, rho u, E (acoustics only; rho v untouched)
    if a.eps4 > 0:
        c = np.sqrt(gam * p / rho)
        for k in (0, 1, 3):
            q = U[k]
            qe = np.r_[q[1], q[0], q, q[-1], q[-2]] if k != 1 else np.r_[-q[1], -q[0], q, -q[-1], -q[-2]]
            d4 = qe[4:] - 4 * qe[3:-1] + 6 * qe[2:-2] - 4 * qe[1:-3] + qe[:-4]
            out[k] -= a.eps4 * c / dr * d4
    return out


cmax = 1.0 + a.mach
dt_a = a.cfl * dr / cmax
dt_v = 0.2 * dr**2 / (4.0 / 3.0 * nu)
dt = min(dt_a, dt_v)
nsteps = int(math.ceil(tend / dt))
dt = tend / nsteps


def report(U, t):
    rho, u, v, p = prim(U)
    # inner-wall shear from a one-sided cubic through the wall value and three cells
    w = np.r_[om, v[:3] / rc[:3]]
    rr = np.r_[R, rc[:3]] - R
    cfit = np.polyfit(rr, w, 3)
    S = R * cfit[-2]
    err = (v - v_ex) / (om * R)
    T = p / rho * gam  # T/T0 with T0 = p0/rho0 * gam = 1
    return dict(
        t=t,
        S_err=S / S_ex - 1.0,
        L2=math.sqrt(np.mean(err**2)),
        maxerr=np.abs(err).max(),
        umax=np.abs(u).max() / (om * R),
        rho_in=rho[0],
        rho_out=rho[-1],
        T_in=T[0],
        T_out=T[-1],
        mass=float(np.sum(rho * vol)) / float(np.sum(vol)),
    )


print(f"M {a.mach} Re {a.re} N {N} tau {tau:.3f} t_end {tend:.2f} steps {nsteps} dt {dt:.3e}")
every = max(nsteps // a.nprint, 1)
t = 0.0
for n in range(nsteps):
    k1 = U + dt * rhs(U)
    k2 = 0.75 * U + 0.25 * (k1 + dt * rhs(k1))
    U = (U + 2.0 * (k2 + dt * rhs(k2))) / 3.0
    t += dt
    if (n + 1) % every == 0 or n + 1 == nsteps:
        r = report(U, t)
        print(
            "t/tau {:6.2f}  S_err {:+.3e}  L2 {:.3e}  max {:.3e}  |u_r|max/OmR {:.2e}  T_in-1 {:.2e} T_out-1 {:.2e}  rho_in {:.5f} rho_out {:.5f}".format(
                t / tau, r["S_err"], r["L2"], r["maxerr"], r["umax"], r["T_in"] - 1, r["T_out"] - 1, r["rho_in"], r["rho_out"]
            )
        )
