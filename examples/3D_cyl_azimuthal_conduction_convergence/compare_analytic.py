#!/usr/bin/env python3
"""
Verification for examples/3D_cyl_azimuthal_conduction_convergence.

At t = 0 the state is quiescent (u = 0) at uniform pressure, so every Euler flux and every
cylindrical geometric source vanishes and the only RHS contribution is Fourier conduction:

    (rho*E(dt) - rho*E(0)) / dt  =  k * grad^2 T + O(dt)
                                 =  3*k*A*T0*cos(theta)/R**2,

a constant times cos(theta), the same at every radius. The radial half of that answer,
4*k*A*T0*cos(theta)/R**2, is reproduced exactly by the discretization (T is quadratic in r on a
uniform r-grid), so subtracting it isolates the azimuthal half, -k*A*T0*cos(theta)/R**2, which is
what the (1/r**2) metric in m_conduction has to produce. Dropping either metric factor multiplies
that half by r**2 -- an r-dependent error that no refinement removes -- so the table reports the
azimuthal half ring by ring as well as in aggregate.

Cells excluded from the bulk: k = 0 reads the ghost across the axis (the documented non-converging
axis error), and the two outermost cells. k = NR-1 reads the mirrored wall ghost, where dT/dr is not
actually zero, giving it an O(1/dr) RHS; k = NR-2 then picks up an O(dt) share of that through the
later Runge-Kutta stages.

    for np_ in 32 64 128 256; do
        NP=$np_ ./mfc.sh run examples/3D_cyl_azimuthal_conduction_convergence/case.py -n 1
        NP=$np_ ./build/venv/bin/python3 examples/3D_cyl_azimuthal_conduction_convergence/compare_analytic.py
    done
"""

import os
import sys

import numpy as np

GAM = 1.4
CV = 1.0
RHO0 = 1.0
P0 = 1.0
AMP = float(os.environ.get("AMP", "0.4"))
K_THERM = float(os.environ.get("K_THERM", "1.0"))
R = 2.0
DT = 1.0e-8

NX = int(os.environ.get("NX", "32"))
NR = int(os.environ.get("NR", "32"))
NP = int(os.environ.get("NP", "32"))

HERE = os.path.dirname(os.path.abspath(__file__))
RESTART = os.path.join(HERE, "restart_data")

NVAR = 6  # alpha_rho(1), mom_x, mom_y, mom_z, E, alpha(1)
E_IDX = 4  # zero-based index of rho*E in the record


def read_step(step):
    path = os.path.join(RESTART, f"lustre_{step}.dat")
    if not os.path.exists(path):
        sys.exit(f"missing {path}; run the case first")
    # Fortran order (i, j, k) = (x, r, theta) per variable, x fastest
    return np.fromfile(path, dtype=np.float64).reshape((NVAR, NP, NR, NX))


def centers(name, n):
    cb = np.fromfile(os.path.join(RESTART, f"lustre_{name}_cb.dat"), dtype=np.float64)
    if cb.size != n + 1:
        sys.exit(f"{name}_cb has {cb.size} entries, expected {n + 1}; check NR/NP")
    return 0.5 * (cb[1:] + cb[:-1])


def main():
    rhs = ((read_step(1)[E_IDX] - read_step(0)[E_IDX]) / DT).mean(axis=2)  # average over x
    r = centers("y", NR)
    theta = centers("z", NP)

    t0 = P0 / ((GAM - 1.0) * RHO0 * CV)
    cos_t = np.cos(theta)[:, None]
    radial = 4.0 * K_THERM * AMP * t0 * cos_t / R**2  # exact, and exact in the discretization
    azimuthal = -1.0 * K_THERM * AMP * t0 * cos_t / R**2
    exact = radial + azimuthal
    scale = np.sqrt(np.mean(exact**2))

    bulk = slice(1, NR - 2)
    rel = (rhs - exact)[:, bulk] / scale
    dtheta = 2.0 * np.pi / NP

    # Azimuthal half alone, recovered ring by ring: amplitude of (rhs - radial) projected on cos.
    got = ((rhs - radial) * cos_t).sum(axis=0) / (cos_t**2).sum()
    want = -K_THERM * AMP * t0 / R**2

    print(f"NP={NP:5d} NR={NR:4d}  bulk L2={np.sqrt(np.mean(rel**2)):.3e} rel  max={np.max(np.abs(rel)):.3e} rel  " f"(predicted dtheta^2/36 = {dtheta**2 / 36:.3e})")
    print("            azimuthal half / exact:  " + "  ".join(f"r={r[k]:.3f}: {got[k] / want:+.6f}" for k in (1, NR // 4, NR // 2, NR - 3)))


if __name__ == "__main__":
    main()
