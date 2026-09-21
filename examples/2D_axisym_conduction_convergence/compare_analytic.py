#!/usr/bin/env python3
"""
Verification for examples/2D_axisym_conduction_convergence.

At t = 0 the state is quiescent (u = 0) at uniform pressure, so every Euler flux and every
cylindrical geometric source vanishes and the only RHS contribution is Fourier conduction:

    (rho*E(dt) - rho*E(0)) / dt  =  k * (1/r) d/dr (r dT/dr) + O(dt)
                                 =  -4*k*A*T0/R**2 - 16*k*B*T0*r**2/R**4.

With the default B = 0 that is a constant in every radial cell, the axis cell k = 0 included.
The axis cell has no face pair, so it is the one cell served by s_compute_conduction_axis_source;
measuring it against a flat exact field makes a wrong sign or factor there unmissable.

Two cells are reported apart from the bulk. Cell 0 is the axis cell itself. Cell 1 is the one
cell whose inner face lies on the half-width cell that the axisymmetric grid places at r = 0
(m_grid.f90 gives that cell half the width of the rest); there the two-point face gradient that
every MFC diffusive flux uses is only first order, so its error does not shrink with NR. The
outermost cells read the mirrored wall ghost.

    # constant exact answer: the axis-cell check
    for nr in 25 50 100 200; do
        NR=$nr ./mfc.sh run examples/2D_axisym_conduction_convergence/case.py -n 1
        NR=$nr ./build/venv/bin/python3 examples/2D_axisym_conduction_convergence/compare_analytic.py
    done

    # r-dependent exact answer: the second-order convergence table
    for nr in 25 50 100 200; do
        NR=$nr BAMP=0.4 ./mfc.sh run examples/2D_axisym_conduction_convergence/case.py -n 1
        NR=$nr BAMP=0.4 ./build/venv/bin/python3 examples/2D_axisym_conduction_convergence/compare_analytic.py
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
BAMP = float(os.environ.get("BAMP", "0.0"))
K_THERM = float(os.environ.get("K_THERM", "1.0e-2"))
R = 1.0
NX = 30
DT = 1.0e-8

NR = int(os.environ.get("NR", "100"))
HERE = os.path.dirname(os.path.abspath(__file__))
RESTART = os.path.join(HERE, "restart_data")

NVAR = 5  # alpha_rho(1), mom_x, mom_y, E, alpha(1)
E_IDX = 3  # zero-based index of rho*E in the record


def read_step(step):
    path = os.path.join(RESTART, f"lustre_{step}.dat")
    if not os.path.exists(path):
        sys.exit(f"missing {path}; run the case first")
    raw = np.fromfile(path, dtype=np.float64)
    return raw.reshape((NVAR, NR, NX))


def radii():
    """Cell centers of the axisymmetric grid: m_grid.f90 halves the cell at r = 0."""
    h = R / (2 * NR - 1)
    return np.concatenate(([0.5 * h], 2.0 * h * np.arange(1, NR)))


def main():
    q0 = read_step(0)
    q1 = read_step(1)

    t0 = P0 / ((GAM - 1.0) * RHO0 * CV)
    r = radii()
    exact = -4.0 * K_THERM * AMP * t0 / R**2 - 16.0 * K_THERM * BAMP * t0 * r**2 / R**4

    # T varies in r only, so average over x; the spread across x is zero by construction.
    rhs = ((q1[E_IDX] - q0[E_IDX]) / DT).mean(axis=1)
    scale = np.sqrt(np.mean(exact**2))
    rel = (rhs - exact) / scale

    bulk = rel[2:-2]
    print(f"NR={NR:5d}  axis={rhs[0]:.8e} exact={exact[0]:.8e} ({rel[0]:+.3e} rel)  " f"cell1={rel[1]:+.3e} rel  " f"bulk L2={np.sqrt(np.mean(bulk**2)):.3e} rel  max={np.max(np.abs(bulk)):.3e} rel")


if __name__ == "__main__":
    main()
