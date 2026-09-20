#!/usr/bin/env python3
"""
Verification for examples/1D_conduction_convergence.

At t = 0 the state is quiescent (u = 0) at uniform pressure, so every Euler flux is
zero and the only RHS contribution is Fourier conduction. One time step therefore gives

    (rho*E(dt) - rho*E(0)) / dt  =  k * d2T/dx2 + O(dt) + O(dx^2)

with the exact right-hand side

    d2T/dx2 = -T0 * A * (2*pi/L)**2 * sin(2*pi*x/L).

Run at several NX and confirm the L2 error falls at second order.

    for nx in 50 100 200 400; do
        NX=$nx ./mfc.sh run examples/1D_conduction_convergence/case.py -n 1
        NX=$nx ./build/venv/bin/python3 examples/1D_conduction_convergence/compare_analytic.py
    done
"""

import os
import sys

import numpy as np

GAM = 1.4
CV = 1.0
RHO0 = 1.0
P0 = 1.0
AMP = 0.1
K_THERM = 1.0e-3
L = 1.0
DT = 1.0e-8

NX = int(os.environ.get("NX", "100"))
HERE = os.path.dirname(os.path.abspath(__file__))
RESTART = os.path.join(HERE, "restart_data")

NVAR = 4  # alpha_rho(1), mom_x, E, alpha(1)
E_IDX = 2  # zero-based index of rho*E in the record


def read_step(step):
    path = os.path.join(RESTART, f"lustre_{step}.dat")
    if not os.path.exists(path):
        sys.exit(f"missing {path}; run the case first")
    raw = np.fromfile(path, dtype=np.float64)
    return raw.reshape((NVAR, NX))


def main():
    q0 = read_step(0)
    q1 = read_step(1)

    dx = L / NX
    x = (np.arange(NX) + 0.5) * dx

    t0 = P0 / ((GAM - 1.0) * RHO0 * CV)
    d2t_exact = -t0 * AMP * (2.0 * np.pi / L) ** 2 * np.sin(2.0 * np.pi * x / L)
    rhs_exact = K_THERM * d2t_exact

    rhs_num = (q1[E_IDX] - q0[E_IDX]) / DT

    err = np.sqrt(np.mean((rhs_num - rhs_exact) ** 2))
    scale = np.sqrt(np.mean(rhs_exact**2))
    print(f"NX={NX:5d}  L2 error={err:.6e}  relative={err / scale:.6e}")


if __name__ == "__main__":
    main()
