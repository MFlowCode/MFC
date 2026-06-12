#!/usr/bin/env python3
"""
Compare the steady MFC velocity profile of the 2D shear-thickening power-law
Poiseuille channel (examples/2D_poiseuille_thickening_nn/case.py) against the
closed-form analytic solution.

For a power-law fluid (Herschel-Bulkley with yield stress tau0 = 0, consistency
index K, flow index n) driven by a constant body acceleration g in a channel of
height L_y (half-height H = L_y/2, centerline at y = H), the steady fully-developed
profile is

    u(y) = (n/(n+1)) * (rho*g/K)^(1/n) * ( H^((n+1)/n) - |y - H|^((n+1)/n) ).

For n > 1 (shear-thickening) the effective viscosity mu = K*|du/dy|^(n-1) -> 0 at
the shear-free centerline, so no regularization cap is needed and the analytic
profile is exact everywhere; the profile is more POINTED than a parabola.

This script reads the raw restart binaries (restart_data/lustre_*.dat), extracts
the x-averaged u(y) from the last save, computes the analytic profile at the same
cell centers, and reports the relative L2 error

    err = sqrt( sum((u_num - u_ana)^2) / sum(u_ana^2) ).

The local momentum balance K*|du/dy|^n = rho*g*(H - y) is also checked. Steady
state is confirmed by comparing the last two saves.

Run:  ./build/venv/bin/python3 examples/2D_poiseuille_thickening_nn/compare_analytic.py
"""

import glob
import os

import numpy as np

# Must match case.py
RHO = 1.0
G_X = 5.0e-2
K = 5.0e-2
NN = 1.5
L_X = 0.2
L_Y = 0.2
H = 0.5 * L_Y

# Grid (m, n in case.py -> m+1, n+1 cells)
NX = 25
NY = 64
# Conservative variables per cell: alpha_rho(1), mom_x, mom_y, E, alpha(1)
NVAR = 5

HERE = os.path.dirname(os.path.abspath(__file__))
RESTART_DIR = os.path.join(HERE, "restart_data")


def u_analytic(y):
    """Closed-form steady power-law Poiseuille profile."""
    expo = (NN + 1.0) / NN
    pref = (NN / (NN + 1.0)) * (RHO * G_X / K) ** (1.0 / NN)
    return pref * (H**expo - np.abs(y - H) ** expo)


def read_profile(path):
    """Return x-averaged u(y) from one restart binary (Fortran order, x fastest)."""
    raw = np.fromfile(path, dtype=np.float64).reshape((NVAR, NY, NX))
    u = raw[1] / raw[0]  # mom_x / alpha_rho -> u[y, x]
    return u.mean(axis=1)  # average over x (flow is x-invariant)


def main():
    dats = sorted(
        glob.glob(os.path.join(RESTART_DIR, "lustre_[0-9]*.dat")),
        key=lambda p: int(os.path.basename(p).split("_")[1].split(".")[0]),
    )
    if not dats:
        raise SystemExit(f"No restart files in {RESTART_DIR}. " "Run: ./mfc.sh run examples/2D_poiseuille_thickening_nn/case.py -n 2")

    ycb = np.fromfile(os.path.join(RESTART_DIR, "lustre_y_cb.dat"), dtype=np.float64)
    y = 0.5 * (ycb[:-1] + ycb[1:])

    u_last = read_profile(dats[-1])
    u_ana = u_analytic(y)

    if len(dats) >= 2:
        u_prev = read_profile(dats[-2])
        drift = np.sqrt(np.sum((u_last - u_prev) ** 2) / np.sum(u_last**2))
        print(f"Steady-state drift between last two saves (rel L2): {drift:.3e}")

    err = np.sqrt(np.sum((u_last - u_ana) ** 2) / np.sum(u_ana**2))
    i_peak = int(np.argmin(np.abs(y - H)))
    print(f"u_max numeric  = {u_last.max():.6e}  (analytic {u_ana.max():.6e})")
    print(f"u at walls     = {u_last[0]:.3e}, {u_last[-1]:.3e}  (analytic 0)")
    print(f"u at centerline= {u_last[i_peak]:.6e}  (analytic {u_ana[i_peak]:.6e})")
    print(f"Relative L2 error vs analytic power-law profile: {err:.4e}")

    # Bluntness: mean/peak. Parabola (Newtonian) = 2/3; a power-law n>1 profile is
    # more pointed, theory (n+1)/(2n+1).
    print(f"mean/peak (bluntness): {u_last.mean()/u_last.max():.3f}  " f"(parabola 0.667; n={NN} theory {(NN + 1) / (2 * NN + 1):.3f})")

    # Local momentum balance: K*|du/dy|^n must equal rho*g*(H - y).
    dudy = np.gradient(u_last, y)
    print("Local momentum balance K|du/dy|^n vs rho*g*(H-y):")
    for frac in (0.10, 0.20, 0.30):
        i = int(np.argmin(np.abs(y - frac * L_Y)))
        tau_num = K * np.abs(dudy[i]) ** NN
        tau_th = RHO * G_X * (H - y[i])
        print(f"  y={y[i]:.4f}  num={tau_num:.4e}  theory={tau_th:.4e}  ratio={tau_num / tau_th:.3f}")


if __name__ == "__main__":
    main()
