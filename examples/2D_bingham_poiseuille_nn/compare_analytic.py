#!/usr/bin/env python3
"""
Compare the steady MFC velocity profile of the 2D Bingham Poiseuille channel
(examples/2D_bingham_poiseuille_nn/case.py) against the closed-form analytic
plug-flow solution, and measure the numerical plug half-width.

For a Bingham fluid (Herschel-Bulkley with flow index n = 1, consistency K = mu,
yield stress tau0) driven by body acceleration g in a channel of height L_y
(half-height H = L_y/2, centerline y = H), with tau_w = rho*g*H > tau0:

  plug half-width  : y0   = tau0/(rho*g)
  sheared region   : u(y) = (1/(2*mu*rho*g)) *
                            [ (tau_w - tau0)^2 - (rho*g*(H-y) - tau0)^2 ]   for |y-H| >= y0
  plug             : u_plug = (1/(2*mu*rho*g)) * (tau_w - tau0)^2           for |y-H| <  y0

The key yield-term signature is a flat PLUG of uniform velocity within |y-H| < y0.

Reads the raw restart binaries (restart_data/lustre_*.dat) written by the run
(conservative variables, Fortran order, x fastest), extracts the x-averaged u(y)
from the last save, computes the analytic profile, and reports the relative L2
error and the measured plug half-width (region near the centerline within 1% of
u_max) versus y0. Steady state is confirmed by comparing the last two saves.

Run:  ./build/venv/bin/python3 examples/2D_bingham_poiseuille_nn/compare_analytic.py
"""

import glob
import os

import numpy as np

# Must match case.py
RHO = 1.0
G_X = 0.1
K = 5.0e-2  # n = 1 -> mu = K
TAU0 = 4.0e-3
L_X = 0.2
L_Y = 0.2
H = 0.5 * L_Y

TAU_W = RHO * G_X * H
Y0 = TAU0 / (RHO * G_X)  # analytic plug half-width

# Grid (m, n in case.py -> m+1, n+1 cells)
NX = 25
NY = 64
# Conservative variables per cell: alpha_rho(1), mom_x, mom_y, E, alpha(1)
NVAR = 5

HERE = os.path.dirname(os.path.abspath(__file__))
RESTART_DIR = os.path.join(HERE, "restart_data")


def u_analytic(y):
    """Closed-form steady Bingham plug-flow profile."""
    s = np.abs(y - H)  # distance from centerline
    tau_local = RHO * G_X * s  # |shear stress| = rho*g*s: 0 at center, tau_w at walls
    u_plug = (TAU_W - TAU0) ** 2 / (2.0 * K * RHO * G_X)
    # Sheared region where s >= y0, i.e. tau_local >= tau0
    u_shear = ((TAU_W - TAU0) ** 2 - np.maximum(tau_local - TAU0, 0.0) ** 2) / (2.0 * K * RHO * G_X)
    return np.where(s < Y0, u_plug, u_shear)


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
        raise SystemExit(f"No restart files in {RESTART_DIR}. " "Run: ./mfc.sh run examples/2D_bingham_poiseuille_nn/case.py -n 2")

    ycb = np.fromfile(os.path.join(RESTART_DIR, "lustre_y_cb.dat"), dtype=np.float64)
    y = 0.5 * (ycb[:-1] + ycb[1:])

    u_last = read_profile(dats[-1])
    u_ana = u_analytic(y)

    if len(dats) >= 2:
        u_prev = read_profile(dats[-2])
        drift = np.sqrt(np.sum((u_last - u_prev) ** 2) / np.sum(u_last**2))
        print(f"Steady-state drift between last two saves (rel L2): {drift:.3e}")

    err = np.sqrt(np.sum((u_last - u_ana) ** 2) / np.sum(u_ana**2))
    umax = u_last.max()
    i_peak = int(np.argmin(np.abs(y - H)))
    print(f"tau_w = rho*g*H = {TAU_W:.4e}   tau0 = {TAU0:.4e}   (need tau_w > tau0: {TAU_W > TAU0})")
    print(f"u_max numeric  = {umax:.6e}  (analytic plug {u_ana.max():.6e})")
    print(f"u at walls     = {u_last[0]:.3e}, {u_last[-1]:.3e}  (analytic 0)")
    print(f"u at centerline= {u_last[i_peak]:.6e}  (analytic {u_ana[i_peak]:.6e})")
    print(f"Relative L2 error vs analytic Bingham profile: {err:.4e}")

    # Measure numerical plug half-width: contiguous band around centerline where
    # u >= 0.99 * u_max. Half-width = max distance from centerline in that band.
    plug_mask = u_last >= 0.99 * umax
    plug_y = y[plug_mask]
    if plug_y.size:
        meas_half = max(abs(plug_y.max() - H), abs(plug_y.min() - H))
    else:
        meas_half = 0.0
    print(f"Plug half-width: measured (>=99% u_max) = {meas_half:.4e}   analytic y0 = {Y0:.4e}   (= {Y0 / H:.2f} H)")
    print(f"  ratio measured/analytic = {meas_half / Y0:.3f}")

    # Near-wall momentum balance for the sheared region: mu*|du/dy| + tau0 = rho*g*(H-y).
    dudy = np.gradient(u_last, y)
    print("Sheared-region balance mu|du/dy|+tau0 vs rho*g*(H-y):")
    for frac in (0.05, 0.10, 0.20):
        i = int(np.argmin(np.abs(y - frac * L_Y)))
        tau_num = K * np.abs(dudy[i]) + TAU0
        tau_th = RHO * G_X * (H - y[i])
        print(f"  y={y[i]:.4f}  num={tau_num:.4e}  theory={tau_th:.4e}  ratio={tau_num / tau_th:.3f}")


if __name__ == "__main__":
    main()
