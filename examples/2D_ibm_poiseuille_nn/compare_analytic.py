#!/usr/bin/env python3
"""
Validation B (and C) for examples/2D_ibm_poiseuille_nn: compare the steady MFC
velocity profile in the gap between the two IBM wall slabs against the
closed-form power-law Poiseuille solution

    u(y) = (n/(n+1)) * (rho*g/K)^(1/n) * ( H^((n+1)/n) - |y - y_c|^((n+1)/n) )

using both the nominal gap geometry (slab faces at y = 0.05, 0.25 -> H = 0.1,
y_c = 0.15) and a fitted geometry (H, y_c from the u -> 0 crossings adjacent to
each wall, since IBM walls are sharp only to ~half a cell).

Validation C: the analytic steady x-force per wall per unit depth is
tau_w * L_x with tau_w = rho*g*H; compared against the IBM volume-integrated
force in restart_data/ib_state_<last>.dat (written because ib_state_wrt = T;
patch_ib%mass = 0 so the record is the pure pressure+viscous integration).

Run:  ./mfc.sh run examples/2D_ibm_poiseuille_nn/case.py -n 2
      ./build/venv/bin/python3 examples/2D_ibm_poiseuille_nn/compare_analytic.py
"""

import glob
import os

import numpy as np

# Must match case.py (powerlaw mode)
RHO = 1.0
G_X = 5.0e-2
K = 5.0e-2
NN = 1.5
L_X = 0.2
Y_LO, Y_HI = 0.05, 0.25  # nominal slab gap faces
NX = 25
NY = 96
NVAR = 5  # alpha_rho(1), mom_x, mom_y, E, alpha(1)

HERE = os.path.dirname(os.path.abspath(__file__))
RESTART_DIR = os.path.join(HERE, "restart_data")


def u_analytic(y, h, y_c):
    expo = (NN + 1.0) / NN
    pref = (NN / (NN + 1.0)) * (RHO * G_X / K) ** (1.0 / NN)
    return pref * (h**expo - np.abs(y - y_c) ** expo)


def read_profile(path):
    raw = np.fromfile(path, dtype=np.float64).reshape((NVAR, NY, NX))
    return (raw[1] / raw[0]).mean(axis=1)  # x-averaged u(y)


def zero_crossing(y0, u0, y1, u1):
    """Linear-interpolated y where u crosses 0 between two samples."""
    return y0 - u0 * (y1 - y0) / (u1 - u0)


def main():
    dats = sorted(
        glob.glob(os.path.join(RESTART_DIR, "lustre_[0-9]*.dat")),
        key=lambda p: int(os.path.basename(p).split("_")[1].split(".")[0]),
    )
    if not dats:
        raise SystemExit(f"No restart files in {RESTART_DIR}. " "Run: ./mfc.sh run examples/2D_ibm_poiseuille_nn/case.py -n 2")

    ycb = np.fromfile(os.path.join(RESTART_DIR, "lustre_y_cb.dat"), dtype=np.float64)
    y = 0.5 * (ycb[:-1] + ycb[1:])
    gap = (y > Y_LO) & (y < Y_HI)
    i_lo, i_hi = np.argmax(gap), NY - 1 - np.argmax(gap[::-1])  # first/last gap row

    u = read_profile(dats[-1])
    if len(dats) >= 2:
        u_prev = read_profile(dats[-2])
        drift = np.sqrt(np.sum((u[gap] - u_prev[gap]) ** 2) / np.sum(u[gap] ** 2))
        print(f"Steady-state drift between last two saves (rel L2, gap): {drift:.3e}")

    # Fitted wall positions from the u -> 0 crossings into the IBM ghost layers
    y_lo_fit = zero_crossing(y[i_lo - 1], u[i_lo - 1], y[i_lo], u[i_lo])
    y_hi_fit = zero_crossing(y[i_hi], u[i_hi], y[i_hi + 1], u[i_hi + 1])
    h_fit, yc_fit = 0.5 * (y_hi_fit - y_lo_fit), 0.5 * (y_hi_fit + y_lo_fit)
    print(f"Fitted walls: y = {y_lo_fit:.5f}, {y_hi_fit:.5f}  (nominal {Y_LO}, {Y_HI}; dy = {y[1] - y[0]:.5f})")

    for label, h, y_c in (
        ("nominal H = 0.1   ", 0.5 * (Y_HI - Y_LO), 0.5 * (Y_LO + Y_HI)),
        (f"fitted  H = {h_fit:.4f}", h_fit, yc_fit),
    ):
        u_ana = u_analytic(y[gap], h, y_c)
        err = np.sqrt(np.sum((u[gap] - u_ana) ** 2) / np.sum(u_ana**2))
        print(f"Relative L2 error vs analytic ({label}): {err:.4e}   " f"(u_max num {u[gap].max():.4e}, ana {u_ana.max():.4e})")

    # Bluntness: parabola (Newtonian) = 2/3; power-law theory (n+1)/(2n+1)
    print(f"mean/peak (bluntness, gap): {u[gap].mean() / u[gap].max():.3f}  " f"(parabola 0.667; n={NN} theory {(NN + 1) / (2 * NN + 1):.3f})")

    # Validation C: IBM-integrated x-force per wall vs analytic tau_w*L_x
    t_last = int(os.path.basename(dats[-1]).split("_")[1].split(".")[0])
    ib_state = os.path.join(RESTART_DIR, f"ib_state_{t_last}.dat")
    if os.path.exists(ib_state):
        rec = np.fromfile(ib_state, dtype=np.float64).reshape(2, 20)
        f_ana = RHO * G_X * 0.5 * (Y_HI - Y_LO) * L_X  # tau_w*L_x, nominal H
        print(f"IBM x-force per wall: {rec[0, 1]:.4e} (bottom), {rec[1, 1]:.4e} (top); " f"analytic tau_w*L_x = {f_ana:.4e}")
        print(f"  force ratios vs analytic: {rec[0, 1] / f_ana:.3f}, {rec[1, 1] / f_ana:.3f}")
    else:
        print(f"No {ib_state}; skipping force check")


if __name__ == "__main__":
    main()
