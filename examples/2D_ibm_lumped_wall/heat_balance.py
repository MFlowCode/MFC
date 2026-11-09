#!/usr/bin/env python3
"""Plot the lumped cylinder's temperature and heat rate from D/ib_heat.dat and check its heat balance.

Run from this directory after the simulation; the file is plain text, so no post_process step is needed.
"""

import math
import os

import numpy as np

COLUMNS = ("id", "time", "mdot", "mdot_face", "qdot", "Twall")


def read_ib_heat(path, num_ibs=1):
    """Rows of D/ib_heat.dat as an array (rows, num_ibs, 6) in the COLUMNS order; bodies share a row in id order."""
    return np.atleast_2d(np.loadtxt(path)).reshape(-1, num_ibs, len(COLUMNS))


def balance_residual(time, qdot, twall, heat_capacity):
    """Per step, m c_s dT - Q dt [J] for a body with no heat_power, radiation or mass loss.

    Row k carries the heat over the step that ended at time[k] and Twall at its end, so the step's dt is time[k] - time[k-1].
    """
    return heat_capacity * np.diff(twall) - qdot[1:] * np.diff(time)


def main():
    import case as c
    import matplotlib.pyplot as plt

    p = c.case
    here = os.path.dirname(os.path.abspath(__file__))
    rows = read_ib_heat(os.path.join(here, "D", "ib_heat.dat"))[:, 0, :]
    t, q, tw = rows[:, 1], rows[:, 4], rows[:, 5]

    radius = p["patch_ib(1)%radius"]
    heat_capacity = p["patch_ib(1)%rho_solid"] * p["patch_ib(1)%cp_solid"] * math.pi * radius**2
    residual = balance_residual(t, q, tw, heat_capacity)
    energy = heat_capacity * (tw[0] - tw[-1])
    print(f"Twall: {tw[0]:.1f} K -> {tw[-1]:.1f} K over {t[-1] - t[0]:.3e} s")
    print(f"heat lost by the body: m c_s dT = {energy:.4e} J/m, -sum(Q dt) = {-np.sum(q[1:] * np.diff(t)):.4e} J/m")
    print(f"largest per-step imbalance: {np.max(np.abs(residual)):.2e} J/m")

    t_inf, t0 = c.T_inf, p["patch_ib(1)%Twall"]
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.5), constrained_layout=True)
    ax[0].plot(t * 1e3, tw, label="MFC")
    ax[0].plot(t * 1e3, t_inf + (t0 - t_inf) * np.exp(-t / c.tau), "--", label=r"$T_\infty + (T_0 - T_\infty)e^{-t/\tau}$ (estimate)")
    ax[0].set(xlabel="t [ms]", ylabel="Twall [K]")
    ax[0].legend()
    ax[1].plot(t * 1e3, -q)
    ax[1].set(xlabel="t [ms]", ylabel="heat lost by the body [W/m]")
    fig.savefig(os.path.join(here, "lumped_wall.png"), dpi=150)


if __name__ == "__main__":
    main()
