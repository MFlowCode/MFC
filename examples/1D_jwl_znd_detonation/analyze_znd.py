#!/usr/bin/env python3
"""ZND validation for the 1D jwl_reactive + jwl_delta_e detonation.

Builds the exact analytic reference for the single-fluid energy-offset JWL
model -- unreacted/product Hugoniots, Rayleigh line, CJ tangency -- from the
case's LX-10-0 parameters (Garno et al. 2020, Table I), then compares the MFC
profile against it:

  1. detonation speed (front trajectory fit)   vs analytic D_CJ (8.821 km/s)
  2. peak (von Neumann) pressure               vs analytic P_vN of this model
  3. pressure at reaction completion           vs analytic P_CJ (37.5 GPa)

The analytic CJ values are rate-independent (EOS + Q only); the spike is
grid-limited from below, so it approaches P_vN monotonically under refinement.

Usage: ./analyze_znd.py [case_dir] [dt]
  case_dir: directory containing the case's D/ output (default .)
  dt:       the case's time step, needed for the front-trajectory fit
            (default 2.0e-10, matching case.py -- pass the actual dt if
            you change it)
"""

import glob
import math
import os
import sys

import numpy as np

# Same constants as case.py (Garno Table I, LX-10-0)
A, B = 880.2e9, 17.437e9
R1, R2, omega = 4.60, 1.20, 0.30
rho0, Q, p0 = 1860.0, 5.59e6, 101325.0
v0 = 1.0 / rho0


def F(V):
    return A * (1 - omega / (R1 * V)) * np.exp(-R1 * V) + B * (1 - omega / (R2 * V)) * np.exp(-R2 * V)


delta_e = (p0 - F(1.0)) / (omega * rho0)


def hugoniot_p(V, lam):
    """p on the partial-reaction Hugoniot: RH energy jump + offset EOS, linear in p."""
    rhs = F(V) + (omega / (V * v0)) * (0.5 * p0 * v0 * (1 - V) + lam * Q + (1 - lam) * delta_e)
    return rhs / (1.0 - (omega / (2.0 * V)) * (1.0 - V))


def rayleigh_p(V, D):
    return p0 + rho0 * D**2 * (1.0 - V)


def bisect(f, a, b, tol=1e-10, it=200):
    fa = f(a)
    for _ in range(it):
        c = 0.5 * (a + b)
        fc = f(c)
        if abs(b - a) < tol:
            return c
        if fa * fc <= 0:
            b = c
        else:
            a, fa = c, fc
    return 0.5 * (a + b)


def cj_state():
    """D_CJ from tangency: min over V of (Hugoniot - Rayleigh) crosses zero in D."""
    Vg = np.linspace(0.55, 0.95, 4001)

    def gap(D):
        g = hugoniot_p(Vg, 1.0) - rayleigh_p(Vg, D)
        i = np.argmin(g)
        return g[i], Vg[i]

    D = bisect(lambda D: gap(D)[0], 6000.0, 12000.0, tol=1e-6)
    Vcj = gap(D)[1]
    return D, Vcj, rayleigh_p(Vcj, D)


def main():
    d = os.path.join(sys.argv[1] if len(sys.argv) > 1 else ".", "D")
    D_cj, V_cj, P_cj = cj_state()
    V_vn = bisect(lambda V: hugoniot_p(V, 0.0) - rayleigh_p(V, D_cj), 0.30, 0.90)
    P_vn = rayleigh_p(V_vn, D_cj)
    print(f"analytic:  D_CJ = {D_cj:.1f} m/s   P_CJ = {P_cj/1e9:.2f} GPa   " f"P_vN = {P_vn/1e9:.2f} GPa   (delta_e = {delta_e:.4e} J/kg)")

    steps = sorted({int(f.rsplit(".", 2)[1]) for f in glob.glob(f"{d}/prim.3.00.*.dat")})
    if len(steps) < 3:
        sys.exit(f"not enough saves in {d} -- run the case first")
    dt = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0e-10

    # Front trajectory from the later saves (front = rightmost cell above 2 GPa),
    # skipping early saves where the driver transient is still relaxing.
    ts, xf = [], []
    for t in steps[len(steps) // 2 :]:
        a = np.loadtxt(f"{d}/prim.3.00.{t:06d}.dat")
        x, p = a[:, 0], a[:, 1]
        ahead = np.where(p > 2.0e9)[0]
        if ahead.size:
            ts.append(t * dt)
            xf.append(x[ahead[-1]])
    D_mfc = np.polyfit(ts, xf, 1)[0]

    # Spike and CJ pressure from the final save
    t = steps[-1]
    a = np.loadtxt(f"{d}/prim.3.00.{t:06d}.dat")
    x, p = a[:, 0], a[:, 1]
    lam = np.loadtxt(f"{d}/prim.5.00.{t:06d}.dat")[:, 1]
    ipk = int(np.argmax(p))
    P_peak = p[ipk]
    burned = np.where((x < x[ipk]) & (lam > 0.99))[0]
    P_burn = p[burned[-1]] if burned.size else float("nan")

    def err(sim, ref):
        return 100.0 * (sim - ref) / ref

    print(f"MFC:       D    = {D_mfc:.1f} m/s   ({err(D_mfc, D_cj):+.2f}%)")
    print(f"           P_pk = {P_peak/1e9:.2f} GPa  ({err(P_peak, P_vn):+.2f}% of P_vN; " f"grid-limited from below)")
    print(f"           P_CJ = {P_burn/1e9:.2f} GPa  ({err(P_burn, P_cj):+.2f}%)")

    spike = P_peak > 1.05 * P_cj
    speed = abs(err(D_mfc, D_cj)) < 2.0
    cjok = abs(err(P_burn, P_cj)) < 5.0
    print(f"spike above CJ: {'PASS' if spike else 'FAIL'}   " f"front speed: {'PASS' if speed else 'FAIL'}   " f"CJ pressure: {'PASS' if cjok else 'FAIL'}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        for tt in steps[2::2]:
            aa = np.loadtxt(f"{d}/prim.3.00.{tt:06d}.dat")
            ax[0].plot(aa[:, 0] * 1e3, aa[:, 1] / 1e9, lw=1)
        ax[0].axhline(P_vn / 1e9, ls="--", c="k", lw=0.8, label=r"analytic $P_{vN}$")
        ax[0].axhline(P_cj / 1e9, ls=":", c="k", lw=0.8, label=r"analytic $P_{CJ}$")
        ax[0].set(xlabel="x [mm]", ylabel="p [GPa]", title="ZND detonation, LX-10-0")
        ax[0].legend()

        Vg = np.linspace(0.40, 1.0, 400)
        ax[1].plot(Vg, hugoniot_p(Vg, 0.0) / 1e9, label=r"Hugoniot $\lambda=0$")
        ax[1].plot(Vg, hugoniot_p(Vg, 1.0) / 1e9, label=r"Hugoniot $\lambda=1$")
        ax[1].plot(Vg, rayleigh_p(Vg, D_cj) / 1e9, "k--", label=r"Rayleigh $D_{CJ}$")
        ax[1].plot([V_cj], [P_cj / 1e9], "ko", ms=5)
        ax[1].plot([V_vn], [P_vn / 1e9], "ks", ms=5)
        ax[1].set(xlabel=r"$V = \rho_0/\rho$", ylabel="p [GPa]", ylim=(0, 90), title="Hugoniots, Rayleigh line, CJ + vN points")
        ax[1].legend()
        fig.tight_layout()
        out = os.path.join(os.path.dirname(d) or ".", "znd_profile.png")
        fig.savefig(out, dpi=150)
        print(f"wrote {out}")
    except ImportError:
        pass

    sys.exit(0 if (spike and speed and cjok) else 1)


if __name__ == "__main__":
    main()
