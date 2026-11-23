#!/usr/bin/env python3
"""
Chapman-Jouguet reference for 1D_cj_detonation: equilibrium CJ state, frozen von Neumann (vN) state, ZND induction length, and the
Taylor-Zeldovich expansion behind a CJ wave from a closed end. Cantera only; writes cj_ref.json and taylor_eq.txt next to this file.

    python3 examples/1D_cj_detonation/cj_ref.py [--p1 1e5 --T1 298 --X "H2:2, O2:1, AR:7"]

CJ: equilibrium Rankine-Hugoniot state with the sonic condition w2 = c_eq (solved for rho2, T2, D). ZND: constant-area steady
reaction zone behind the frozen shock at 1.001 D_CJ (the CJ solution is singular at the sonic point); the induction length is the
distance to the maximum temperature gradient. Taylor wave: centred C+ fan x/t = u + c along the equilibrium isentrope through the
CJ state, u - int dp/(rho c) = const, down to the u = 0 plateau.
"""

import argparse
import json
import os

import cantera as ct
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

HERE = os.path.dirname(os.path.abspath(__file__))
MECH = "h2o2.yaml"


def c_eq(gas):
    """Equilibrium sound speed: central difference of p(rho) at fixed s with equilibrium composition."""
    s, v = gas.entropy_mass, 1.0 / gas.density
    out = []
    for f in (1.0 - 1e-5, 1.0 + 1e-5):
        g = ct.Solution(MECH)
        g.TPX = gas.T, gas.P, gas.X
        g.SV = s, v * f
        g.equilibrate("SV")
        out.append((g.P, g.density))
    gas.SV = s, v
    gas.equilibrate("SV")
    return np.sqrt((out[1][0] - out[0][0]) / (out[1][1] - out[0][1]))


def eq_at(gas, X1, T, rho):
    gas.TDX = T, rho, X1
    gas.equilibrate("TV")
    return gas


def cj_state(g1):
    """Equilibrium CJ state (gas), and the wave speed D."""
    rho1, h1, p1, X1 = g1.density, g1.enthalpy_mass, g1.P, g1.X
    g = ct.Solution(MECH)

    def res(z):
        r, T, D = z[0] * rho1, z[1] * 1000.0, z[2] * 1000.0
        eq_at(g, X1, T, r)
        w2 = D * rho1 / r
        return [(g.P - p1 - rho1 * D**2 * (1.0 - rho1 / r)) / p1, (g.enthalpy_mass + 0.5 * w2**2 - h1 - 0.5 * D**2) / 1e6, (w2 - c_eq(eq_at(g, X1, T, r))) / 1000.0]

    z = fsolve(res, [1.8, 2.9, 1.9], xtol=1e-13)
    eq_at(g, X1, z[1] * 1000.0, z[0] * rho1)
    return g, z[2] * 1000.0


def frozen_shock(g1, D):
    """Frozen (vN) post-shock state for shock speed D."""
    rho1, h1, p1 = g1.density, g1.enthalpy_mass, g1.P
    g = ct.Solution(MECH)

    def res(z):
        r, T = z[0] * rho1, z[1] * 1000.0
        g.TDX = T, r, g1.X
        w2 = D * rho1 / r
        return [(g.P - p1 - rho1 * D**2 * (1 - rho1 / r)) / p1, (g.enthalpy_mass + 0.5 * w2**2 - h1 - 0.5 * D**2) / 1e6]

    z = fsolve(res, [5.0, 1.5], xtol=1e-13)
    g.TDX = z[1] * 1000.0, z[0] * rho1, g1.X
    return g


def znd_lengths(g1, D):
    """Induction and exothermic length of the ZND wave at speed D (thermicity form of the steady reaction-zone equations)."""
    vn = frozen_shock(g1, D)
    g = ct.Solution(MECH)
    W, rho1 = g.molecular_weights, g1.density
    w0 = D * rho1 / vn.density

    def rhs(t, y):
        w, T = y[0], y[1]
        Y = np.clip(y[2:], 0.0, None)
        g.TDY = T, rho1 * D / w, Y
        yd = g.net_production_rates * W / g.density
        hk = g.partial_molar_enthalpies / W
        sig = np.sum((g.mean_molecular_weight / W - hk / (g.cp_mass * T)) * yd)
        dw = w * sig / (1.0 - (w / g.sound_speed) ** 2)
        dT = -(w * dw + np.sum(hk * yd)) / g.cp_mass
        return np.concatenate([[dw, dT], yd])

    sol = solve_ivp(rhs, (0.0, 2e-5), np.concatenate([[w0, vn.T], vn.Y]), method="BDF", rtol=1e-9, atol=1e-13, dense_output=True)
    tt = np.linspace(0.0, 0.5e-5, 200001)
    y = sol.sol(tt)
    x = np.concatenate([[0.0], np.cumsum(0.5 * (y[0, 1:] + y[0, :-1]) * np.diff(tt))])
    dTdx = np.gradient(y[1], x, edge_order=2)
    i = int(np.argmax(dTdx))
    return float(x[i]), vn


def taylor(g2, uCJ, n=60):
    """Equilibrium isentrope from the CJ state to u = 0: rows p, T, rho, u, c, xi = u + c."""
    g = ct.Solution(MECH)
    g.TDX = g2.T, g2.density, g2.X
    s, Xf = g.entropy_mass, g.X
    rows, u, v0 = [], uCJ, 1.0 / g.density
    for v in v0 * np.linspace(1.0, 4.0, 4000):
        g.SVX = s, v, Xf
        g.equilibrate("SV")
        c = c_eq(g)
        if rows:
            pp, rp, cp = rows[-1][0], rows[-1][2], rows[-1][4]
            u -= 0.5 * (1 / (rp * cp) + 1 / (g.density * c)) * (pp - g.P)
        rows.append([g.P, g.T, g.density, u, c, u + c])
        if u <= 0:
            break
    a = np.array(rows)
    f = -a[-2, 3] / (a[-1, 3] - a[-2, 3])
    a[-1] = a[-2] + f * (a[-1] - a[-2])
    keep = np.unique(np.linspace(0, len(a) - 1, n).astype(int))
    return a[keep]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p1", type=float, default=1.0e5)
    ap.add_argument("--T1", type=float, default=298.0)
    ap.add_argument("--X", default="H2:2, O2:1, AR:7")
    args = ap.parse_args()
    g1 = ct.Solution(MECH)
    g1.TPX = args.T1, args.p1, args.X
    g2, D = cj_state(g1)
    uCJ = D * (1 - g1.density / g2.density)
    ind, vn = znd_lengths(g1, 1.001 * D)
    vn = frozen_shock(g1, D)
    tz = taylor(g2, uCJ)
    out = {
        "X": args.X,
        "T1": args.T1,
        "p1": args.p1,
        "rho1": g1.density,
        "D_CJ": D,
        "p_CJ": g2.P,
        "T_CJ": g2.T,
        "rho_CJ": g2.density,
        "u_CJ": uCJ,
        "Y_CJ": dict(zip(g2.species_names, g2.Y.tolist())),
        "p_vN": vn.P,
        "T_vN": vn.T,
        "rho_vN": vn.density,
        "u_vN": D * (1 - g1.density / vn.density),
        "c_vN": vn.sound_speed,
        "ind_len": ind,
        "plateau": {"p": tz[-1, 0], "T": tz[-1, 1], "rho": tz[-1, 2], "c": tz[-1, 4]},
    }
    with open(os.path.join(HERE, "cj_ref.json"), "w") as f:
        json.dump(out, f, indent=1)
    np.savetxt(os.path.join(HERE, "taylor_eq.txt"), tz, header="p T rho u c xi", comments="", fmt="%.6e")
    print(
        f"D_CJ = {D:.3f} m/s, p_CJ = {g2.P / 1e3:.1f} kPa, T_CJ = {g2.T:.0f} K, u_CJ = {uCJ:.1f} m/s; vN p {vn.P / 1e3:.0f} kPa T {vn.T:.0f} K; ind_len {ind * 1e3:.4f} mm; plateau p {tz[-1, 0] / 1e3:.1f} kPa T {tz[-1, 1]:.0f} K"
    )


if __name__ == "__main__":
    main()
