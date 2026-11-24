#!/usr/bin/env python3
"""
Flame + vortex pair: integral metrics per save and grid, observed order and Richardson value; flat-flame speeds against Cantera.

    analyze_fv.py RUN64 RUN128 RUN256 [--times 40 70 100]      (run directories of case.py at N = 64, 128, 256; Silo output)
    analyze_fv.py RUN --flat                                   (--gamma 0 strip run: S_d and S_c against the Cantera S_L)

Per save: heat release Q = int qdot dA [W/m] (Cantera rates at each cell state), H2O mass [kg/m], length of the T = 1200 K isoline [m]
(periodic seam closed). Order from the three finest grids (ratio 2): p = log2((f2 - f1)/(f3 - f2)), Richardson value f3 + (f3 - f2)/(2^p - 1).
Flat flame: displacement speed S_d = (u_a - dx_f/dt) rho_a/rho_u (x_f the y-mean T = 1200 K position, linear fit over the window; u_a and rho_a at T_up + 10 K, window means) against
S_L, and the consumption speed S_c = -int w_H2 dA / (rho_u Y_H2,u L_y) against the Cantera flame over the same window.
"""

import argparse
import os
import re

import cantera as ct
import matplotlib
import numpy as np
from fv_common import MECH, L, flame_on, gas_unburnt, profile

from mfc.viz import assemble_silo, discover_timesteps

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def inp(run, key):
    txt = open(os.path.join(run, "simulation.inp")).read()
    return float(re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", txt, re.M).group(1).replace("d", "e").replace("D", "e"))


def state(run, s, gas):
    a = assemble_silo(run, s)
    Y = np.array([a.variables[f"Y_{n}"] for n in gas.species_names])
    return a.x_cc, a.y_cc, a.variables["temp"] if "temp" in a.variables else a.variables["T"], a.variables["pres"], a.variables["vel1"], a.variables["alpha_rho1"], Y


def metrics(run, s, gas):
    x, y, T, p, u, rho, Y = state(run, s, gas)
    dA = (x[1] - x[0]) * (y[1] - y[0])
    sa = ct.SolutionArray(gas, T.size)
    sa.TPY = T.ravel(), p.ravel(), Y.reshape(gas.n_species, -1).T
    q = sa.heat_release_rate.reshape(T.shape)
    yy, TT = np.concatenate([y, y[:1] + y[-1] + y[0]]), np.concatenate([T, T[:, :1]], axis=1)  # close the periodic seam
    cs = plt.contour(x, yy, TT.T, [1200.0])
    Lf = sum(np.hypot(*np.diff(seg, axis=0).T).sum() for seg in cs.allsegs[0])
    plt.close()
    return dict(Q=q.sum() * dA, mH2O=(rho * Y[gas.species_index("H2O")]).sum() * dA, Lf=Lf)


def order(f):
    e1, e2 = f[1] - f[0], f[2] - f[1]
    if e1 * e2 <= 0:
        return np.nan, np.nan
    p = np.log2(e1 / e2)
    return p, f[2] + e2 / (2**p - 1)


def convergence(runs, times):
    gas = ct.Solution(MECH)
    res = {}
    for r in runs:
        N, dt, nsave = int(inp(r, "m")) + 1, inp(r, "dt"), int(inp(r, "t_step_save"))
        steps = {round(s * dt * 1e6): s for s in discover_timesteps(r, "silo")}
        res[N] = {t: metrics(r, steps[t], gas) for t in times if t in steps}
    Ns = sorted(res)
    print("t [us]  quantity  " + "  ".join(f"N={N}" for N in Ns) + "  p  Richardson")
    for t in times:
        for k in ("Q", "mH2O", "Lf"):
            f = [res[N][t][k] for N in Ns if t in res[N]]
            if len(f) >= 3:
                p, fr = order(f[-3:])
                print(f"{t:5d}  {k:5s}  " + "  ".join(f"{v:.6g}" for v in f) + f"  p = {p:.2f}  R = {fr:.6g}  (finest vs R: {f[-1] / fr - 1:+.1e})")
            else:
                print(f"{t:5d}  {k:5s}  " + "  ".join(f"{v:.6g}" for v in f))


def flat(run, window):
    gas, g = ct.Solution(MECH), gas_unburnt()
    f = profile()
    rho_u, Yu, k = g.density, g.Y[g.species_index("H2")], gas.species_index("H2")
    dt = inp(run, "dt")
    rows = []
    for s in discover_timesteps(run, "silo"):
        x, y, T, p, u, rho, Y = state(run, s, gas)
        sa = ct.SolutionArray(gas, T.size)
        sa.TPY = T.ravel(), p.ravel(), Y.reshape(gas.n_species, -1).T
        w = (sa.net_production_rates[:, k] * gas.molecular_weights[k]).reshape(T.shape)
        Tm, um, rm = T.mean(1), u.mean(1), rho.mean(1)
        xf, xa = np.interp(1200.0, Tm, x), np.interp(Tm[2] + 10.0, Tm, x)
        Sc = -w.sum() * (x[1] - x[0]) / (rho_u * Yu * len(y))
        rows.append([s * dt, p.mean(1)[2], xf, np.interp(xa, x, um), np.interp(xa, x, rm), Sc])
    a = np.array(rows)
    t = a[:, 0]
    sel = (t >= window[0] * 1e-6) & (t <= window[1] * 1e-6)
    slope = np.polyfit(t[sel], a[sel, 2], 1)[0]  # dx_f/dt over the window
    Sd = (a[sel, 3].mean() - slope) * a[sel, 4].mean() / rho_u
    print(f"S_L (Cantera) = {f['SL']:.4f} m/s")
    print(f"S_d over {window[0]:g}-{window[1]:g} us ({sel.sum()} saves): {Sd:.4f} m/s ({Sd / f['SL'] - 1:+.2%})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--flat", action="store_true")
    ap.add_argument("--times", type=int, nargs="+", default=[40, 70, 100], help="Save times [us] for the convergence table.")
    ap.add_argument("--window", type=float, nargs=2, default=[30.0, 70.0], help="Fit window [us] of the flat-flame speed.")
    a = ap.parse_args()
    flat(a.runs[0], a.window) if a.flat else convergence(a.runs, a.times)
