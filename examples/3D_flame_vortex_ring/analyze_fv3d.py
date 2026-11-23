#!/usr/bin/env python3
"""
3D flame + vortex ring: integral metrics per grid at given times, observed order and Richardson value.

    analyze_fv3d.py RUN64 RUN128 [RUN256] [--times 40 70 100]      (run directories of case.py, Silo output)

Per time: H2O mass [kg], flame area A_f of the T = 1200 K isosurface [m^2] (marching cubes, periodic seams closed), heat release
Q = int qdot dV [W] (Cantera rates at each cell state). With three grids of ratio 2: p = log2((f2 - f1)/(f3 - f2)) and the Richardson value
f3 + (f3 - f2)/(2^p - 1). Needs scikit-image.
"""

import argparse
import importlib.util
import os
import re

import cantera as ct
import numpy as np
from skimage.measure import marching_cubes

from mfc.viz import assemble_silo, discover_timesteps

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("fv_common", os.path.join(HERE, "..", "2D_flame_vortex_pair", "fv_common.py"))
fv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fv)


def inp(run, key):
    txt = open(os.path.join(run, "simulation.inp")).read()
    return float(re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", txt, re.M).group(1).replace("d", "e").replace("D", "e"))


def metrics(run, s, gas):
    a = assemble_silo(run, s)
    dx = a.x_cc[1] - a.x_cc[0]
    T, p = a.variables["T"], a.variables["pres"]
    Y = np.array([a.variables[f"Y_{n}"] for n in gas.species_names])
    sa = ct.SolutionArray(gas, T.size)
    sa.TPY = T.ravel(), p.ravel(), Y.reshape(gas.n_species, -1).T
    q = sa.heat_release_rate
    Tp = np.concatenate([T, T[:, :1]], axis=1)
    Tp = np.concatenate([Tp, Tp[:, :, :1]], axis=2)  # close the periodic y and z seams
    v, f, _, _ = marching_cubes(Tp, 1200.0, spacing=(dx, dx, dx))
    u, w, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(w - u, c - u), axis=1).sum()
    return dict(mH2O=(a.variables["alpha_rho1"] * Y[gas.species_index("H2O")]).sum() * dx**3, Af=area, Q=q.sum() * dx**3)


def order(f):
    e1, e2 = f[1] - f[0], f[2] - f[1]
    if e1 * e2 <= 0:
        return np.nan, np.nan
    p = np.log2(e1 / e2)
    return p, f[2] + e2 / (2**p - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--times", type=int, nargs="+", default=[40, 70, 100], help="Save times [us].")
    args = ap.parse_args()
    gas = ct.Solution(fv.MECH)
    res = {}
    for r in args.runs:
        N, dt = int(inp(r, "m")) + 1, inp(r, "dt")
        steps = {round(s * dt * 1e6): s for s in discover_timesteps(r, "silo")}
        res[N] = {t: metrics(r, steps[t], gas) for t in args.times if t in steps}
    Ns = sorted(res)
    print("t [us]  quantity  " + "  ".join(f"N={N}" for N in Ns))
    for t in args.times:
        for k in ("mH2O", "Af", "Q"):
            f = [res[N][t][k] for N in Ns if t in res[N]]
            line = f"{t:5d}  {k:5s}  " + "  ".join(f"{v:.6g}" for v in f)
            if len(f) >= 3:
                p, fr = order(f[-3:])
                line += f"  p = {p:.2f}  R = {fr:.6g}  (finest vs R: {f[-1] / fr - 1:+.1e})"
            print(line)


if __name__ == "__main__":
    main()
