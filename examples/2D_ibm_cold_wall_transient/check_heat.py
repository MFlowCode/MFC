#!/usr/bin/env python3
"""
Heat absorbed by the cold cylinder of 2D_ibm_cold_wall_transient against a body-fitted 1D reference (ref1d.py).

    check_heat.py CASE_DIR [--times 0.003 0.01 0.03 0.1]

Run pre_process and simulation in CASE_DIR first. Q(t) is the running sum of column 5 of D/ib_heat.dat (energy into the body per unit
depth, the face-flux heat including the enthalpy of any exchanged mass) times the step. The reference is the low-Mach conduction
solution with the Cantera lambda(T) and c_p(T) of the same frozen air; times are in R^2/alpha_inf of the hot gas. It prints
the error (Q_MFC / Q_ref - 1) at each time the run reaches.
"""

import argparse
import os
import re

import cantera as ct
import numpy as np
from ref1d import solve_props


def read_inp(path):
    out = {}
    for line in open(path):
        m = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
        if m:
            out[m.group(1)] = m.group(2).strip("'\"")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--times", type=float, nargs="+", default=[0.003, 0.01, 0.03, 0.1])
    args = ap.parse_args()

    inp = {**read_inp(os.path.join(args.case_dir, "pre_process.inp")), **read_inp(os.path.join(args.case_dir, "simulation.inp"))}
    f = lambda k: float(inp[k].replace("d", "e"))  # noqa: E731
    R, Tw, p, delta0 = f("patch_ib(1)%radius"), f("patch_ib(1)%Twall"), f("patch_icpp(1)%pres"), f("patch_ib(1)%thermal_layer")
    dx = (f("x_domain%end") - f("x_domain%beg")) / (int(inp["m"]) + 1)

    g = ct.Solution("h2o2.yaml")
    y = [f(f"patch_icpp(1)%Y({i + 1})") for i in range(g.n_species)]
    g.DPY = f("patch_icpp(1)%alpha_rho(1)"), p, y
    Tinf = g.T
    Ts = np.linspace(150.0, 2500.0, 941)
    lam_t, cp_t = [], []
    for T in Ts:
        g.TPY = T, p, y
        lam_t.append(g.thermal_conductivity)
        cp_t.append(g.cp_mass)
    g.TPY = Tinf, p, y
    tth = R**2 * g.density * g.cp_mass / g.thermal_conductivity
    W = g.mean_molecular_weight / 1000.0

    h = np.loadtxt(os.path.join(args.case_dir, "D", "ib_heat.dat"), ndmin=2)
    h = h[np.argsort(h[:, 1])]
    t, Qd = h[:, 1], h[:, 4]
    Q = np.cumsum(Qd * np.diff(np.concatenate([[0.0], t])))
    fr = [x for x in args.times if x * tth <= t[-1] * (1 + 1e-9)]
    ref = solve_props(Tw, Tinf, lambda T: np.interp(T, Ts, lam_t), lambda T: np.interp(T, Ts, cp_t), W, p, [x * tth for x in fr], delta0, R, 8 * R)
    wc = inp.get("patch_ib(1)%wall_column", "F")
    print(f"{R / dx:.0f} cells/R, wall_column = {wc}, Tw = {Tw:.0f} K, Tinf = {Tinf:.0f} K, R^2/alpha = {tth:.3e} s, {t.size} steps")
    print("t [R^2/alpha]   Q_MFC [J/m]   Q_ref [J/m]   error")
    for x, (_, _, Qr) in zip(fr, ref):
        Qm = np.interp(x * tth, t, Q)
        print(f"{x:8.3f}   {Qm:12.5e}   {Qr:12.5e}   {Qm / Qr - 1:+.2%}")


if __name__ == "__main__":
    main()
