#!/usr/bin/env python3
"""
Thermodynamic equilibrium of the relaxed state in 1D_exp_tube_phasechange (relax_model = 6, pTg equilibrium).

    check_saturation.py CASE_DIR [--case-args '-me 3'] [--tol TOL]

CASE_DIR is a finished run (pre_process, simulation, post_process). Fluids 1, 2, 3 are liquid water, water vapor and air,
each a stiffened gas with p + pi = (gamma - 1) cv rho T. In every cell of a save it takes p and each phase density
rho_k = alpha_rho_k/alpha_k and prints
  dT     max |T_k/T_vapor - 1| over the phases (thermal equilibrium)
  dTsat  max |T - T_sat(p)| [K] where liquid and vapor coexist, T_sat(p) the root of g_liquid(p, T) = g_vapor(p, T)
         with g = gamma cv T + q - T (cv ln(T^gamma/(p + pi)^(gamma - 1)) + q') (chemical equilibrium)
  sum    max |sum alpha_k - 1|
for the first and the last save, with the vapor mass per unit area and the pressure range. T_sat(1 atm) of the fluid
parameters is printed against the measured 373.12 K. Exits 1 if dTsat at the last save exceeds --tol.
"""

import argparse
import json
import math
import subprocess
import sys

import numpy as np
from scipy.optimize import brentq

from mfc.viz import assemble_silo, discover_timesteps


def fluids(case):
    out = []
    for i in (1, 2, 3):
        g = 1.0 + 1.0 / case[f"fluid_pp({i})%gamma"]
        out.append(dict(g=g, pi=case[f"fluid_pp({i})%pi_inf"] * (g - 1.0) / g, cv=case[f"fluid_pp({i})%cv"], q=case[f"fluid_pp({i})%qv"], qp=case[f"fluid_pp({i})%qvp"]))
    return out


def gibbs(f, p, T):
    return f["g"] * f["cv"] * T + f["q"] - T * (f["cv"] * math.log(T ** f["g"] / (p + f["pi"]) ** (f["g"] - 1.0)) + f["qp"])


def t_sat(fl, fv, p):
    return brentq(lambda T: gibbs(fl, p, T) - gibbs(fv, p, T), 200.0, 1000.0, xtol=1e-12)


def state(d, step, fl):
    V = {k: np.asarray(v).ravel() for k, v in assemble_silo(d, step).variables.items()}
    p = V["pres"]
    T = [(p + f["pi"]) * V[f"alpha{i + 1}"] / ((f["g"] - 1.0) * f["cv"] * V[f"alpha_rho{i + 1}"]) for i, f in enumerate(fl)]
    both = (V["alpha1"] > 1e-6) & (V["alpha2"] > 1e-6)
    dsat = max((abs(T[1][j] - t_sat(fl[0], fl[1], p[j])) for j in np.flatnonzero(both)), default=float("nan"))
    return dict(
        dT=float(max(np.abs(T[i] / T[1] - 1.0).max() for i in (0, 2))),
        dsat=dsat,
        sum=float(np.abs(V["alpha1"] + V["alpha2"] + V["alpha3"] - 1.0).max()),
        mv=V["alpha_rho2"],
        p=(float(p.min()), float(p.max())),
        T=(float(T[1].min()), float(T[1].max())),
        both=int(both.sum()),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--case-args", default="", help="Arguments given to case.py for the run.")
    ap.add_argument("--tol", type=float, default=0.0, help="Fail if dTsat [K] at the last save exceeds this (0: no check).")
    args = ap.parse_args()

    d = args.case_dir
    case = json.loads(subprocess.run([sys.executable, "case.py", *args.case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)
    fl = fluids(case)
    dx = (case["x_domain%end"] - case["x_domain%beg"]) / (case["m"] + 1)
    print(f"T_sat(101325 Pa) of the fluid parameters: {t_sat(fl[0], fl[1], 101325.0):.2f} K (measured 373.12 K)")
    steps = discover_timesteps(d, "silo")
    for s in (steps[0], steps[-1]):
        r = state(d, s, fl)
        print(
            f"step {s:6d} (t = {s * case['dt']:.3e} s): dT {r['dT']:.2e}, dTsat {r['dsat']:.3e} K over {r['both']} two-phase cells, sum alpha - 1 {r['sum']:.1e}, "
            f"p {r['p'][0]:.5e} .. {r['p'][1]:.5e} Pa, T {r['T'][0]:.3f} .. {r['T'][1]:.3f} K, vapor mass {r['mv'].sum() * dx:.6e} kg/m^2"
        )
    ok = args.tol <= 0 or r["dsat"] <= args.tol
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
