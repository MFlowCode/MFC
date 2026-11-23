#!/usr/bin/env python3
"""
Compare 1D_h2_air_ignition with a Cantera constant-volume reactor.

    compare_cantera.py CASE_DIR [--case-args '--T0 1200'] [--tol-tau TOL] [--tol-T TOL]

Reads the Silo output of CASE_DIR (run pre_process, simulation and post_process first) and the case parameters from case.py, advances an
IdealGasReactor on the same mechanism and initial state, and prints
  tau   ignition delay (time of max dT/dt, parabola through the peak), MFC vs Cantera sampled at the same times
  dT    max |T_MFC - T_Cantera| over the saved times
  Teq   (T_end - T_UV-equilibrium) at the last save, and the largest relative difference of the major species from equilibrium
  sumY  max |sum(Y) - 1|
Exits 1 if |tau error| > --tol-tau or dT > --tol-T. The defaults hold for the split integrator at dt = 1e-8 s.
"""

import argparse
import json
import os
import subprocess
import sys

import cantera as ct
import numpy as np

from mfc.viz import assemble_silo, discover_timesteps

MAJORS = ("H2O", "H2", "O2", "OH", "N2")


def tau_max_slope(t, T):
    """Time of max dT/dt from samples: central differences and a parabola through the peak."""
    tm = 0.5 * (t[1:] + t[:-1])
    s = np.diff(T) / np.diff(t)
    i = int(np.clip(np.argmax(s), 1, s.size - 2))
    y0, y1, y2 = s[i - 1 : i + 2]
    den = y0 - 2 * y1 + y2
    off = 0.5 * (y0 - y2) / den if den != 0 else 0.0
    return tm[i] + off * (tm[i + 1] - tm[i])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--case-args", default="", help="Arguments given to case.py for the run, e.g. '--T0 1200'.")
    ap.add_argument("--tol-tau", type=float, default=5e-5, help="Max relative ignition-delay error.")
    ap.add_argument("--tol-T", type=float, default=1.0, help="Max |dT| [K].")
    args = ap.parse_args()

    # dt and the initial state come from case.py run with the same arguments as the simulation
    d = args.case_dir
    case = json.loads(subprocess.run([sys.executable, "case.py", *args.case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)
    dt = case["dt"]
    gas = ct.Solution(case["cantera_file"])
    y0 = np.array([case[f"patch_icpp(1)%Y({i + 1})"] for i in range(gas.n_species)])
    gas.DPY = case["patch_icpp(1)%alpha_rho(1)"], case["patch_icpp(1)%pres"], y0
    T0 = gas.T

    steps = discover_timesteps(d, "silo")
    mid = 12
    t, T, Y = [], [], []
    for s in steps:
        a = assemble_silo(d, s)
        t.append(s * dt)
        T.append(float(a.variables["T"][mid]))
        Y.append([float(a.variables[f"Y_{k}"][mid]) for k in gas.species_names])
    t, T, Y = np.array(t), np.array(T), np.array(Y)

    r = ct.IdealGasReactor(gas, clone=True)
    net = ct.ReactorNet([r])
    net.rtol, net.atol = 1e-12, 1e-22
    Tc = []
    for ti in t:
        if ti > 0:
            net.advance(ti)
        Tc.append(r.phase.T)
    Tc = np.array(Tc)

    eq = ct.Solution(case["cantera_file"])
    eq.DPY = case["patch_icpp(1)%alpha_rho(1)"], case["patch_icpp(1)%pres"], y0
    eq.equilibrate("UV")
    tau, tau_c = tau_max_slope(t, T), tau_max_slope(t, Tc)
    e_tau = tau / tau_c - 1
    dT = float(np.max(np.abs(T - Tc)))
    ysp = {k: Y[-1, gas.species_index(k)] for k in MAJORS}
    dY = max(abs(ysp[k] / eq.Y[gas.species_index(k)] - 1) for k in MAJORS)
    print(f"T0 = {T0:.0f} K, dt = {dt:.3g} s, {t.size} saves to t = {t[-1] * 1e6:.1f} us, substeps = {case['chem_params%reaction_substeps']}")
    print(f"tau: MFC {tau * 1e6:.4f} us, Cantera {tau_c * 1e6:.4f} us, error {e_tau:+.2e}")
    print(f"max |dT| = {dT:.4f} K; T_end = {T[-1]:.3f} K, UV equilibrium {eq.T:.3f} K, difference {T[-1] - eq.T:+.4f} K")
    print(f"majors vs equilibrium (max rel. diff): {dY:.2e}; max |sum(Y) - 1| = {np.max(np.abs(Y.sum(axis=1) - 1)):.2e}")
    ok = abs(e_tau) <= args.tol_tau and dT <= args.tol_T
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
