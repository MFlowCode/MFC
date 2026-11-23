#!/usr/bin/env python3
"""
Laminar flame speed of 1D_h2_air_flame_speed against the Cantera FreeFlame value.

    compare_flame_speed.py CASE_DIR [--case-args '--dx_um 20'] [--t-avg 0.75e-3] [--tol TOL]

Reads the Silo output of CASE_DIR and the Cantera reference written to IC/reference.dat by case.py. Prints
  S_c  consumption speed, int(omega_H2) dx / (rho_u (Y_H2,u - Y_H2,out)), with int(omega_H2) = -dM_H2/dt + F_left - F_right from the H2
       mass in the domain and the convective H2 fluxes at the end cells. Y_H2,out is the outlet value: the domain ends inside the
       slow recombination zone, so using the Cantera burnt value 2 cm downstream would read about 3 % low.
  S_f  front speed, u(x_f - 1 mm) - dx_f/dt, x_f the c = 0.5 isoline of the H2 progress variable
averaged over t >= --t-avg (+- is the standard deviation in time). Exits 1 if |S_c/S_L - 1| > --tol.
"""

import argparse
import json
import os
import subprocess
import sys

import numpy as np

from mfc.viz import assemble_silo, discover_timesteps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--t-avg", type=float, default=0.75e-3, help="Average the speeds over t >= T_AVG [s].")
    ap.add_argument("--tol", type=float, default=0.01, help="Max |S_c / S_L - 1|.")
    ap.add_argument("--case-args", default="", help="Arguments given to case.py for the run, e.g. '--dx_um 20' (gives dt).")
    args = ap.parse_args()
    dt = json.loads(subprocess.run([sys.executable, "case.py", *args.case_args.split()], cwd=args.case_dir, capture_output=True, text=True, check=True).stdout)["dt"]

    SL, delta, rho_u, Yu, Yb = np.loadtxt(os.path.join(args.case_dir, "IC", "reference.dat"))
    steps = discover_timesteps(args.case_dir, "silo")
    t, M, FL, FR, Y_out, xf, uu = [], [], [], [], [], [], []
    for s in steps:
        a = assemble_silo(args.case_dir, s)
        x, rho, u, Y = a.x_cc, a.variables["alpha_rho1"], a.variables["vel1"], a.variables["Y_H2"]
        dx = x[1] - x[0]
        t.append(s * dt)
        M.append(float(np.sum(rho * Y) * dx))
        FL.append(float(rho[0] * u[0] * Y[0]))
        FR.append(float(rho[-1] * u[-1] * Y[-1]))
        Y_out.append(float(Y[-1]))
        c = (Yu - Y) / (Yu - Yb)
        i = int(np.argmax(c > 0.5))
        xf.append(x[i - 1] + (0.5 - c[i - 1]) / (c[i] - c[i - 1]) * dx)
        uu.append(float(np.interp(xf[-1] - 1.0e-3, x, u)))
    t, M, FL, FR, Y_out, xf, uu = map(np.array, (t, M, FL, FR, Y_out, xf, uu))
    wint = -np.diff(M) / np.diff(t) + 0.5 * (FL[1:] + FL[:-1]) - 0.5 * (FR[1:] + FR[:-1])
    Sc = wint / (rho_u * (Yu - 0.5 * (Y_out[1:] + Y_out[:-1])))
    Sf = uu[1:] - np.diff(xf) / np.diff(t)
    sel = 0.5 * (t[1:] + t[:-1]) >= args.t_avg
    cells = dx * 1e6
    Sc_m, Sf_m = Sc[sel].mean(), Sf[sel].mean()
    print(f"dx = {cells:.1f} um ({delta / dx:.1f} cells per Cantera thermal thickness), {t.size} saves to t = {t[-1] * 1e3:.3f} ms")
    print(f"S_L (Cantera) = {SL:.4f} m/s")
    print(f"S_c = {Sc_m:.4f} +- {Sc[sel].std():.4f} m/s, error {Sc_m / SL - 1:+.2%}")
    print(f"S_f = {Sf_m:.4f} +- {Sf[sel].std():.4f} m/s, error {Sf_m / SL - 1:+.2%}")
    ok = abs(Sc_m / SL - 1) <= args.tol
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
