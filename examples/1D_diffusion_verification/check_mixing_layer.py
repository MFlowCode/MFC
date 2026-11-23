#!/usr/bin/env python3
"""
Error of the mixing layer of case.py --mode mixing against the reference of ref_mixing_layer.py.

    check_mixing_layer.py CASE_DIR REF.npz [--case-args '--mode mixing --n 300']

Run pre_process, simulation and post_process in CASE_DIR first. At the 10 and 20 us saves the MFC profiles are compared in the mass
coordinate psi = int rho dx from the wall at x = 0, which removes the uniform compression a small pressure offset of the outflow
causes (it is not a diffusion error). The reference is interpolated linearly in psi. For species k, L1 = mean|e| L / (|dY_k| L_layer)
and Linf = max|e| / |dY_k|, with dY_k the jump of Y_k across the layer, L = 3 mm, L_layer = 1 mm (the same for T, jump 700 K). MFC's T is first corrected
isentropically to the reference pressure. Also prints the pressure deviation and the closure max|sum Y - 1|.
"""

import argparse
import json
import subprocess
import sys

import cantera as ct
import numpy as np
from ref_mixing_layer import LEFT, LX, P0, RIGHT, X0

from mfc.viz import assemble_silo, discover_timesteps

SP = ("H2", "O2", "H2O", "N2")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("ref")
    ap.add_argument("--case-args", default="--mode mixing")
    args = ap.parse_args()
    case = json.loads(subprocess.run([sys.executable, "case.py", *args.case_args.split()], cwd=args.case_dir, capture_output=True, text=True, check=True).stdout)
    ref = np.load(args.ref)
    g = ct.Solution(case["cantera_file"])
    YL = ct.Solution(case["cantera_file"])
    YR = ct.Solution(case["cantera_file"])
    YL.TPX, YR.TPX = (LEFT[0], P0, LEFT[1]), (RIGHT[0], P0, RIGHT[1])
    steps = discover_timesteps(args.case_dir, "silo")[1:]
    n = case["m"] + 1
    print(f"N = {n}, saves {steps}, reference {ref['psi'].size} cells")
    print("t [us]   " + "  ".join(f"L1_{s}" for s in SP) + "      L1_T       Linf_T     max|p/p0-1|  max|sumY-1|")
    for i, s in enumerate(steps):
        a = assemble_silo(args.case_dir, s)
        rho = a.variables["alpha_rho1"]
        dx = a.x_cc[1] - a.x_cc[0]
        psi = np.cumsum(rho * dx) - 0.5 * rho * dx
        Y = np.stack([a.variables[f"Y_{k}"] for k in g.species_names], axis=1)
        gg = ct.SolutionArray(g, n)
        gg.TPY = a.variables["T"], a.variables["pres"], np.clip(Y, 0, 1)
        Tc = a.variables["T"] * (P0 / a.variables["pres"]) ** (ct.gas_constant / gg.mean_molecular_weight / gg.cp_mass)
        row = f"{ref['t'][i] * 1e6:6.1f}   "
        for k in SP:
            j = g.species_index(k)
            e = Y[:, j] - np.interp(psi, ref["psi"], ref["Y"][i][:, j])
            row += f"{np.abs(e).mean() * LX / (abs(YR.Y[j] - YL.Y[j]) * X0):.3e}  "
        e = Tc - np.interp(psi, ref["psi"], ref["T"][i])
        dT = RIGHT[0] - LEFT[0]
        row += f"{np.abs(e).mean() * LX / (dT * X0):.3e}  {np.abs(e).max() / dT:.3e}  {np.abs(a.variables['pres'] / P0 - 1).max():.3e}  {np.abs(Y.sum(1) - 1).max():.2e}"
        print(row)


if __name__ == "__main__":
    main()
