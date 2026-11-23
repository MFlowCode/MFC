#!/usr/bin/env python3
"""
Check a run of 1D_diffusion_verification (erf_species, erf_heat or stab; the mixing layer has check_mixing_layer.py).

    check_diffusion.py CASE_DIR [--case-args '--mode erf_heat --n 256'] [--tol TOL]

Run pre_process, simulation and post_process in CASE_DIR first (the output is written next to case.py). The case parameters are
read by running case.py with the same arguments as the run.

erf modes: the exact solution of the constant-coefficient diffusion equation from a step, (A/2) erfc((x - x0)/sqrt(4 D t)),
averaged over each cell, with D from Cantera at the uniform state (the mixture-averaged D of H2 in N2, or the thermal
diffusivity lambda/(rho cp) of N2). Prints L1 = sum|e| dx / (A s) and Linf = max|e| / A (s = sqrt(4 D t_end)) and exits 1 above --tol.
stab mode: amplitude of the grid-scale (Nyquist) mode of T at the last save over its first value; growth above 1 is unstable.
"""

import argparse
import json
import math
import subprocess
import sys

import cantera as ct
import numpy as np

from mfc.viz import assemble_silo, discover_timesteps


def erfc_cell_average(x_cb, x0, s, amp):
    """Cell averages of (amp/2) erfc((x - x0)/s), using the antiderivative of erfc."""
    z = (x_cb - x0) / s
    g = z * np.array([math.erfc(v) for v in z]) - np.exp(-(z**2)) / math.sqrt(math.pi)
    return 0.5 * amp * s * np.diff(g) / np.diff(x_cb)


def nyquist(f):
    return abs(np.sum((-1.0) ** np.arange(f.size) * (f - f.mean()))) / f.size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--case-args", default="", help="Arguments given to case.py for the run.")
    ap.add_argument("--tol", type=float, default=0.0, help="Fail above this L1 (erf modes), or at growth above 1 (stab).")
    args = ap.parse_args()

    d = args.case_dir
    case = json.loads(subprocess.run([sys.executable, "case.py", *args.case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)
    mode = args.case_args.split("--mode")[1].split()[0] if "--mode" in args.case_args else "erf_species"
    steps = discover_timesteps(d, "silo")
    last = assemble_silo(d, steps[-1])
    gas = ct.Solution(case["cantera_file"])

    if mode == "stab":
        first = assemble_silo(d, steps[0])
        a0, a1 = nyquist(first.variables["T"]), nyquist(last.variables["T"])
        print(f"Nyquist amplitude of T: {a0:.3e} -> {a1:.3e}, growth {a1 / a0:.3e}")
        ok = a1 / a0 < 1.0
    else:
        x_cb = np.linspace(case["x_domain%beg"], case["x_domain%end"], case["m"] + 2)
        x0 = 0.5 * (x_cb[0] + x_cb[-1])
        t = steps[-1] * case["dt"]
        gas.transport_model = "mixture-averaged"
        if mode == "erf_species":
            gas.TPY = 300.0, case["patch_icpp(1)%pres"], {"N2": 1.0}
            ih2 = gas.species_index("H2")
            amp, D = case[f"patch_icpp(1)%Y({ih2 + 1})"], gas.mix_diff_coeffs[ih2]
            field = last.variables["Y_H2"]
        else:
            gas.TPY = 300.0, case["patch_icpp(1)%pres"], {"N2": 1.0}
            amp, D = 0.01, gas.thermal_conductivity / (gas.density * gas.cp_mass)
            field = last.variables["T"] - 300.0
        s = math.sqrt(4.0 * D * t)
        e = field - erfc_cell_average(x_cb, x0, s, amp)
        dx = np.diff(x_cb)
        l1, linf = float(np.sum(np.abs(e) * dx) / (amp * s)), float(np.max(np.abs(e)) / amp)
        print(f"{mode}: N = {case['m'] + 1}, t = {t:.4e} s, front width sqrt(4Dt) = {s / dx[0]:.2f} cells, D = {D:.5e} m^2/s")
        print(f"L1 = {l1:.4e}, Linf = {linf:.4e} (fractions of the step, as the L1 of the profile over A s)")
        ok = l1 <= args.tol if args.tol > 0 else True
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
