#!/usr/bin/env python3
"""
Density error of 1D_smooth_wave_order runs against the exact cell averages, and the observed order.

    check_smooth.py RUN_DIR [RUN_DIR ...] --case-args "ARGS_1;ARGS_2;..."

Each RUN_DIR is a finished run (pre_process, simulation, post_process) of a copy of this example; --case-args are the
case.py arguments of each run (one string for all, or one per run separated by ';'). The exact density of the gamma = 3
wave follows from the characteristics: u + c = rho0(xi) with xi + rho0(xi) t = x, u - c = -rho0(eta) with
eta - rho0(eta) t = x, and rho = c. Prints the L1 and max errors of the density (10-point Gauss cell averages).
"""

import argparse
import json
import math
import subprocess
import sys

import numpy as np

from mfc.viz import assemble_silo, discover_timesteps

A = 0.2


def rho0(x):
    return 1 + A * np.sin(2 * np.pi * x)


def exact_rho(x, t, problem):
    if problem == "advect":
        return rho0(x - t)
    xi, eta = x.copy(), x.copy()
    for _ in range(60):  # Newton on the characteristic feet
        xi -= (xi + rho0(xi) * t - x) / (1 + 2 * np.pi * A * np.cos(2 * np.pi * xi) * t)
        eta -= (eta - rho0(eta) * t - x) / (1 - 2 * np.pi * A * np.cos(2 * np.pi * eta) * t)
    return 0.5 * (rho0(xi) + rho0(eta))


def check(d, case_args):
    case = json.loads(subprocess.run([sys.executable, "case.py", *case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)
    problem = "advect" if "advect" in case_args else "gamma3"
    step = discover_timesteps(d, "silo")[-1]
    rho = np.asarray(assemble_silo(d, step).variables["alpha_rho1"]).squeeze()
    n = rho.size
    h = 1.0 / n
    xc = (np.arange(n) + 0.5) * h
    xg, wg = np.polynomial.legendre.leggauss(10)
    exact = (exact_rho((xc[:, None] + 0.5 * h * xg).ravel(), step * case["dt"], problem).reshape(n, -1) * wg).sum(-1) / 2
    return n, float(np.abs(rho - exact).mean()), float(np.abs(rho - exact).max())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--case-args", default="", help="case.py arguments of the runs; ';' separates one entry per run")
    args = ap.parse_args()
    cargs = args.case_args.split(";")
    if len(cargs) == 1:
        cargs *= len(args.runs)
    prev = None
    for d, ca in zip(args.runs, cargs):
        n, l1, linf = check(d, ca)
        order = f"  order {math.log(prev[1] / l1) / math.log(n / prev[0]):.2f}" if prev else ""
        print(f"{ca:50s} N = {n:4d}  L1 = {l1:.3e}  max = {linf:.3e}{order}")
        prev = (n, l1)


if __name__ == "__main__":
    main()
