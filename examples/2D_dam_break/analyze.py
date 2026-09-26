#!/usr/bin/env python3
"""
Surge front Z = z/a against T = t*sqrt(2g/a) for the 2D dam break, compared with
Martin & Moyce (1952), table 2 (n^2 = 2, a = 2.25 in). Also reports the residual column
height H = eta/(2a) at the left wall. Reads the restart files of a run of case.py.

usage: analyze.py CASE_DIR [--plot FILE]
"""

import argparse
import glob
import math
import os
import re

import numpy as np

a, g = 0.05715, 9.81
T_exp = [0.0, 0.41, 0.84, 1.19, 1.43, 1.63, 1.82, 1.97, 2.2, 2.32, 2.5, 2.64, 2.82, 2.96]
Z_exp = [1.0, 1.11, 1.23, 1.44, 1.67, 1.89, 2.11, 2.33, 2.56, 2.78, 3.0, 3.22, 3.44, 3.67]

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("dir", help="case directory holding simulation.inp and restart_data/")
parser.add_argument("--plot", help="write a Z(T) figure to this file")
args = parser.parse_args()

with open(os.path.join(args.dir, "simulation.inp")) as f:
    inp = f.read()
# Output is numbered by step at fixed dt, by save index under cfl_adap_dt
adaptive = re.search(r"^\s*cfl_adap_dt\s*=\s*T", inp, re.M) is not None
dt = float(re.search(r"^\s*" + ("t_save" if adaptive else "dt") + r"\s*=\s*(\S+)", inp, re.M).group(1))
nx, ny = (int(re.search(rf"^\s*{v}\s*=\s*(\S+)", inp, re.M).group(1)) + 1 for v in "mn")
x, y = (np.arange(nx) + 0.5) * 5 * a / nx, (np.arange(ny) + 0.5) * 3 * a / ny


def alpha_water(step):
    """Water volume fraction from the parallel-I/O restart file: the 7 conserved fields (alpha_rho(1:2), mom(1:2), E,
    alpha(1:2)) over the global grid, x fastest."""
    q = np.fromfile(os.path.join(args.dir, "restart_data", f"lustre_{step}.dat"), dtype=np.float64)
    return q.reshape(7, ny, nx)[5]


steps = sorted(int(f.split("_")[-1][:-4]) for f in glob.glob(os.path.join(args.dir, "restart_data", "lustre_[0-9]*.dat")))
T, Z, H = [], [], []
print(f"{'T':>6} {'Z':>6} {'Z_exp':>6} {'H':>6}")
for s in steps:
    al = alpha_water(s)
    wet = al > 0.5
    T.append(s * dt * math.sqrt(2 * g / a))
    Z.append((x[wet.any(axis=0)].max() + 0.5 * (x[1] - x[0])) / a)
    H.append((y[wet[:, 0]].max() + 0.5 * (y[1] - y[0])) / (2 * a))
    print(f"{T[-1]:6.2f} {Z[-1]:6.2f} {np.interp(T[-1], T_exp, Z_exp, right=np.nan):6.2f} {H[-1]:6.2f}")

if args.plot:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(T, Z, "-", label="MFC, all-Mach projection")
    ax.plot(T_exp, Z_exp, "ko", mfc="none", label="Martin & Moyce (1952)")
    ax.set_xlabel(r"$T = t\sqrt{2g/a}$")
    ax.set_ylabel(r"$Z = z/a$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.plot, dpi=150)
