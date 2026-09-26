#!/usr/bin/env python3
"""
Laplace pressure jump and spurious currents of the static drop. The jump is the mean
pressure over the 5x5 cells at the drop center less that at the far corner; the exact
value is sigma/R. Reads the restart files of a run of case.py.

usage: analyze.py CASE_DIR
"""

import argparse
import glob
import os
import re

import numpy as np

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("dir", help="case directory holding simulation.inp and restart_data/")
args = parser.parse_args()

with open(os.path.join(args.dir, "simulation.inp")) as f:
    inp = f.read()


def param(name):
    return float(re.search(rf"^\s*{re.escape(name)}\s*=\s*(\S+)", inp, re.M).group(1))


dt, N, sigma = param("dt"), int(param("m")) + 1, param("sigma")
gw, pw, ga = param("fluid_pp(1)%gamma"), param("fluid_pp(1)%pi_inf"), param("fluid_pp(2)%gamma")

print(f"{'t':>7} {'dp':>8} {'sigma/R':>8} {'max|u|':>9}")
for f in sorted(glob.glob(os.path.join(args.dir, "restart_data", "lustre_[0-9]*.dat")), key=lambda f: int(f.split("_")[-1][:-4])):
    # alpha_rho(1:2), mom(1:2), E, alpha(1:2), color function; x fastest
    q = np.fromfile(f, dtype=np.float64).reshape(-1, N, N)
    rho = q[0] + q[1]
    u, v = q[2] / rho, q[3] / rho
    p = (q[4] - 0.5 * rho * (u**2 + v**2) - q[5] * pw) / (q[5] * gw + q[6] * ga)
    dp = p[:5, :5].mean() - p[-5:, -5:].mean()
    print(f"{int(f.split('_')[-1][:-4]) * dt:7.3f} {dp:8.3f} {sigma / 0.15:8.3f} {np.hypot(u, v).max():9.3e}")
