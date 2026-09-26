#!/usr/bin/env python3
"""
Kinetic energy and velocity error of the 2D Taylor-Green vortex against the exact decay
u(t) = u(0) exp(-2 nu t). Reads the restart files of a run of case.py.

usage: analyze.py CASE_DIR
"""

import argparse
import glob
import math
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


dt, N = param("dt"), int(param("m")) + 1


def state(step):
    """alpha_rho, mom_x, mom_y, E, alpha over the global grid, x fastest."""
    q = np.fromfile(os.path.join(args.dir, "restart_data", f"lustre_{step}.dat"), dtype=np.float64).reshape(5, N, N)
    return q[1] / q[0], q[2] / q[0], q[0]


steps = sorted(int(f.split("_")[-1][:-4]) for f in glob.glob(os.path.join(args.dir, "restart_data", "lustre_[0-9]*.dat")))
u0, v0, r0 = state(steps[0])
nu = 1.0 / param("fluid_pp(1)%Re(1)") / r0.mean()
ke0 = np.sum(r0 * (u0**2 + v0**2))
print(f"{'t':>7} {'KE/KE0':>9} {'exact':>9} {'rel err':>9} {'u L2 err':>9}")
for s in steps:
    u, v, r = state(s)
    t = s * dt
    decay = math.exp(-2 * nu * t)
    ke = np.sum(r * (u**2 + v**2)) / ke0
    err = math.sqrt(np.sum((u - u0 * decay) ** 2 + (v - v0 * decay) ** 2) / np.sum(u0**2 + v0**2)) / decay
    print(f"{t:7.3f} {ke:9.5f} {decay**2:9.5f} {ke / decay**2 - 1:9.2e} {err:9.2e}")
