#!/usr/bin/env python3
"""
Compare the final IB state of 2D_mibm_constant_acceleration with x = g t^2 / 2, v = g t.

    check_kinematics.py CASE_DIR [--g G] [--tol TOL]

Reads restart_data/ib_state_<last>.dat (written with ib_state_wrt) and the time step from simulation.inp. Prints the relative
position and velocity errors and exits 1 if either exceeds --tol.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

NFIELDS = 22  # ib_state_nfields: time, force, torque, velocity (8:10), angular velocity, angles, centroid (17:19), ...

parser = argparse.ArgumentParser()
parser.add_argument("case_dir")
parser.add_argument("--g", type=float, default=1.0e5)
parser.add_argument("--tol", type=float, default=1.0e-10)
args = parser.parse_args()

saves = glob.glob(os.path.join(args.case_dir, "restart_data", "ib_state_*.dat"))
last = max(saves, key=lambda f: int(re.search(r"ib_state_(\d+)", f).group(1)))
rec = np.fromfile(last, dtype="<f8")[:NFIELDS]
t, v, x = rec[0], rec[7], rec[16]
x_ex, v_ex = 0.5 * args.g * t**2, args.g * t
ex, ev = x / x_ex - 1.0, v / v_ex - 1.0
print(f"t = {t:.6e} s  x = {x:.9e} (exact {x_ex:.9e}, rel err {ex:+.2e})  v = {v:.9e} (exact {v_ex:.9e}, rel err {ev:+.2e})")
sys.exit(0 if max(abs(ex), abs(ev)) <= args.tol else 1)
