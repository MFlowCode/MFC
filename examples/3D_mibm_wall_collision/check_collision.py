#!/usr/bin/env python3
"""
Rebound of 3D_mibm_wall_collision: normal restitution and tangential velocity ratio from the last IB state.

    check_collision.py CASE_DIR [--angle DEG] [--e E] [--mu MU] [--speed V] [--tol TOL]

Reports e_n = -vy_out / vy_in against the prescribed e, and the tangential speed against the rigid-body result (sliding or
sticking, whichever applies). Exits 1 if |e_n/e - 1| > --tol.
"""

import argparse
import glob
import math
import os
import re
import sys

import numpy as np

NFIELDS = 22  # ib_state_nfields; velocity is fields 8:10

parser = argparse.ArgumentParser()
parser.add_argument("case_dir")
parser.add_argument("--angle", type=float, default=0.0)
parser.add_argument("--e", type=float, default=0.8)
parser.add_argument("--mu", type=float, default=0.0)
parser.add_argument("--speed", type=float, default=3.9)
parser.add_argument("--reflective", action="store_true", help="y = 0 is a reflective face (--wall-bc -2): frictionless, v_t' = v_t.")
parser.add_argument("--tol", type=float, default=0.02)
args = parser.parse_args()

saves = glob.glob(os.path.join(args.case_dir, "restart_data", "ib_state_*.dat"))
last = max(saves, key=lambda f: int(re.search(r"ib_state_(\d+)", f).group(1)))
vx, vy = np.fromfile(last, dtype="<f8")[:NFIELDS][7:9]

th = math.radians(args.angle)
vt_in, vn_in = args.speed * math.sin(th), args.speed * math.cos(th)
en = -vy / -vn_in
err = en / args.e - 1.0
print(f"e_n = {en:.5f} (prescribed {args.e}, error {100 * err:+.2f} %)")
if vt_in > 0.0:
    mu = 0.0 if args.reflective else args.mu
    sliding = math.tan(th) >= 3.5 * mu * (1.0 + args.e)
    ref = vt_in - mu * (1.0 + args.e) * vn_in if sliding else 5.0 / 7.0 * vt_in
    print(f"v_t out = {vx:.5f} (rigid-body {'sliding' if sliding else 'sticking'}: {ref:.5f}, error {100 * (vx / ref - 1.0):+.2f} %)")
sys.exit(0 if abs(err) <= args.tol else 1)
