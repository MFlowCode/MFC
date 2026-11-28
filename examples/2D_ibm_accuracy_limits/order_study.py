#!/usr/bin/env python3
"""
Wall-temperature order study of docs/documentation/ib_accuracy_limits.md: the annulus at 10, 20 and 40 cells per radius with
inverse-distance (ib_ip_interp = 1) and TSC (2) image points, the observed orders of the wall error between successive grids, and
a check against the page's rows.

    python3 examples/2D_ibm_accuracy_limits/order_study.py [--cells 10 20 40] [--tol 0.1] [-- <./mfc.sh run options>]

Each run goes through ./mfc.sh run in this directory (default options: -n 1) and is measured with
../2D_ibm_conduction_annulus/check_annulus.py. Exits 1 if an observed order differs from the page by more than --tol.
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
CHECK = os.path.join(ROOT, "examples", "2D_ibm_conduction_annulus", "check_annulus.py")

# Observed orders (10-20, 20-40 cells/R) on the page, per image-point weighting and wall
PAGE = {1: {"inner": (0.6, 1.2), "outer": (0.6, 1.2)}, 2: {"inner": (2.3, 2.0), "outer": (1.1, 1.0)}}
NAME = {1: "inverse distance", 2: "TSC"}

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--cells", type=int, nargs=3, default=[10, 20, 40], help="Cells per radius of the three grids.")
parser.add_argument("--tol", type=float, default=0.1, help="Largest accepted difference from the page's orders.")
parser.add_argument("mfc_args", nargs="*", help="Options for ./mfc.sh run (after --).")
args = parser.parse_args()
mfc_args = args.mfc_args or ["-n", "1"]


def wall_errors(cells, interp):
    for d in ("D", "restart_data", "p_all", "silo_hdf5"):
        shutil.rmtree(os.path.join(HERE, d), ignore_errors=True)
    cmd = [os.path.join(ROOT, "mfc.sh"), "run", os.path.join(HERE, "case.py"), "-t", "pre_process", "simulation", *mfc_args]
    subprocess.run([*cmd, "--", "--cells_per_R", str(cells), "--interp", str(interp)], cwd=ROOT, check=True)
    out = subprocess.run([sys.executable, CHECK, HERE], check=True, capture_output=True, text=True).stdout
    print(out.strip(), flush=True)
    match = re.search(r"wall error inner (\S+) outer (\S+)", out)
    return abs(float(match.group(1))), abs(float(match.group(2)))


ok = True
for interp in (1, 2):
    err = [wall_errors(c, interp) for c in args.cells]
    print(f"\n{NAME[interp]} (ib_ip_interp = {interp})")
    print("| cells / R | inner | outer |\n|---:|---:|---:|")
    for c, (e_in, e_out) in zip(args.cells, err):
        print(f"| {c} | {e_in:.2e} | {e_out:.2e} |")
    for w, wall in enumerate(("inner", "outer")):
        orders = [math.log(err[k][w] / err[k + 1][w]) / math.log(args.cells[k + 1] / args.cells[k]) for k in (0, 1)]
        good = all(abs(o - e) <= args.tol for o, e in zip(orders, PAGE[interp][wall]))
        ok = ok and good
        print(f"{wall} wall order {orders[0]:.2f}, {orders[1]:.2f} (page {PAGE[interp][wall][0]}, {PAGE[interp][wall][1]}): {'matches' if good else 'DIFFERS'}")
sys.exit(0 if ok else 1)
