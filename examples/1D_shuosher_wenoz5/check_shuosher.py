#!/usr/bin/env python3
"""
Density error of Shu-Osher runs against a fine reference run.

    check_shuosher.py REF_DIR RUN_DIR [RUN_DIR ...] [--tol TOL]

Every directory is a finished 1D run of the same problem (any examples/1D_shuosher_* case) whose last save holds the
conservative fields (restart_data/lustre_<step>.dat or p_all/p0/<step>/). The reference (case.py -N 16016 here) is averaged
onto each run's cells, so its cell count must be a multiple of theirs. Prints the L1 density error (sum |rho - rho_ref| dx) over
the domain and over the shocked entropy waves, 5.5 <= x - x_beg <= 7.5, and the order between successive runs of one scheme.
Exit 1 if the domain L1 error of the last run exceeds --tol.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np


def read_inp(d):
    out = {}
    for name in ("pre_process.inp", "simulation.inp"):
        for line in open(os.path.join(d, name)):
            match = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
            if match:
                out[match.group(1)] = match.group(2).strip("'\",")
    return out


def last_density(d):
    inp = read_inp(d)
    n = int(inp["m"]) + 1
    lustre = [f for f in glob.glob(os.path.join(d, "restart_data", "lustre_*.dat")) if re.search(r"lustre_\d+\.dat$", f)]
    if lustre:
        rho = np.fromfile(max(lustre, key=lambda f: int(re.findall(r"\d+", os.path.basename(f))[0])), dtype=np.float64)[:n]
    else:
        last = max(glob.glob(os.path.join(d, "p_all", "p0", "[0-9]*")), key=lambda s: int(os.path.basename(s)))
        with open(os.path.join(last, "q_cons_vf1.dat"), "rb") as f:
            f.read(4)
            rho = np.fromfile(f, dtype=np.float64, count=n)
    return rho, float(inp["x_domain%beg"]), float(inp["x_domain%end"]), inp


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("ref", help="Reference run directory (fine grid).")
parser.add_argument("runs", nargs="+", help="Run directories.")
parser.add_argument("--tol", type=float, default=None, help="Exit 1 if the domain L1 error of the last run exceeds this.")
args = parser.parse_args()

ref, x0, x1, _ = last_density(args.ref)
print(f"reference {os.path.basename(os.path.abspath(args.ref))}: {ref.size} cells on [{x0:g}, {x1:g}]")
print(f"{'run':>24} {'cells':>6} {'L1 domain':>10} {'L1 waves':>10} {'order':>6}")
prev = None
for d in args.runs:
    rho, b0, b1, _ = last_density(d)
    if (b0, b1) != (x0, x1) or ref.size % rho.size:
        sys.exit(f"{d}: domain or cell count ({rho.size}) does not fit the reference")
    avg = ref.reshape(rho.size, -1).mean(axis=1)
    dx = (x1 - x0) / rho.size
    xc = x0 + (np.arange(rho.size) + 0.5) * dx - x0
    err = np.abs(rho - avg) * dx
    l1, l1w = err.sum(), err[(xc >= 5.5) & (xc <= 7.5)].sum()
    order = ""
    if prev is not None and prev[0] != rho.size:
        order = f"{np.log(prev[1] / l1) / np.log(rho.size / prev[0]):.2f}"
    print(f"{os.path.basename(os.path.abspath(d)):>24} {rho.size:>6} {l1:>10.3e} {l1w:>10.3e} {order:>6}")
    prev = (rho.size, l1)
if args.tol is not None:
    print(f"{'PASS' if prev[1] <= args.tol else 'FAIL'}: L1 error {prev[1]:.3e} (tolerance {args.tol:g})")
    sys.exit(0 if prev[1] <= args.tol else 1)
