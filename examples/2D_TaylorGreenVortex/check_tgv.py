#!/usr/bin/env python3
"""
Kinetic-energy decay of the 2D Taylor-Green vortex against the exact viscous rate.

    check_tgv.py RUN_DIR [RUN_DIR ...] [--tol TOL]

Each directory is a finished run of case.py (restart_data/lustre_<step>.dat, or p_all/ of a one-rank run). The exact solution of the
incompressible Navier-Stokes equations has u = U0 sin(x/L0) cos(y/L0) exp(-2 nu t / L0^2), so the kinetic energy decays as
exp(-4 nu t / L0^2), nu = mu / rho0. Per run the script prints the decay rate fitted to ln KE over all saves, its error against
4 nu / L0^2, and the error of KE at the last save against the exact value. Exit 1 if the rate error of the last run exceeds --tol.
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


def saves(d, nvar, shape):
    """{step: conservative fields 1..nvar} of a one-rank or parallel_io run (restart_data/lustre_<step>.dat or p_all/p0/<step>/)."""
    n = int(np.prod(shape))
    out = {}
    for f in glob.glob(os.path.join(d, "restart_data", "lustre_*.dat")):
        if re.search(r"lustre_\d+\.dat$", f):
            out[int(re.findall(r"\d+", os.path.basename(f))[0])] = np.fromfile(f, dtype=np.float64)[: nvar * n].reshape(nvar, *shape)
    if out:
        return out
    if len(glob.glob(os.path.join(d, "p_all", "p*"))) != 1:
        sys.exit(f"{d}: run on one rank or with parallel_io = T")
    for s in glob.glob(os.path.join(d, "p_all", "p0", "[0-9]*")):
        q = []
        for i in range(1, nvar + 1):
            with open(os.path.join(s, f"q_cons_vf{i}.dat"), "rb") as f:
                f.read(4)
                q.append(np.fromfile(f, dtype=np.float64, count=n).reshape(shape))
        out[int(os.path.basename(s))] = np.array(q)
    return out


def analyze(d):
    inp = read_inp(d)
    nx, ny = int(inp["m"]) + 1, int(inp["n"]) + 1
    data = saves(d, 3, (ny, nx))
    t, ke = [], []
    for s in sorted(data):
        q = data[s]
        t.append(s * float(inp["dt"]))
        ke.append(0.5 * np.mean((q[1] ** 2 + q[2] ** 2) / q[0]))
    t, ke = np.array(t), np.array(ke)
    rho0, mu, l0 = float(inp["patch_icpp(1)%alpha_rho(1)"]), 1.0 / float(inp["fluid_pp(1)%Re(1)"]), float(inp["patch_icpp(1)%vel(2)"])
    rate_exact = 4.0 * mu / rho0 / l0**2
    rate = -np.polyfit(t, np.log(ke), 1)[0]
    ke_err = ke[-1] / (ke[0] * np.exp(-rate_exact * (t[-1] - t[0]))) - 1.0
    return dict(name=os.path.basename(os.path.abspath(d)), n=nx, saves=len(t), decay=rate_exact * (t[-1] - t[0]), rate=rate, exact=rate_exact, ke_err=ke_err)


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("runs", nargs="+", help="Run directories, coarse to fine.")
parser.add_argument("--tol", type=float, default=None, help="Exit 1 if |rate / exact - 1| of the last run exceeds this.")
args = parser.parse_args()

res = [analyze(d) for d in args.runs]
print(f"{'run':>20} {'N':>5} {'saves':>5} {'4 nu t':>7} {'rate':>11} {'exact':>11} {'rate err':>10} {'KE(t) err':>10} {'order':>6}")
for i, r in enumerate(res):
    order = ""
    if i > 0:
        order = f"{np.log(abs(res[i - 1]['rate'] / res[i - 1]['exact'] - 1.0) / abs(r['rate'] / r['exact'] - 1.0)) / np.log(r['n'] / res[i - 1]['n']):.2f}"
    print(f"{r['name']:>20} {r['n']:>5} {r['saves']:>5} {r['decay']:>7.3f} {r['rate']:>11.5e} {r['exact']:>11.5e} {r['rate'] / r['exact'] - 1:>+10.2e} {r['ke_err']:>+10.2e} {order:>6}")
if args.tol is not None:
    err = abs(res[-1]["rate"] / res[-1]["exact"] - 1.0)
    print(f"{'PASS' if err <= args.tol else 'FAIL'}: decay-rate error {err:.2e} (tolerance {args.tol:g})")
    sys.exit(0 if err <= args.tol else 1)
