#!/usr/bin/env python3
"""
L1 error of a 1D two-state Riemann problem (Sod, Lax, 123, ...) against the exact solution.

    check_riemann.py RUN_DIR [RUN_DIR ...] [--tol TOL]

Each directory is a finished single-fluid ideal-gas run (1D, or a 2D strip uniform in y) with two patches (left and right state) whose last save holds the
conservative fields (p_all/p0/<step>/ or restart_data/lustre_<step>.dat). The states, the interface, gamma and the time are read
from pre_process.inp and simulation.inp, so the script also checks examples/1D_sodshocktube and examples/1D_laxshocktube.
Prints the L1 norms (sum |q - q_exact| dx over the domain) of density, velocity and pressure against exact cell averages and, for
several runs, the observed order between successive grids. Exit 1 if the density L1 error of the finest run exceeds --tol.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
from riemann_exact import cell_averages


def read_inp(d):
    out = {}
    for name in ("pre_process.inp", "simulation.inp"):
        for line in open(os.path.join(d, name)):
            match = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
            if match:
                out[match.group(1)] = match.group(2).strip("'\",")
    return out


def last_save(d, nvar, ncell):
    """Conservative fields (nvar, ncell) and the step of the last save of a 1D run (or of a quasi-1D strip, ncell = nx * ny)."""
    lustre = [f for f in glob.glob(os.path.join(d, "restart_data", "lustre_*.dat")) if re.search(r"lustre_\d+\.dat$", f)]
    if lustre:
        path = max(lustre, key=lambda f: int(re.findall(r"\d+", os.path.basename(f))[0]))
        return np.fromfile(path, dtype=np.float64)[: nvar * ncell].reshape(nvar, ncell), int(re.findall(r"\d+", os.path.basename(path))[0])
    steps = [s for s in glob.glob(os.path.join(d, "p_all", "p0", "*")) if os.path.basename(s).isdigit()]
    if len(glob.glob(os.path.join(d, "p_all", "p*"))) != 1:
        sys.exit(f"{d}: run on one rank (or with parallel_io = T)")
    last = max(steps, key=lambda s: int(os.path.basename(s)))
    q = []
    for i in range(1, nvar + 1):
        with open(os.path.join(last, f"q_cons_vf{i}.dat"), "rb") as f:
            f.read(4)
            q.append(np.fromfile(f, dtype=np.float64, count=ncell))
    return np.array(q), int(os.path.basename(last))


def analyze(d):
    inp = read_inp(d)
    ncell, ny, x0, x1 = int(inp["m"]) + 1, int(inp.get("n", 0)) + 1, float(inp["x_domain%beg"]), float(inp["x_domain%end"])
    g = 1.0 + 1.0 / float(inp["fluid_pp(1)%gamma"])
    states = [tuple(float(inp[f"patch_icpp({k})%{v}"]) for v in ("alpha_rho(1)", "vel(1)", "pres")) for k in (1, 2)]
    xc1, len1 = float(inp["patch_icpp(1)%x_centroid"]), float(inp["patch_icpp(1)%length_x"])
    xc2 = float(inp["patch_icpp(2)%x_centroid"])
    xi = xc1 + 0.5 * len1 if xc1 < xc2 else xc1 - 0.5 * len1
    left, right = (states[0], states[1]) if xc1 < xc2 else (states[1], states[0])
    nvar = 3 if ny == 1 else 4  # a quasi-1D 2D strip (IGR) carries y momentum
    q, step = last_save(d, nvar, ncell * ny)
    if step == 0:
        sys.exit(f"{d}: no save after step 0 (did the run finish?)")
    q = q.reshape(nvar, ny, ncell).mean(axis=1)
    t = step * float(inp["dt"])
    rho = q[0]
    u = q[1] / rho
    p = (g - 1.0) * (q[nvar - 1] - 0.5 * rho * u**2)
    edges = np.linspace(x0, x1, ncell + 1)
    ex = cell_averages(left, right, g, xi, t, edges)
    dx = (x1 - x0) / ncell
    err = [float(np.sum(np.abs(a - b)) * dx) for a, b in zip((rho, u, p), ex)]
    return dict(name=os.path.basename(os.path.abspath(d)), n=ncell, t=t, err=err)


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("runs", nargs="+", help="Run directories, coarse to fine for an order estimate.")
parser.add_argument("--tol", type=float, default=None, help="Exit 1 if the density L1 error of the last run exceeds this.")
args = parser.parse_args()

res = [analyze(d) for d in args.runs]
print(f"{'run':>24} {'cells':>6} {'t':>7} {'L1 rho':>10} {'L1 u':>10} {'L1 p':>10} {'order rho':>9}")
for i, r in enumerate(res):
    order = ""
    if i > 0 and r["n"] != res[i - 1]["n"]:
        order = f"{np.log(res[i - 1]['err'][0] / r['err'][0]) / np.log(r['n'] / res[i - 1]['n']):.2f}"
    print(f"{r['name']:>24} {r['n']:>6} {r['t']:>7.4f} {r['err'][0]:>10.3e} {r['err'][1]:>10.3e} {r['err'][2]:>10.3e} {order:>9}")
if args.tol is not None:
    ok = res[-1]["err"][0] <= args.tol
    print(f"{'PASS' if ok else 'FAIL'}: density L1 error {res[-1]['err'][0]:.3e} (tolerance {args.tol:g})")
    sys.exit(0 if ok else 1)
