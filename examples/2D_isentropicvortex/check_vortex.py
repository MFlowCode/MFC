#!/usr/bin/env python3
"""
Error of the translating isentropic vortex (case.py --convect) against the exact solution, and the observed order.

    check_vortex.py RUN_DIR [RUN_DIR ...] [--min-order P]

Each directory is a finished one-rank run of case.py --convect. The exact solution is the initial vortex (eps from the case,
default 5) moved by (u0, v0) t in the periodic box; the script integrates its density over each cell with a 5 x 5 Gauss rule and
prints the L1, L2 and max-norm density errors (L1 and L2 normalized by the box area) and the order between successive runs.
Exit 1 if the L1 order between the last two runs is below --min-order.
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


def exact_density(xc, yc, edges_x, edges_y, eps, box):
    xg, wg = np.polynomial.legendre.leggauss(5)
    avg = 0.0
    for a, wa in zip(xg, wg):
        x = 0.5 * (edges_x[:-1] + edges_x[1:]) + 0.5 * a * np.diff(edges_x)
        for b, wb in zip(xg, wg):
            y = 0.5 * (edges_y[:-1] + edges_y[1:]) + 0.5 * b * np.diff(edges_y)
            dxv = (x[None, :] - xc + 0.5 * box) % box - 0.5 * box  # nearest periodic image
            dyv = (y[:, None] - yc + 0.5 * box) % box - 0.5 * box
            t_fac = 1.0 - eps**2 / (16.0 * 2.4 * np.pi**2) * np.exp(2.0 * (1.0 - dxv**2 - dyv**2))
            avg = avg + 0.25 * wa * wb * t_fac**1.4
    return avg


def analyze(d):
    inp = read_inp(d)
    nx, ny = int(inp["m"]) + 1, int(inp["n"]) + 1
    x0, x1, y0, y1 = (float(inp[k]) for k in ("x_domain%beg", "x_domain%end", "y_domain%beg", "y_domain%end"))
    steps = [s for s in glob.glob(os.path.join(d, "p_all", "p0", "*")) if os.path.basename(s).isdigit()]
    if len(glob.glob(os.path.join(d, "p_all", "p*"))) != 1:
        sys.exit(f"{d}: run on one rank")
    last = max(steps, key=lambda s: int(os.path.basename(s)))
    step = int(os.path.basename(last))
    with open(os.path.join(last, "q_cons_vf1.dat"), "rb") as f:
        f.read(4)
        rho = np.fromfile(f, dtype=np.float64, count=nx * ny).reshape(ny, nx)
    t = step * float(inp["dt"])
    eps = float(inp.get("patch_icpp(1)%epsilon", 5.0))
    eps = eps if eps > 0 else 5.0
    xc = float(inp["patch_icpp(1)%x_centroid"]) + float(inp["patch_icpp(1)%vel(1)"]) * t
    yc = float(inp["patch_icpp(1)%y_centroid"]) + float(inp["patch_icpp(1)%vel(2)"]) * t
    ex = exact_density(xc, yc, np.linspace(x0, x1, nx + 1), np.linspace(y0, y1, ny + 1), eps, x1 - x0)
    e = np.abs(rho - ex)
    return dict(name=os.path.basename(os.path.abspath(d)), n=nx, t=t, steps=step, err=[e.mean(), np.sqrt((e**2).mean()), e.max()])


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("runs", nargs="+", help="Run directories, coarse to fine.")
parser.add_argument("--min-order", type=float, default=None, help="Exit 1 if the last L1 order is below this.")
args = parser.parse_args()

res = [analyze(d) for d in args.runs]
print(f"{'run':>20} {'N':>5} {'steps':>6} {'t':>6} {'L1':>10} {'L2':>10} {'Linf':>10} {'o(L1)':>6} {'o(Linf)':>7}")
order = None
for i, r in enumerate(res):
    o1 = oi = ""
    if i > 0 and r["n"] != res[i - 1]["n"]:
        ratio = np.log(r["n"] / res[i - 1]["n"])
        order = np.log(res[i - 1]["err"][0] / r["err"][0]) / ratio
        o1, oi = f"{order:.2f}", f"{np.log(res[i - 1]['err'][2] / r['err'][2]) / ratio:.2f}"
    print(f"{r['name']:>20} {r['n']:>5} {r['steps']:>6} {r['t']:>6.3f} {r['err'][0]:>10.3e} {r['err'][1]:>10.3e} {r['err'][2]:>10.3e} {o1:>6} {oi:>7}")
if args.min_order is not None and order is not None:
    print(f"{'PASS' if order >= args.min_order else 'FAIL'}: L1 order {order:.2f} (minimum {args.min_order:g})")
    sys.exit(0 if order >= args.min_order else 1)
