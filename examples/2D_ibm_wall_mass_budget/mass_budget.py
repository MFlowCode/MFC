#!/usr/bin/env python3
"""
Close the gas-mass budget of the 2D_ibm_wall_mass_budget case from its output.

Reads pre_process.inp and simulation.inp (grid, domain, cylinder, save interval), restart_data/lustre_<k>.dat (conservative
variables, parallel_io = T, one rank) and D/ib_surface_*_<k>.dat (ib_surface_wrt), and for every save k >= --skip reports, with t_1 the first save used:

    gain    = M(t_k) - M(t_1)          gas mass in the fluid cells (kg per unit depth)
    gasif   = integral of sum(area*mdot) over the body's surface records from t_1 to t_k (trapezoid rule)
    defect  = (gain - gasif)/gasif     0 when the wall passes exactly what the surface gasifies
    leak_i  = (M_i(t_k) - M_i(t_1))/gasif for each --inert species (0 when the wall passes none of it)

Without ib_flux_wall the records are a two-cell band of ghost points whose area sums to the wall area within about 1 %,
so the defect there carries that error; with it they are the face areas that set the face fluxes, and the defect is the
budget's closure (plus the trapezoid rule between saves).
"""

import argparse
import glob
import os
import re

import numpy as np

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("case_dir", help="Directory of the run (holds simulation.inp, D/ and restart_data/).")
parser.add_argument("--skip", type=int, default=1, help="First save used (the start-up step is excluded; default 1).")
parser.add_argument("--inert", nargs="*", default=["N2"], help="Species that must not cross the wall: names from the gas mechanism (needs cantera) or 0-based indices.")
parser.add_argument(
    "--mech",
    default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "2D_ibm_reacting_surface", "carbon_gasphase_reduced_gri11.yaml"),
    help="Gas mechanism, for the species names of --inert.",
)
parser.add_argument("--dtype", default="float64", help="Storage type of restart_data (float32 for single-precision builds).")
parser.add_argument("--tol", type=float, default=None, help="Exit 1 if the final |defect| exceeds this.")
args = parser.parse_args()


def read_inp(path):
    out = {}
    for line in open(path):
        match = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
        if match:
            out[match.group(1)] = match.group(2).strip("'\"")
    return out


inp = {**read_inp(os.path.join(args.case_dir, "pre_process.inp")), **read_inp(os.path.join(args.case_dir, "simulation.inp"))}
m, n = int(inp["m"]) + 1, int(inp["n"]) + 1
x0, x1 = float(inp.get("x_domain%beg", 0.0)), float(inp["x_domain%end"])
y0, y1 = float(inp.get("y_domain%beg", 0.0)), float(inp["y_domain%end"])
dx, dy = (x1 - x0) / m, (y1 - y0) / n
t_save = float(inp["t_save"])
xc, yc, rad = (float(inp[f"patch_ib(1)%{k}"]) for k in ("x_centroid", "y_centroid", "radius"))
X, Y = np.meshgrid(x0 + (np.arange(m) + 0.5) * dx, y0 + (np.arange(n) + 0.5) * dy, indexing="ij")
fluid = np.hypot(X - xc, Y - yc) > rad

names = None
if any(not s.isdigit() for s in args.inert):
    import cantera as ct

    names = ct.Solution(args.mech).species_names

saves = sorted(int(re.search(r"lustre_(\d+)\.dat$", f).group(1)) for f in glob.glob(os.path.join(args.case_dir, "restart_data", "lustre_*.dat")))
saves = [k for k in saves if k >= args.skip]
if len(saves) < 2:
    raise SystemExit(f"need at least two saves from {args.skip} in {args.case_dir}/restart_data")


def masses(k):
    """Total gas mass and per-species gas masses in the fluid cells, kg per unit depth."""
    raw = np.fromfile(os.path.join(args.case_dir, "restart_data", f"lustre_{k}.dat"), dtype=args.dtype)
    q = raw.reshape(-1, n, m)  # variable, y, x (x fastest)
    q = np.swapaxes(q, 1, 2)  # variable, x, y
    cell = dx * dy
    # variables: rho, rho*u, rho*v, E, alpha, then rho*Y_k
    return float(q[0][fluid].sum() * cell), np.array([float(s[fluid].sum() * cell) for s in q[5:]])


def rate(k):
    """Gasification rate of the body, sum(area * mdot) over its surface records, kg/s per unit depth."""
    total = 0.0
    for f in glob.glob(os.path.join(args.case_dir, "D", f"ib_surface_*_{k}.dat")):
        rec = np.loadtxt(f, ndmin=2)
        if rec.size:
            total += float((rec[:, 7] * rec[:, 9]).sum())
    return total


inert = [int(s) if s.isdigit() else names.index(s) for s in args.inert]
M0, Ms0 = masses(saves[0])
R_prev, gasif = rate(saves[0]), 0.0
print(f"{'t [s]':>12} {'gain [kg/m]':>14} {'gasif [kg/m]':>14} {'defect':>11} " + " ".join(f"leak_{s:>6}" for s in args.inert))
for k_prev, k in zip(saves[:-1], saves[1:]):
    R = rate(k)
    gasif += 0.5 * (R_prev + R) * (k - k_prev) * t_save
    R_prev = R
    M, Ms = masses(k)
    gain = M - M0
    line = f"{k * t_save:12.5e} {gain:14.6e} {gasif:14.6e} {(gain - gasif) / gasif:11.3e} "
    print(line + " ".join(f"{(Ms[i] - Ms0[i]) / gasif:12.3e}" for i in inert))
defect = (gain - gasif) / gasif
if args.tol is not None and abs(defect) > args.tol:
    raise SystemExit(f"|defect| = {abs(defect):.3e} > {args.tol:g}")
