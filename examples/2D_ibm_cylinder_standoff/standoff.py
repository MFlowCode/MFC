#!/usr/bin/env python3
"""
Bow-shock stand-off of a run of 2D_ibm_cylinder_standoff, against Billig's correlation.

    standoff.py RUN_DIR [--gamma 1.4]

Reads pre_process.inp, simulation.inp and the saves (restart_data/lustre_<step>.dat with parallel_io = T, or p_all/p0/<step>/
q_cons_vf<k>.dat without; one rank). On the stagnation line (the cell row through the body's centre) it finds, scanning inward from
the inflow, the first place where the pressure reaches the mean of the free-stream pressure p1 and the post-shock pressure p2 of
a normal shock, and prints for every save

  delta/R   distance of that point from the body surface, in body radii
  p_nose    pressure of the last cell before the body, over p2 (the Rayleigh pitot value over p2 is printed for reference)

and, for the mean of the last three saves, Billig's correlation (1967) delta/R = 0.143 exp(3.24/M^2) (sphere) or 0.386 exp(4.67/M^2) (cylinder) with
the signed error. Billig is an inviscid, high-Reynolds-number fit; the converged body-fitted value is a separate check.
"""

import argparse
import glob
import math
import os
import re

import numpy as np


def read_inp(path):
    out = {}
    for line in open(path):
        m = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
        if m:
            out[m.group(1)] = m.group(2).strip("'\"")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--gamma", type=float, default=1.4)
    args = ap.parse_args()
    g = args.gamma

    inp = {**read_inp(os.path.join(args.run, "pre_process.inp")), **read_inp(os.path.join(args.run, "simulation.inp"))}
    f = lambda k: float(inp[k].replace("d", "e"))  # noqa: E731
    m, n, p = int(inp["m"]) + 1, int(inp["n"]) + 1, int(inp["p"]) + 1
    dim = 3 if p > 1 else 2
    x0, dx = f("x_domain%beg"), (f("x_domain%end") - f("x_domain%beg")) / m
    mach, R = f("patch_icpp(1)%vel(1)"), f("patch_ib(1)%radius")
    xc = f("patch_ib(1)%x_centroid")
    shape = (p, n, m) if dim == 3 else (n, m)
    idx = [int((f(f"patch_ib(1)%{c}_centroid") - f(f"{c}_domain%beg")) / dx) for c in ("z", "y")[-(dim - 1) :]]

    steps = sorted({int(s.group(1)) for q in glob.glob(os.path.join(args.run, "restart_data", "lustre_*.dat")) if (s := re.search(r"lustre_(\d+)\.dat$", q))})
    if not steps:
        steps = sorted(int(os.path.basename(d)) for d in glob.glob(os.path.join(args.run, "p_all", "p0", "[0-9]*")))

    def load(step, k):
        path = os.path.join(args.run, "restart_data", f"lustre_{step}.dat")
        if os.path.exists(path):
            return np.fromfile(path, dtype="<f8").reshape(-1, *shape)[k]
        raw = np.fromfile(os.path.join(args.run, "p_all", "p0", str(step), f"q_cons_vf{k + 1}.dat"), dtype=np.uint8)
        return np.frombuffer(raw[4:-4].tobytes(), dtype="<f8").reshape(shape)

    p1 = 1.0 / g
    p2 = p1 * (1.0 + 2.0 * g / (g + 1.0) * (mach**2 - 1.0))
    pitot = p1 * (((g + 1.0) ** 2 * mach**2 / (4.0 * g * mach**2 - 2.0 * (g - 1.0))) ** (g / (g - 1.0))) * (1.0 - g + 2.0 * g * mach**2) / (g + 1.0)
    pmid = 0.5 * (p1 + p2)
    x = x0 + (np.arange(m) + 0.5) * dx
    surface = xc - R
    print(f"M {mach:g}  {'sphere' if dim == 3 else 'cylinder'}  {2 * R / dx:.1f} cells/D  p2/p1 {p2 / p1:.3f}  pitot/p2 {pitot / p2:.3f}")
    deltas = []
    for step in steps[1:]:
        rho, e = load(step, 0), load(step, dim + 1)
        mom2 = sum(load(step, 1 + d) ** 2 for d in range(dim))
        press = (g - 1.0) * (e - 0.5 * mom2 / rho)
        line = press[tuple(idx)] if dim == 3 else press[idx[0]]
        up = np.flatnonzero(x < surface)
        j = next((i for i in up if line[i] >= pmid), None)
        if j is None or j == 0:
            print(f"  step {step}: no shock found")
            continue
        xs = x[j - 1] + (pmid - line[j - 1]) / (line[j] - line[j - 1]) * dx
        delta = (surface - xs) / R
        nose = line[up[-1]] / p2
        deltas.append(delta)
        print(f"  step {step:>7d}  delta/R {delta:.4f}  p_nose/p2 {nose:.3f}")
    ref = 0.143 * math.exp(3.24 / mach**2) if dim == 3 else 0.386 * math.exp(4.67 / mach**2)
    late = np.mean(deltas[-3:])
    print(f"mean of the last three saves {late:.4f} (range {max(deltas[-3:]) - min(deltas[-3:]):.4f}); Billig delta/R {ref:.4f}: {100 * (late / ref - 1):+.1f} %")


if __name__ == "__main__":
    main()
