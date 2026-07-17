#!/usr/bin/env python3
"""Grid-convergence reduction for the 3D octant PETN spherical ZND case.

Reads the radial probe histories (probe{i}_prim.dat, written by MFC's
probe_wrt) from one or more grid_<N>/D/ run directories, fits a front-arrival
speed D_MFC = dr/dt through the probes' shock-arrival times, and (with >= 2
resolutions) reports the observed order of grid convergence on that speed
against the calibrated CJ eigenvalue (P8 protocol: p = log2((f_c-f_m)/(f_m-f_f))
for 3 grids; a simple ratio-based estimate for 2).

Usage:
    ./analyze_convergence.py grid_39 grid_79 [grid_159 ...]
"""

import glob
import math
import os
import sys

import numpy as np

D_CJ = 6245.0  # measured, G-invariant CJ speed for PETN_KUHL (see case.py docstring)
P_AMB = 101325.0
ARRIVAL_MULT = 2.0  # shock arrival: pressure exceeds ARRIVAL_MULT * P_AMB
TIME_COL, PRES_COL = 0, 5

# Must match case.py's L, charge_radius, target_probe_r (not re-imported: case.py
# needs the mfc toolchain on sys.path, which isn't guaranteed from this script).
L = 0.04
CHARGE_RADIUS = 0.015
TARGET_PROBE_R = tuple(f * CHARGE_RADIUS for f in (0.3, 0.6, 1.0, 1.5, 2.0))


def grid_index(run_dir):
    return int("".join(filter(str.isdigit, os.path.basename(os.path.normpath(run_dir)))))


def probe_radius(run_dir):
    """Recompute each probe's actual (nearest-cell-center) radius at this grid's dx."""
    grid = grid_index(run_dir)
    dx = L / (grid + 1)
    yz = 0.5 * dx

    def nearest_cell_center(x):
        i = int(round(x / dx - 0.5))
        i = max(0, min(grid, i))
        return (i + 0.5) * dx

    return [math.sqrt(nearest_cell_center(r) ** 2 + yz**2 + yz**2) for r in TARGET_PROBE_R]


def arrival_time(path):
    rows = np.loadtxt(path)
    if rows.ndim == 1:
        rows = rows[None, :]
    t, p = rows[:, TIME_COL], rows[:, PRES_COL]
    hit = np.where(p > ARRIVAL_MULT * P_AMB)[0]
    return t[hit[0]] if hit.size else math.nan


def front_speed(run_dir):
    d = os.path.join(run_dir, "D")
    radii = probe_radius(run_dir)
    probe_files = sorted(glob.glob(os.path.join(d, "probe*_prim.dat")), key=lambda p: int("".join(filter(str.isdigit, os.path.basename(p)))))
    if len(probe_files) != len(radii):
        sys.exit(f"{run_dir}: found {len(probe_files)} probe files but case.py defines {len(radii)} probes")
    arrivals = [arrival_time(f) for f in probe_files]
    r, t = np.array(radii), np.array(arrivals)
    ok = ~np.isnan(t)
    if ok.sum() < 2:
        return math.nan, radii, arrivals
    D = np.polyfit(t[ok], r[ok], 1)[0]
    return D, radii, arrivals


def main():
    run_dirs = sys.argv[1:]
    if not run_dirs:
        sys.exit(__doc__)

    results = []
    for run_dir in run_dirs:
        D, radii, arrivals = front_speed(run_dir)
        results.append((run_dir, D))
        err = 100.0 * (D - D_CJ) / D_CJ
        print(f"{run_dir}: D_front = {D:.1f} m/s  ({err:+.2f}% vs D_CJ = {D_CJ:.0f} m/s)")
        for r, t in zip(radii, arrivals):
            print(f"    r = {r * 1e3:6.2f} mm   t_arrival = {t:.4e} s" if not math.isnan(t) else f"    r = {r * 1e3:6.2f} mm   no arrival")

    Ds = [D for _, D in results]
    if len(Ds) == 3 and all(not math.isnan(D) for D in Ds):
        # coarse, medium, fine ordering assumed from argv order (dx halves each step)
        f_c, f_m, f_f = Ds
        p_obs = math.log2(abs(f_c - f_m) / abs(f_m - f_f))
        print(f"\nobserved order of convergence (front speed): p = {p_obs:.2f}")
    elif len(Ds) == 2 and all(not math.isnan(D) for D in Ds):
        print(f"\ntwo-grid change in front speed: {100.0 * abs(Ds[1] - Ds[0]) / D_CJ:.2f}% of D_CJ")


if __name__ == "__main__":
    main()
