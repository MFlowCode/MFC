#!/usr/bin/env python3
"""
Temperature profiles of 1D_reactive_shocktube against Fig. 4.9 of Martinez Ferrer (2013), the case of Martinez Ferrer et al. (2014).

    compare_reference.py CASE_DIR [--case-args '--t_save_us 5']

Run the case with `--t_save_us 5` (saves every 5 us; the default saves are 2.3 us apart, too coarse to compare front positions) and
post_process first. At t = 170, 190 and 230 us it prints
  x_s   leading shock, the rightmost crossing of T = 975 K (748 K ahead, 1194 K induction plateau behind)
  x_r   reaction front, the rightmost crossing of T = 1600 K
  L1    (1/L) int |T_MFC - T_ref| dx as a percentage of the 2400 K axis range of the reference panel
  T_w   temperature at the wall cell, and the mean MFC - reference temperature in the burnt gas behind the merged front
against the curves in ref_T.txt.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from mfc.viz import assemble_silo, discover_timesteps

L = 0.12
TIMES = (170, 190, 230)
BURNT = {190: (0.01, 0.05), 230: (0.02, 0.08)}  # burnt-gas windows [m]
XF = np.linspace(0, L, 240001)


def crossing(x, y, level):
    """Rightmost x where y crosses level (linear interpolation)."""
    s = np.sign(y - level)
    i = np.where(s[:-1] != s[1:])[0]
    if len(i) == 0:
        return np.nan
    i = i[-1]
    return x[i] + (level - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--case-args", default="--t_save_us 5", help="Arguments given to case.py for the run (gives dt).")
    args = ap.parse_args()
    d = Path(args.case_dir)
    dt = json.loads(subprocess.run([sys.executable, "case.py", *args.case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)["dt"]
    ref = np.genfromtxt(Path(__file__).with_name("ref_T.txt"), delimiter=",", skip_header=4)
    steps = {round(s * dt * 1e6): s for s in discover_timesteps(str(d), "silo")}
    print(f"{'t [us]':>7} {'x_s [mm]':>9} {'ref':>7} {'x_r [mm]':>9} {'ref':>7} {'L1 T [%]':>9} {'T_w [K]':>8} {'ref':>7} {'burnt dT [K]':>12}")
    for t in TIMES:
        if t not in steps:
            print(f"{t:7d} not saved (run with --t_save_us 5)")
            continue
        a = assemble_silo(str(d), steps[t])
        x, T = a.x_cc, a.variables["T"]
        r = ref[ref[:, 0] == t]
        xr, Tr = r[:, 1] * 1e-3, r[:, 2]
        l1 = 100 * np.trapezoid(np.abs(np.interp(XF, x, T) - np.interp(XF, xr, Tr)), XF) / L / 2400.0
        burnt = ""
        if t in BURNT:
            xs = np.linspace(*BURNT[t], 401)
            burnt = f"{np.mean(np.interp(xs, x, T) - np.interp(xs, xr, Tr)):+12.1f}"
        print(
            f"{t:7d} {1e3 * crossing(x, T, 975.0):9.2f} {1e3 * crossing(xr, Tr, 975.0):7.2f} {1e3 * crossing(x, T, 1600.0):9.2f} {1e3 * crossing(xr, Tr, 1600.0):7.2f} {l1:9.2f} {T[0]:8.0f} {Tr[0]:7.0f} {burnt}"
        )


if __name__ == "__main__":
    main()
