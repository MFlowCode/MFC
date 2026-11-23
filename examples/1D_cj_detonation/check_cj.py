#!/usr/bin/env python3
"""
Front speed and post-detonation state of 1D_cj_detonation against the equilibrium CJ solution (cj_ref.json).

    check_cj.py CASE_DIR [--case-args '--ic tz --ppi 5 --L 0.5'] [--tol 0.02]

Reads the Silo output of CASE_DIR. Prints
  D            least-squares slope of the shock position over the second half of the run, against D_CJ; the shock is the rightmost
               crossing of p = p1 + (p_CJ - p1)/2, linearly interpolated between cell centres
  D_end        mean local speed (5-save fits) over the last 5 saves
  p_max/p_vN   mean peak pressure within 20 induction lengths of the shock, over the second half, against the vN pressure
  plateau      (--ic tz) mean p and T over the u = 0 plateau against the Taylor-Zeldovich plateau, at x < 0.9 c_plateau t_eff,
               t_eff = x_shock / D_CJ
Exits 1 if |D_end / D_CJ - 1| > --tol.
"""

import argparse
import json
import os
import subprocess
import sys

import numpy as np

from mfc.viz import assemble_silo, discover_timesteps

HERE = os.path.dirname(os.path.abspath(__file__))


def front(x, p, thr):
    i = np.nonzero(p > thr)[0][-1]
    return np.nan if i == p.size - 1 else x[i] + (thr - p[i]) / (p[i + 1] - p[i]) * (x[i + 1] - x[i])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--case-args", default="", help="Arguments given to case.py for the run (gives dt, the grid and the start).")
    ap.add_argument("--tol", type=float, default=0.02, help="Max |D_end / D_CJ - 1|.")
    args = ap.parse_args()
    cargs = args.case_args.split()
    dt = json.loads(subprocess.run([sys.executable, os.path.join(HERE, "case.py"), *cargs], cwd=args.case_dir, capture_output=True, text=True, check=True).stdout)["dt"]
    ref = json.load(open(os.path.join(HERE, "cj_ref.json")))
    thr = ref["p1"] + 0.5 * (ref["p_CJ"] - ref["p1"])
    t, xs, pmax = [], [], []
    for s in discover_timesteps(args.case_dir, "silo"):
        a = assemble_silo(args.case_dir, s)
        x, p = a.x_cc, a.variables["pres"]
        if p.max() <= thr:
            continue
        xf = front(x, p, thr)
        if np.isfinite(xf):
            t.append(s * dt)
            xs.append(xf)
            pmax.append(p[np.abs(x - xf) < 20 * ref["ind_len"]].max())
    t, xs, pmax = map(np.array, (t, xs, pmax))
    w = t >= 0.5 * t[-1]
    D = np.polyfit(t[w], xs[w], 1)[0]
    loc = np.array([np.polyfit(t[i - 2 : i + 3], xs[i - 2 : i + 3], 1)[0] for i in range(2, len(t) - 2)])
    D_end = loc[-5:].mean()
    print(f"{t.size} saves to t = {t[-1] * 1e6:.2f} us, shock at x = {xs[-1] * 1e3:.2f} mm; D_CJ = {ref['D_CJ']:.2f} m/s")
    print(f"D = {D:.2f} m/s, error {D / ref['D_CJ'] - 1:+.2%}; D_end = {D_end:.2f} m/s, error {D_end / ref['D_CJ'] - 1:+.2%}; local sd {loc[len(loc) // 2 :].std() / ref['D_CJ']:.2%}")
    print(f"p_max / p_vN = {pmax[w].mean() / ref['p_vN']:.3f}")
    if "tz" in cargs:
        a = assemble_silo(args.case_dir, discover_timesteps(args.case_dir, "silo")[-1])
        pl = ref["plateau"]
        xi = a.x_cc / (xs[-1] / ref["D_CJ"])
        sel = xi < 0.9 * pl["c"]
        p = a.variables["pres"][sel].mean()
        print(f"plateau p = {p / 1e3:.1f} kPa, error {p / pl['p'] - 1:+.2%} ({sel.sum()} cells)")
    ok = abs(D_end / ref["D_CJ"] - 1) <= args.tol
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
