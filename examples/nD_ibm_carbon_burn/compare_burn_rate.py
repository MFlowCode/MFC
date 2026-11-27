#!/usr/bin/env python3
"""
Burn rate and flame of nD_ibm_carbon_burn against the 1D radial reference (ref_cylinder.txt or ref_sphere.txt).

    compare_burn_rate.py CASE_DIR [--ndim 3] [--times 10 20 30 60] [--tol 0.05]

Burn rate: CASE_DIR/D/ib_heat.dat column 4 (face mass outflow) and column 3 (surface-solve gasification rate), averaged over +-0.5 us
around each time, against the reference 2 pi R sum W_k s_k (4 pi R^2 sum W_k s_k for the sphere); error = MFC/reference - 1.
Compare at times at least 0.5 us before the end of the run: the window is cut off at the last row and the burn rate falls steeply, so a final-time value reads high.
Flame (needs the Silo output of post_process; skipped without it): T_max of the temperature binned in r (bin width dx) and its radius,
at the last save, against the reference at the same time. Exits 1 if |column 4 error| > --tol at any time.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def inp(run, key):
    txt = open(os.path.join(run, "simulation.inp")).read()
    return float(re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", txt, re.M).group(1).replace("d", "e").replace("D", "e"))


def flame(run, ndim, R, ref):
    """T_max and its radius over R of the temperature binned in r, at the last save, with the reference at the same time."""
    from mfc.viz import assemble_silo, discover_timesteps

    s = discover_timesteps(run, "silo")[-1]
    t_us = s * inp(run, "t_save") * 1e6  # adaptive time step: save k is at k t_save
    a = assemble_silo(run, s)
    grids = np.meshgrid(*[a.x_cc, a.y_cc] + ([a.z_cc] if ndim == 3 else []), indexing="ij")
    r = np.sqrt(sum(g**2 for g in grids)).ravel() / R
    T = a.variables["T"].ravel()
    dx = (a.x_cc[1] - a.x_cc[0]) / R
    idx = np.digitize(r, np.arange(1.0, 5.0 + dx, dx))
    bins = [i for i in range(1, idx.max() + 1) if (idx == i).any()]
    rb, Tb = np.array([r[idx == i].mean() for i in bins]), np.array([T[idx == i].mean() for i in bins])
    k, j = np.argmax(Tb), np.argmin(np.abs(ref[:, 0] - t_us))
    return t_us, Tb[k], rb[k], ref[j, 2], ref[j, 3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--ndim", type=int, default=2, choices=(2, 3))
    ap.add_argument("--times", type=float, nargs="+", default=[10.0, 20.0, 30.0, 60.0], help="Times [us].")
    ap.add_argument("--tol", type=float, default=0.05)
    args = ap.parse_args()
    ref = np.loadtxt(os.path.join(HERE, "ref_cylinder.txt" if args.ndim == 2 else "ref_sphere.txt"))
    rows = [ln.split() for ln in open(glob.glob(os.path.join(args.case_dir, "D", "ib_heat.dat"))[0])]
    h = np.array([r for r in rows if len(r) == 6], dtype=float)  # a run in progress may end on a partial row
    t_us, c3, c4 = h[:, 1] * 1e6, h[:, 2], h[:, 3]
    worst = 0.0
    print("t [us]   reference     column 4     error     column 3     error")
    for t in args.times:
        w = np.abs(t_us - t) <= 0.5
        if not w.any():
            continue
        mr = np.interp(t, ref[:, 0], ref[:, 1])
        e4, e3 = c4[w].mean() / mr - 1, c3[w].mean() / mr - 1
        worst = max(worst, abs(e4))
        print(f"{t:6.1f}  {mr:.6e}  {c4[w].mean():.6e}  {e4:+8.3%}  {c3[w].mean():.6e}  {e3:+8.3%}")
    try:
        R = float(re.search(r"patch_ib\(1\)%radius\s*=\s*(\S+)", open(os.path.join(args.case_dir, "simulation.inp")).read()).group(1).replace("d", "e"))
        t_us, Tm, rm, Tr, rr = flame(args.case_dir, args.ndim, R, ref)
        print(f"flame at {t_us:.1f} us: T_max {Tm:.1f} K at r = {rm:.2f} R (reference: {Tr:.1f} K at {rr:.2f} R)")
    except (ImportError, IndexError, KeyError, OSError, ValueError, AttributeError) as err:
        print(f"flame comparison skipped ({type(err).__name__}: {err})")
    ok = worst <= args.tol
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
