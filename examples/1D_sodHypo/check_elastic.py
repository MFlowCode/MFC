#!/usr/bin/env python3
"""
Hypoelastic Riemann problems (1D_sodHypo, 2D_hypo_shear_contact) against linear elastodynamics.

    check_elastic.py RUN_DIR [RUN_DIR ...] [--case-args ARGS] [--tol TOL]

Each RUN_DIR is a finished run (pre_process, simulation, post_process); --case-args are the arguments its case.py
was run with (';' separates one entry per run). Patch 1 is the left state and patch 2 the right state, both at
rest in x, with a jump in the normal stress sigma = -p + tau_xx and/or the shear stress tau_xy at the patch
boundary x0. Linear elastodynamics splits the jump into
  longitudinal waves at c_L = sqrt((K + 4G/3) / rho), K = gamma (p + pi_inf), impedance Z_L = rho c_L:
      u* = (sigma_R - sigma_L) / (Z_L,l + Z_L,r);  on each side dp = -K/M dsigma, dtau_xx = (4G/3)/M dsigma, M = K + 4G/3
  shear waves at c_S = sqrt(G / rho), impedance Z_S = rho c_S:
      v* - v0 = (tau_xy,R - tau_xy,L) / (Z_S,l + Z_S,r),  tau_xy* = tau_xy,L + Z_S,l (v* - v0).
The script averages the last save over y, measures the plateau values (median over the middle 60 % of each
plateau) and the wave fronts (where u or v - v0 crosses half its star value), and prints their errors relative
to the jump and to the wave speeds. These are first-order-in-amplitude predictions: the errors of a finite jump
include O(jump / M) nonlinear terms. Exit 1 if any speed or plateau error exceeds --tol.
"""

import argparse
import json
import math
import os
import subprocess
import sys

import numpy as np

from mfc.viz import assemble, assemble_silo, discover_format, discover_timesteps


def state(case, i):
    g = f"patch_icpp({i})%"
    nf = case["num_fluids"]
    rho = sum(case.get(f"{g}alpha_rho({q})", 0.0) for q in range(1, nf + 1))
    al = [case.get(f"{g}alpha({q})", 0.0) for q in range(1, nf + 1)]
    G = sum(a * case.get(f"fluid_pp({q})%G", 0.0) for q, a in zip(range(1, nf + 1), al))
    gm = sum(a * case[f"fluid_pp({q})%gamma"] for q, a in zip(range(1, nf + 1), al))  # 1/(gamma-1) of the mixture
    pim = sum(a * case.get(f"fluid_pp({q})%pi_inf", 0.0) for q, a in zip(range(1, nf + 1), al))
    p = case[f"{g}pres"]
    gam = 1 + 1 / gm
    K = gam * (p + pim / (gam * gm))
    M = K + 4 * G / 3
    return {
        "rho": rho,
        "p": p,
        "txx": case.get(f"{g}tau_e(1)", 0.0),
        "txy": case.get(f"{g}tau_e(2)", 0.0) if case["n"] > 0 else 0.0,
        "v0": case.get(f"{g}vel(2)", 0.0),
        "K": K,
        "M": M,
        "G": G,
        "cL": math.sqrt(M / rho),
        "cS": math.sqrt(G / rho),
    }


def front(x, y, target, xpred):
    """x where y crosses target, the crossing nearest xpred (linear interpolation)."""
    s = np.sign(y - target)
    idx = np.nonzero(s[1:] != s[:-1])[0]
    if idx.size == 0:
        return math.nan
    i = idx[np.argmin(np.abs(x[idx] - xpred))]
    return x[i] + (target - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])


def plateau(x, y, a, b):
    lo, hi = a + 0.2 * (b - a), b - 0.2 * (b - a)
    return float(np.median(y[(x > lo) & (x < hi)]))


def check(d, case_args):
    case = json.loads(subprocess.run([sys.executable, "case.py", *case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)
    L, R = state(case, 1), state(case, 2)
    x0 = case["patch_icpp(2)%x_centroid"] - 0.5 * case["patch_icpp(2)%length_x"]
    fmt = discover_format(d)
    step = discover_timesteps(d, fmt)[-1]
    t = step * case["dt"]
    a = assemble_silo(d, step) if fmt == "silo" else assemble(d, step, fmt)
    v = {k: (np.asarray(f).mean(axis=1) if np.ndim(f) == 2 else np.asarray(f)) for k, f in a.variables.items()}
    x = a.x_cc
    sL, sR = -L["p"] + L["txx"], -R["p"] + R["txx"]
    us = (sR - sL) / (L["rho"] * L["cL"] + R["rho"] * R["cL"])
    ss = sL + L["rho"] * L["cL"] * us
    xl, xr = x0 - L["cL"] * t, x0 + R["cL"] * t
    print(f"{d}: N = {x.size}, t = {t:.4e}, c_L = {L['cL']:.5g} | {R['cL']:.5g}, c_S = {L['cS']:.5g} | {R['cS']:.5g}")
    errs = []
    if abs(sR - sL) > 0:
        rows = [("u*", plateau(x, v["vel1"], xl, xr), us, abs(us))]
        for side, S, a_, b_ in (("left", L, xl, x0), ("right", R, x0, xr)):
            ds = ss - (sL if side == "left" else sR)
            rows.append((f"p* {side}", plateau(x, v["pres"], a_, b_), S["p"] - S["K"] / S["M"] * ds, abs(ds)))
            rows.append((f"tau_xx* {side}", plateau(x, v["tau1"], a_, b_), S["txx"] + 4 * S["G"] / 3 / S["M"] * ds, abs(ds)))
        for name, m, e, scale in rows:
            errs.append((m - e) / scale)
            print(f"  {name:14s} MFC {m:+.6e}  linear {e:+.6e}  error / jump {(m - e) / scale:+.2e}")
        for side, S, xp in (("left", L, xl), ("right", R, xr)):
            xf = front(x, v["vel1"], 0.5 * us, xp)
            c = abs(xf - x0) / t
            errs.append(c / S["cL"] - 1)
            print(f"  longitudinal front {side}: x = {xf:.5f} (linear {xp:.5f}), speed {c:.5g} vs c_L {S['cL']:.5g}: {c / S['cL'] - 1:+.2e}")
    dtxy = R["txy"] - L["txy"]
    if dtxy != 0:
        zl, zr = L["rho"] * L["cS"], R["rho"] * R["cS"]
        dv = dtxy / (zl + zr)
        tstar = L["txy"] + zl * dv
        sl, sr = x0 - L["cS"] * t, x0 + R["cS"] * t
        for name, m, e in (("v* - v0", plateau(x, v["vel2"], sl, sr) - L["v0"], dv), ("tau_xy*", plateau(x, v["tau2"], sl, sr), tstar)):
            errs.append((m - e) / abs(e))
            print(f"  {name:14s} MFC {m:+.6e}  linear {e:+.6e}  rel. error {errs[-1]:+.2e}")
        for side, S, xp in (("left", L, sl), ("right", R, sr)):
            xf = front(x, v["vel2"] - L["v0"], 0.5 * dv, xp)
            c = abs(xf - x0) / t
            errs.append(c / S["cS"] - 1)
            print(f"  shear front {side}: x = {xf:.5f} (linear {xp:.5f}), speed {c:.5g} vs c_S {S['cS']:.5g}: {c / S['cS'] - 1:+.2e}")
    return max(abs(e) for e in errs)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--case-args", default="", help="case.py arguments of the runs; ';' separates one entry per run")
    ap.add_argument("--tol", type=float, default=0.02, help="Max |error| of any speed or plateau value")
    args = ap.parse_args()
    cargs = args.case_args.split(";")
    if len(cargs) == 1:
        cargs *= len(args.runs)
    worst = max(check(os.path.abspath(d), c) for d, c in zip(args.runs, cargs))
    if not worst <= args.tol:
        print(f"FAIL: worst error {worst:.2e} > {args.tol}")
        sys.exit(1)
    print(f"PASS (worst error {worst:.2e})")


if __name__ == "__main__":
    main()
