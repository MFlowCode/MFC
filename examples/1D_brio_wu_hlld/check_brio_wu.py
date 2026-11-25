#!/usr/bin/env python3
"""
Brio-Wu MHD shock tube (1D_brio_wu, 1D_brio_wu_hlld): grid convergence and the slow compound wave.

    check_brio_wu.py --ref REF_DIR RUN_DIR [RUN_DIR ...] [--rj] [--tol TOL]

RUN_DIRs and REF_DIR are finished runs of the same case (same gamma and end time); the cell count of REF_DIR must
be an integer multiple of each RUN_DIR's, and REF_DIR is block-averaged onto each coarser grid. Prints the L1
error (mean |error|) of rho, u, v, p and B_y at the last save against REF_DIR and the observed order of the rho
error between successive runs. There is no exact solution: with the compound wave the Riemann problem has
non-unique solutions, and dissipative schemes converge to the one with the compound wave.

--rj (for runs with --gamma 5/3) also prints the states on both sides of the intermediate shock of the slow
compound wave and compares them with the values Ryu and Jones (1995, section 5) read from their 8192-cell
calculation: preshock rho 0.647, p 0.484, B_y 0.536; postshock rho 0.827, p 0.758. The shock is where B_y changes
sign; preshock = median over cells 12..4 before it, postshock = the rho maximum within 10 cells after it (the shock is
followed directly by the slow rarefaction, so the postshock state is a peak; their postshock B_y = -0.0738 sits on the
steep front of that rarefaction and is not compared).
Exit 1 if the finest run's L1(rho) exceeds --tol.
"""

import argparse
import math
import os
import sys

import numpy as np

from mfc.viz import assemble, assemble_silo, discover_format, discover_timesteps

VARS = ("rho", "vel1", "vel2", "pres", "By")
RJ = {"pre": {"rho": 0.647, "pres": 0.484, "By": 0.536}, "post": {"rho": 0.827, "pres": 0.758}}


def last(d):
    fmt = discover_format(d)
    step = discover_timesteps(d, fmt)[-1]
    a = assemble_silo(d, step) if fmt == "silo" else assemble(d, step, fmt)
    return a.x_cc, {k: np.asarray(a.variables[k]) for k in VARS}


def compound(x, v):
    """States on both sides of the intermediate shock (B_y sign change left of the contact)."""
    i0 = int(np.argmax(v["By"] < 0))
    pre = {k: float(np.median(v[k][i0 - 12 : i0 - 3])) for k in RJ["pre"]}
    j = i0 + int(np.argmax(v["rho"][i0 : i0 + 10]))
    post = {k: float(v[k][j]) for k in RJ["post"]}
    return x[i0], pre, post


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--rj", action="store_true", help="Compare the compound wave with Ryu and Jones (gamma = 5/3)")
    ap.add_argument("--tol", type=float, default=1e-2, help="Max L1(rho) of the finest run")
    args = ap.parse_args()
    xr, vr = last(os.path.abspath(args.ref))
    res = []
    for d in args.runs:
        x, v = last(os.path.abspath(d))
        k = xr.size // x.size
        if k * x.size != xr.size:
            sys.exit(f"{d}: {x.size} cells does not divide the reference's {xr.size}")
        l1 = {n: float(np.mean(np.abs(v[n] - vr[n].reshape(-1, k).mean(axis=1)))) for n in VARS}
        print(f"{d}: N = {x.size}, " + ", ".join(f"L1({n}) = {e:.3e}" for n, e in l1.items()))
        res.append((x.size, l1["rho"]))
    res.sort()
    for (n1, e1), (n2, e2) in ((a, b) for a, b in zip(res, res[1:]) if a[0] != b[0]):
        print(f"observed order of L1(rho) {n1} -> {n2}: {math.log(e1 / e2) / math.log(n2 / n1):.2f}")
    if args.rj:
        for name, (x, v) in [("ref", (xr, vr))] + [(d, last(os.path.abspath(d))) for d in args.runs]:
            xs, pre, post = compound(x, v)
            row = [f"{side} {n} {val:.4f} ({val / RJ[side][n] - 1:+.1%})" for side, st in (("pre", pre), ("post", post)) for n, val in st.items()]
            print(f"compound wave, {name} (N = {x.size}, shock at x = {xs:.4f}): " + ", ".join(row))
    if res[-1][1] > args.tol:
        print(f"FAIL: L1(rho) = {res[-1][1]:.3e} > {args.tol}")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
