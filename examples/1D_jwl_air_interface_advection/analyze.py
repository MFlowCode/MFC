#!/usr/bin/env python3
"""Report the Abgrall interface-advection metrics for the JWL/air closure.

The case initializes one uniform pressure p0 and one uniform velocity u0 across
a JWL-products slab embedded in air, so the exact solution holds p == p0 and
u == u0 for all time while the two material contacts merely translate. This
script reads the primitive-variable output and reports, for each saved step,

    max|p - p0| / p0     and     max|u - u0|

split into the pure regions (products volume fraction alpha < 0.01 or > 0.99)
and the WENO-reconstructed mixed interface band (0.01 <= alpha <= 0.99). p0 and
u0 are the uniform reference, taken from the median of the pure-region field at
each step (robust to the few transition cells).

Note on step 0: MFC's pre_process writes its t=0 primitive diagnostic by
recovering pressure with the products mass fraction defaulted to 1, so the air
region is written at omega_products/Gamma_air = 0.75x the true value. The
CONSERVED energy is correct (from step 1 on both regions recover the uniform
p0), so step 0 is reported separately and excluded from the pass metric.
"""

import glob
import os
import sys

import numpy as np

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "D")


def series(var):
    files = sorted(glob.glob(os.path.join(D, f"prim.{var}.*.dat")))
    if not files:
        sys.exit(f"no prim.{var}.*.dat files found in {D}; run the case first")
    return files


def col2(path):
    return np.loadtxt(path)[:, 1]


def main():
    pres_files = series(4)
    vel_files = series(3)
    alpha_files = series(5)

    print(f"{'step':>10} {'region':>7} {'p0 (Pa)':>14} " f"{'max|dp|/p0':>12} {'max|du| (m/s)':>14}")

    worst_pure_dp = worst_pure_du = 0.0
    for i, (pf, uf, af) in enumerate(zip(pres_files, vel_files, alpha_files)):
        p, u, a = col2(pf), col2(uf), col2(af)
        mixed = (a >= 0.01) & (a <= 0.99)
        pure = ~mixed
        p0 = float(np.median(p[pure]))
        u0 = float(np.median(u[pure]))
        step = os.path.basename(pf).split(".")[3]

        for label, mask in (("pure", pure), ("mixed", mixed)):
            if not np.any(mask):
                continue
            dp = float(np.max(np.abs(p[mask] - p0))) / abs(p0)
            du = float(np.max(np.abs(u[mask] - u0)))
            tag = " (t=0 diag)" if i == 0 else ""
            print(f"{step:>10} {label:>7} {p0:>14.6e} {dp:>12.3e} {du:>14.3e}{tag}")
            if label == "pure" and i > 0:
                worst_pure_dp = max(worst_pure_dp, dp)
                worst_pure_du = max(worst_pure_du, du)

    print()
    print(f"pure-region worst (steps > 0): max|dp|/p0 = {worst_pure_dp:.3e}, " f"max|du| = {worst_pure_du:.3e} m/s")


if __name__ == "__main__":
    main()
