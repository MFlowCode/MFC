#!/usr/bin/env python3
"""Mean drag coefficient and Strouhal number of the cylinder wake from D/ib_forces.dat.

usage: check_wake.py <case dir> [<reference case dir>] [--window T]

Cd is the mean of Fx/(rho U^2 D / 2) and St = D/(U T) with T the mean period between upward zero
crossings of Cl - <Cl>, both over the whole shedding periods in the last T time units (default 250,
about 8 periods). A second case directory (e.g. the uniform fine grid) is compared against as well.
"""

import argparse
import os

import numpy as np

RHO, U, D = 1.0, 0.2, 1.0  # free stream of case.py
REFS = {  # Re = 100
    "Cd": [("Henderson 1995", 1.35), ("Park et al. 1998", 1.33)],
    "St": [("Williamson 1989 fit", -3.3265 / 100 + 0.1816 + 1.6e-4 * 100), ("Park et al. 1998", 0.165)],
}


def wake(case_dir, window):
    f = np.loadtxt(os.path.join(case_dir, "D", "ib_forces.dat"))
    t, cd, cl = f[:, 1], f[:, 2] / (0.5 * RHO * U**2 * D), f[:, 3] / (0.5 * RHO * U**2 * D)
    w = t >= t[-1] - window
    t, cd, cl = t[w], cd[w], cl[w]
    s = cl - cl.mean()
    i = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
    if len(i) < 3:
        raise SystemExit(f"{case_dir}: fewer than two shedding periods in the last {window} time units")
    tc = t[i] - s[i] * (t[i + 1] - t[i]) / (s[i + 1] - s[i])
    per = (tc[-1] - tc[0]) / (len(tc) - 1)
    k = (t >= tc[0]) & (t <= tc[-1])
    return {"Cd": cd[k].mean(), "St": D / (U * per), "Cl_amp": 0.5 * (cl[k].max() - cl[k].min()), "periods": len(tc) - 1, "spread": np.diff(tc).std() / per}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("ref_dir", nargs="?")
    ap.add_argument("--window", type=float, default=250.0)
    a = ap.parse_args()
    r = wake(a.case_dir, a.window)
    print(f"{a.case_dir}: {r['periods']} periods (period spread {100 * r['spread']:.2f} %), Cl amplitude {r['Cl_amp']:.4f}")
    other = wake(a.ref_dir, a.window) if a.ref_dir else None
    for q in ("Cd", "St"):
        line = f"  {q} = {r[q]:.4f}"
        for name, v in REFS[q]:
            line += f"; {100 * (r[q] - v) / v:+.2f} % vs {name} ({v:.4f})"
        if other:
            line += f"; {100 * (r[q] - other[q]) / other[q]:+.2f} % vs {a.ref_dir} ({other[q]:.4f})"
        print(line)


if __name__ == "__main__":
    main()
