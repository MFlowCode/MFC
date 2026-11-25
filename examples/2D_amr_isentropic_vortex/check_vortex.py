#!/usr/bin/env python3
"""
Density error of the isentropic vortex against the exact cell averages at the last save, on the composite grid (the fine cells
where a block covers a point, level 0 elsewhere), and the observed order between successive run directories of one kind.

    python3 check_vortex.py <run dir> [<run dir> ...]

Each run directory holds a case.py run at one --N. L1 is the area-weighted mean |error|, Linf the largest; "seam" is the
largest error in the cells within 3 level-0 cells of a block edge, "away" in all other cells. The vortex ends at x = 3.
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amr_composite import _num, leaves, read_inp  # noqa: E402

GAM, EPS, X0, U0 = 1.4, 5.0, -3.0, 1.0


def rho_exact(x, y, t):
    r2 = (x - X0 - U0 * t) ** 2 + y**2
    return (1.0 - (GAM - 1.0) * EPS**2 / (8.0 * GAM * math.pi**2) * np.exp(1.0 - r2)) ** (1.0 / (GAM - 1.0))


def rho_avg(x, y, dx, t):
    """Exact cell average (4x4-point Gauss, eighth order), the quantity a finite-volume cell holds."""
    g, w = np.polynomial.legendre.leggauss(4)
    return sum(wa * wb / 4.0 * rho_exact(x + a * dx[0] / 2, y + b * dx[1] / 2, t) for a, wa in zip(g, w) for b, wb in zip(g, w))


def kind(inp):
    """Run label and base N: a uniform run at 2N (--fine) takes the base-N step, twice the CFL number of one at N."""
    cells = int(inp["m"]) + 1
    if inp.get("amr", "F") != "T":
        fine = _num(inp["dt"]) * cells / (_num(inp["x_domain%end"]) - _num(inp["x_domain%beg"])) > 0.1
        return ("uniform 2N", cells // 2) if fine else ("uniform N", cells)
    mode = "regrid" if int(inp.get("amr_regrid_int", "0")) > 0 else "static"
    sub = "subcycle" if inp.get("amr_subcycle", "F") == "T" else "lock-step"
    return f"{mode} {sub} {'quartic' if inp.get('amr_prolong', '1') == '2' else 'linear'}", cells


def errors(run_dir):
    inp = read_inp(run_dir)
    step = int(inp["t_step_stop"])
    t = _num(inp["dt"]) * step
    parts, _ = leaves(run_dir, step)
    edges = [p[3] for p in parts[1:]]
    dx0 = parts[0][1]
    x0, y0 = _num(inp["x_domain%beg"]), _num(inp["y_domain%beg"])
    l1 = area = 0.0
    seam = away = 0.0
    for (x, y), dx, q, _ in parts:
        e = np.abs(q - rho_avg(x, y, dx, t))
        a = float(np.prod(dx))
        l1 += e.sum() * a
        area += e.size * a
        near = np.zeros(e.size, dtype=bool)
        for blo, bhi in edges:
            # distance outside (positive) or inside (negative) the block box, in level-0 cells
            ix = np.maximum(x0 + blo[0] * dx0[0] - x, x - x0 - (bhi[0] + 1) * dx0[0]) / dx0[0]
            iy = np.maximum(y0 + blo[1] * dx0[1] - y, y - y0 - (bhi[1] + 1) * dx0[1]) / dx0[1]
            d = np.maximum(ix, iy)
            near |= np.abs(d) < 3.0
        seam = max(seam, float(e[near].max(initial=0.0)))
        away = max(away, float(e[~near].max(initial=0.0)))
    fine = sum(p[2].size for p in parts[1:]) / (4.0 * np.prod([int(inp["m"]) + 1, int(inp["n"]) + 1]))
    label, nbase = kind(inp)
    return {"kind": label, "N": nbase, "L1": l1 / area, "Linf": max(seam, away), "seam": seam, "away": away, "cover": fine}


def main(dirs):
    rows = [errors(d) for d in dirs]
    print(f"{'run':28s} {'N':>5s} {'cover':>6s} {'L1':>10s} {'order':>6s} {'Linf':>10s} {'order':>6s} {'seam':>10s} {'away':>10s}")
    prev = {}
    for r in rows:
        p = prev.get(r["kind"])
        o1 = oi = ""
        if p is not None and r["N"] != p["N"]:
            f = math.log(r["N"] / p["N"])
            o1 = f"{math.log(p['L1'] / r['L1']) / f:6.2f}"
            oi = f"{math.log(p['Linf'] / r['Linf']) / f:6.2f}"
        print(f"{r['kind']:28s} {r['N']:5d} {r['cover']:6.3f} {r['L1']:10.3e} {o1:>6s} {r['Linf']:10.3e} {oi:>6s} {r['seam']:10.3e} {r['away']:10.3e}")
        prev[r["kind"]] = r


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
