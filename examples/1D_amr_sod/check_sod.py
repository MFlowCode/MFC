#!/usr/bin/env python3
"""
L1 density error and shock and contact positions of 1D_amr_sod runs against the exact Riemann solution.

    check_sod.py RUN_DIR [RUN_DIR ...] [--max-l1 X]

Each RUN_DIR holds a finished run of case.py (parallel_io = F, one rank). The density is the composite solution: on each
part of the domain the finest level covering it, from the restart files (level 0 in q_cons_vf*.dat, the blocks in
amr_fine.dat). L1 = sum over those cells of |rho - exact cell average| * dx. The shock and the contact are where rho
crosses the middle of its exact jump, interpolated between cell centres. Exits 1 if any L1 exceeds --max-l1.
"""

import argparse
import os
import re

import numpy as np

from mfc.viz import read_amr_restart

GAMMA = 1.4
LEFT, RIGHT, X0 = (1.0, 0.0, 1.0), (0.125, 0.0, 0.1), 0.5


def star_state(left, right, g=GAMMA):
    """Pressure and velocity between the waves (Toro, ch. 4), by Newton on the pressure function."""

    def f(p, rho, pk):
        a = np.sqrt(g * pk / rho)
        if p > pk:  # shock
            A, B = 2 / ((g + 1) * rho), (g - 1) / (g + 1) * pk
            return (p - pk) * np.sqrt(A / (p + B)), np.sqrt(A / (p + B)) * (1 - (p - pk) / (2 * (B + p)))
        return 2 * a / (g - 1) * ((p / pk) ** ((g - 1) / (2 * g)) - 1), (p / pk) ** (-(g + 1) / (2 * g)) / (rho * a)

    (rl, ul, pl), (rr, ur, pr) = left, right
    p = 0.5 * (pl + pr)
    for _ in range(50):
        fl, dl = f(p, rl, pl)
        fr, dr = f(p, rr, pr)
        p = max(p - (fl + fr + ur - ul) / (dl + dr), 1e-12)
    fl, _ = f(p, rl, pl)
    fr, _ = f(p, rr, pr)
    return p, 0.5 * (ul + ur + fr - fl)


def exact(x, t, g=GAMMA):
    """Exact density of the Sod problem (left rarefaction, contact, right shock) and the contact and shock positions."""
    (rl, ul, pl), (rr, _, pr) = LEFT, RIGHT
    ps, us = star_state(LEFT, RIGHT)
    al = np.sqrt(g * pl / rl)
    rsl = rl * (ps / pl) ** (1 / g)
    asl = al * (ps / pl) ** ((g - 1) / (2 * g))
    rsr = rr * (ps / pr + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * ps / pr + 1)
    s = us * rsr / (rsr - rr)  # shock speed from mass conservation
    xi = (x - X0) / t
    fan = rl * np.maximum(2 / (g + 1) + (g - 1) / ((g + 1) * al) * (ul - xi), 0) ** (2 / (g - 1))
    rho = np.select([xi < ul - al, xi < us - asl, xi < us, xi < s], [rl, fan, rsl, rsr], rr)
    return rho, X0 + us * t, X0 + s * t, (rsl, rsr)


def read_inp(run_dir):
    with open(os.path.join(run_dir, "simulation.inp")) as f:
        txt = f.read()
    get = lambda k, d=None: (re.search(rf"^\s*{re.escape(k)}\s*=\s*([^\s,]+)", txt, re.M | re.I) or [None, d])[1]
    return {k: get(k, d) for k, d in (("m", None), ("dt", None), ("t_step_stop", None), ("amr", "F"))}


def composite(run_dir, inp):
    """Leaf cells (centre, width, rho) of the composite grid."""
    m, step = int(inp["m"]), int(inp["t_step_stop"])
    dx0 = 1.0 / (m + 1)
    if "T" in inp["amr"].upper():
        q0, blocks = read_amr_restart(run_dir, step)
    else:  # one Fortran record: 4-byte marker, then the data
        q0, blocks = [np.fromfile(os.path.join(run_dir, "p_all", "p0", str(step), "q_cons_vf1.dat"), "<f8", m + 1, offset=4)], []
    cells = {0: [((np.arange(m + 1) + 0.5) * dx0, dx0, q0[0])]}
    edges = {}  # level -> list of (x_lo, x_hi) of its blocks
    for b in blocks:
        nf = b.q.shape[1]
        h = dx0 / 2**b.level
        x_lo = b.lo[0] * h * nf / (b.hi[0] - b.lo[0] + 1)  # region in level-0 or parent-level cells
        edges.setdefault(b.level, []).append((x_lo, x_lo + nf * h))
        cells.setdefault(b.level, []).append((x_lo + (np.arange(nf) + 0.5) * h, h, b.q[0, :, 0, 0]))
    xs, ws, rs = [], [], []
    for lvl, lst in cells.items():
        for x, h, rho in lst:
            keep = np.ones(x.size, bool)
            for a, c in edges.get(lvl + 1, []):
                keep &= ~((x > a) & (x < c))
            xs.append(x[keep])
            ws.append(np.full(keep.sum(), h))
            rs.append(rho[keep])
    x, w, r = (np.concatenate(v) for v in (xs, ws, rs))
    o = np.argsort(x)
    assert abs(w.sum() - 1.0) < 1e-12, f"composite cells cover {w.sum()} of the domain"
    return x[o], w[o], r[o]


def crossing(x, r, level, lo, hi):
    """Rightmost crossing of rho = level in lo < x < hi, linear between cell centres."""
    i = [k for k in range(x.size - 1) if lo < x[k] < hi and (r[k] - level) * (r[k + 1] - level) <= 0][-1]
    return x[i] + (level - r[i]) * (x[i + 1] - x[i]) / (r[i + 1] - r[i])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--max-l1", type=float, default=None)
    a = ap.parse_args()
    ok = True
    print(f"{'run':40s} {'cells':>6s} {'dx_min':>9s} {'L1(rho)':>10s} {'contact err':>12s} {'shock err':>10s}")
    for d in a.run_dirs:
        inp = read_inp(d)
        t = float(inp["dt"]) * int(inp["t_step_stop"])
        x, w, r = composite(d, inp)
        sub = (np.arange(64) + 0.5) / 64 - 0.5
        avg = exact(x[:, None] + w[:, None] * sub[None, :], t)[0].mean(axis=1)
        _, xc, xs, (rsl, rsr) = exact(np.array([0.0]), t)
        l1 = float(np.sum(np.abs(r - avg) * w))
        ec = crossing(x, r, 0.5 * (rsl + rsr), xc - 0.1, 0.5 * (xc + xs)) - xc
        es = crossing(x, r, 0.5 * (rsr + RIGHT[0]), 0.5 * (xc + xs), xs + 0.1) - xs
        print(f"{os.path.basename(os.path.normpath(d)):40s} {x.size:6d} {w.min():9.3e} {l1:10.3e} {ec:+12.2e} {es:+10.2e}")
        ok &= a.max_l1 is None or l1 <= a.max_l1
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
