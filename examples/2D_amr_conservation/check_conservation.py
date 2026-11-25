#!/usr/bin/env python3
"""
Drift of the total mass, momentum and energy of a 2D_amr_conservation run.

    check_conservation.py RUN_DIR [--tol 1e-13]

Reads the serial restart files of every save of RUN_DIR (parallel_io = F, one rank). The totals are level-0 sums: at the end
of each step every covered level-0 cell is set to the average of its children, so they equal the composite sums. The
checker confirms that on each save by summing the leaf cells of every level (blocks from amr_fine.dat, the uncovered level-0
cells). Sums are exact (math.fsum), so the drift is the solver's. Prints (S(t) - S(0))/S(0) per variable and save, the largest
relative gap between the level-0 and the composite totals, and exits 1 if any |drift| exceeds --tol.
"""

import argparse
import math
import os

import numpy as np

from mfc.viz import read_amr_restart

NAMES = ("mass", "x-momentum", "y-momentum", "energy")


def fsum(a):
    """Exact sum of each variable of a (nvar, ...) array."""
    return np.array([math.fsum(v.ravel()) for v in a])


def totals(run_dir, step, n0):
    """Level-0 and composite totals (per unit area) of the conserved variables, and the refined cell count."""
    q0, blocks = read_amr_restart(run_dir, step)
    q0 = q0.reshape((q0.shape[0], n0, n0), order="F")
    # per level: [(lower corner in that level's cells, fields)]
    levels = {0: [((0, 0), q0)]}
    for b in blocks:
        levels.setdefault(b.level, []).append(((b.lo[0] * 2**b.level, b.lo[1] * 2**b.level), b.q[..., 0]))
    comp = np.zeros(q0.shape[0])
    for lvl, patches in levels.items():
        nl = n0 * 2**lvl
        covered = np.zeros((nl, nl), bool)
        for (i, j), q in levels.get(lvl + 1, []):
            covered[i // 2 : (i + q.shape[1]) // 2, j // 2 : (j + q.shape[2]) // 2] = True
        for (i, j), q in patches:
            leaf = ~covered[i : i + q.shape[1], j : j + q.shape[2]]
            comp += fsum(q[:, leaf]) / nl**2
    return fsum(q0) / n0**2, comp, sum(q.shape[1] * q.shape[2] for lvl, p in levels.items() if lvl for _, q in p)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--tol", type=float, default=1e-13)
    a = ap.parse_args()
    p0 = os.path.join(a.run_dir, "p_all", "p0")
    steps = sorted(int(s) for s in os.listdir(p0) if s.isdigit())
    n0 = math.isqrt(int(np.fromfile(os.path.join(p0, "0", "q_cons_vf1.dat"), "<i4", 1)[0]) // 8)
    # step 0 is pre_process's level-0 data; the square's edges lie on base faces, so it is also the blocks' initial average
    ref = fsum([np.fromfile(os.path.join(p0, "0", f"q_cons_vf{i}.dat"), "<f8", n0 * n0, offset=4) for i in (1, 2, 3, 4)]) / n0**2
    worst = 0.0
    print(f"{'step':>6s} {'fine cells':>10s} " + " ".join(f"{n:>12s}" for n in NAMES) + f" {'lev0 - comp':>12s}")
    for s in steps[1:]:
        lev0, comp, nf = totals(a.run_dir, s, n0)
        drift = (lev0[:4] - ref) / ref
        worst = max(worst, float(np.max(np.abs(drift))))
        gap = float(np.max(np.abs(lev0[:4] - comp[:4]) / np.abs(ref)))
        print(f"{s:6d} {nf:10d} " + " ".join(f"{d:+12.2e}" for d in drift) + f" {gap:12.2e}")
    print(f"largest |drift| {worst:.2e} (tol {a.tol:.0e})")
    raise SystemExit(0 if worst <= a.tol else 1)


if __name__ == "__main__":
    main()
