#!/usr/bin/env python3
"""Run and check examples/3D_cyl_viscous_mms (see README.md).

    python3 examples/3D_cyl_viscous_mms/run_mms.py --case A3 -N 32 64 128 [-- extra ./mfc.sh run flags]

For each N: pre_process, overwrite restart_data/lustre_0.dat with the manufactured state, one RK3 step with and
without viscosity. Rate = (q_visc - q_inviscid)/dt (the inviscid terms cancel), compared at cell centres with
exact_rates.RATES over the interior rings (the axis rings, the outer extrapolated rings and, for B3, r < 0.5
where the Couette field is singular, are excluded). Errors are max norms relative to the largest exact momentum
(energy) rate; B3, whose exact force is zero, is normalised by the pre-fix hoop error.
"""

import argparse
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from exact_rates import FIELDS, MASTER, RATES  # noqa: E402

MFC = os.path.join(HERE, "..", "..", "mfc.sh")
A, NX, R, DT, GAMP = 0.1, 26, 2.0, 1.0e-8, 1.0 / 0.4


def run_one(case, N, visc, flags):
    rd = os.path.join(HERE, "runs", f"{case}_{N}_{visc}")
    os.makedirs(rd, exist_ok=True)
    write_case(case, N, visc, rd)
    for t in ("pre_process", "simulation"):
        if t == "simulation":
            write_ic(case, N, rd)
        out = subprocess.run([MFC, "run", os.path.join(rd, "case.py"), "-t", t, "-n", "1", *flags], capture_output=True, text=True)
        if out.returncode:
            sys.exit(f"{t} failed in {rd}:\n{out.stdout[-4000:]}\n{out.stderr[-2000:]}")


def write_case(case, N, visc, rd):
    out = subprocess.run([sys.executable, os.path.join(HERE, "case.py"), "--case", case, "-N", str(N), "--viscous", visc], check=True, capture_output=True, text=True)
    open(os.path.join(rd, "case.py"), "w").write(f"print({out.stdout.strip()!r})\n")


def grid(case, N, rd):
    d3 = case.endswith("3")
    NP = N if d3 else 1
    yb = np.fromfile(os.path.join(rd, "restart_data", "lustre_y_cb.dat"))  # the 2D axisymmetric grid is not uniform at the axis
    rc = 0.5 * (yb[1:] + yb[:-1])
    thc = (np.arange(NP) + 0.5) * 2 * np.pi / NP if d3 else np.zeros(1)
    TH, RR = np.meshgrid(thc, rc, indexing="ij")
    return d3, NP, (6 if d3 else 5), RR, TH


def write_ic(case, N, rd):
    d3, NP, nv, RR, TH = grid(case, N, rd)
    f0 = os.path.join(rd, "restart_data", "lustre_0.dat")
    q = np.fromfile(f0).reshape(nv, NP, N, NX)
    vel = FIELDS[case](RR, TH, A)[: 3 if d3 else 2]
    q[0] = 1.0
    for i, v in enumerate(vel):
        q[1 + i] = v[:, :, None]
    q[1 + len(vel)] = (GAMP + 0.5 * sum(v * v for v in vel))[:, :, None]
    q[-1] = 1.0
    q.tofile(f0)


def errors(case, N):
    d3, NP, nv, RR, TH = grid(case, N, os.path.join(HERE, "runs", f"{case}_{N}_T"))
    mub = 1.0 if case in ("C3", "E2") else 0.0

    def ld(v):
        return np.fromfile(os.path.join(HERE, "runs", f"{case}_{N}_{v}", "restart_data", "lustre_1.dat")).reshape(nv, NP, N, NX)

    rate = ((ld("T") - ld("F")) / DT).mean(axis=3)
    ex = [np.asarray(e) for e in RATES[case](RR, TH, A, 1.0, mub)]
    ms = [np.asarray(e) for e in MASTER[case](RR, TH, A, 1.0, mub)]
    kb = max(3, N // 16)
    k0 = max(kb, int(np.searchsorted(RR[0], 0.5))) if case == "B3" else kb
    sl = (slice(None), slice(k0, N - kb))
    rows = {"x": 1, "r": 2, "th": 3, "E": 4} if d3 else {"x": 1, "r": 2, "E": 3}
    comps = {"x": 0, "r": 1, "th": 2, "E": 3}
    moms = [c for c in rows if c != "E"]
    smom = max(np.abs(ex[comps[c]][sl]).max() for c in moms) or max(np.abs(ms[comps[c]][sl]).max() for c in moms)
    res = {}
    for c, row in rows.items():
        scale = np.abs(ex[3][sl]).max() if c == "E" else smom
        res[c] = (np.abs(rate[row][sl] - ex[comps[c]][sl]).max() / scale, np.abs(rate[row][sl] - ms[comps[c]][sl]).max() / scale)
    return res


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--case", default="A3", choices=list(RATES))
    p.add_argument("-N", type=int, nargs="+", default=[32, 64, 128])
    p.add_argument("--no-run", action="store_true", help="only evaluate existing runs")
    a, flags = p.parse_known_args()
    flags = [f for f in flags if f != "--"]
    for N in a.N:
        if not a.no_run:
            for v in ("T", "F"):
                run_one(a.case, N, v, flags)
        res = errors(a.case, N)
        print(f"{a.case} N={N}: " + "; ".join(f"{c} err {e:.3e} (from pre-fix model {m:.3e})" for c, (e, m) in res.items()))
