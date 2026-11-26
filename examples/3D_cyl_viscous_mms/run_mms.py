#!/usr/bin/env python3
"""Run and check examples/3D_cyl_viscous_mms (see README.md).

    python3 examples/3D_cyl_viscous_mms/run_mms.py --case A3 -N 32 64 128 [-- extra ./mfc.sh run flags]

For each N: pre_process, overwrite restart_data/lustre_0.dat with the manufactured state, one RK3 step with and
without viscosity. Rate = (q_visc - q_inviscid)/dt (the inviscid terms cancel), compared at cell centres with
exact_rates.RATES over the interior rings (the axis rings, the outer extrapolated rings and, for B3, r < 0.5
where the Couette field is singular, are excluded). Errors are max norms relative to the largest exact momentum
(energy) rate; B3, whose exact force is zero, is normalised by the pre-fix hoop error.

S2 / S3 (species diffusion, h2o2, unity Lewis, at rest at uniform T and p): Y_H2 = Y0 + eps (r/R)^2 (times
cos(theta) in 3D), with N2 compensating. Rate = (q with diffusion - q without)/dt of rho Y_H2; c(r) is the rate
(3D: its cos(theta) coefficient) over eps/R^2, against its exact value with rho D = k/cp from Cantera (see
species_coefficient). The error, relative to the exact mean, includes the difference between Cantera's and MFC's
transport fits.
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
        out = subprocess.run([MFC, "run", os.path.join(rd, "case.py"), "-t", t, "-n", "1", *flags], capture_output=True, text=True, check=False)
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
    nv = os.path.getsize(os.path.join(rd, "restart_data", "lustre_0.dat")) // (8 * NP * N * NX)
    return d3, NP, nv, RR, TH


def species_state(RR, TH, d3):
    """Cantera state of the species case per (theta, r) cell: rho, E, Y; rho D = k/cp(f) as a function of the H2 shift
    f; the index of H2; the base-state energy."""
    import cantera as ct
    from case import SPECIES

    gas = ct.Solution(SPECIES["mech"])
    gas.TPY = SPECIES["T"], SPECIES["p"], SPECIES["Y0"]
    ih2, in2 = gas.species_index("H2"), gas.species_index("N2")
    y0, E0 = gas.Y.copy(), gas.density * gas.int_energy_mass

    def at(f):
        y = y0.copy()
        y[ih2] += f
        y[in2] -= f
        gas.TPY = SPECIES["T"], SPECIES["p"], y
        return gas

    def rho_d(f):
        g = at(f)
        return g.thermal_conductivity / g.cp_mass

    f = SPECIES["eps"] * (RR / SPECIES["R"]) ** 2 * (np.cos(TH) if d3 else 1.0)
    rho, E, Y = np.empty(RR.shape), np.empty(RR.shape), np.empty(RR.shape + (gas.n_species,))
    for idx in np.ndindex(RR.shape):
        g = at(f[idx])
        rho[idx], E[idx], Y[idx] = g.density, g.density * g.int_energy_mass, g.Y
    return rho, E, Y, rho_d, ih2, E0


def write_ic(case, N, rd):
    d3, NP, nv, RR, TH = grid(case, N, rd)
    f0 = os.path.join(rd, "restart_data", "lustre_0.dat")
    q = np.fromfile(f0).reshape(nv, NP, N, NX)
    if case.startswith("S"):
        nd = 3 if d3 else 2
        rho, E, Y, _, _, E0 = species_state(RR, TH, d3)
        assert nv == 3 + nd + Y.shape[-1], f"unexpected restart layout: {nv} variables"
        # pre_process wrote the uniform base state: Cantera's energy must match MFC's for the state to be at uniform T
        if abs(q[1 + nd].mean() / E0 - 1) > 1e-8:
            sys.exit(f"MFC and Cantera energies differ: {q[1 + nd].mean():.12e} vs {E0:.12e}")
        q[0] = rho[:, :, None]
        q[1 : 1 + nd] = 0.0
        q[1 + nd] = E[:, :, None]
        q[2 + nd] = 1.0
        for k in range(Y.shape[-1]):
            q[3 + nd + k] = (rho * Y[..., k])[:, :, None]
        q.tofile(f0)
        return
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


def species_coefficient(case, N):
    """c(r) over the interior rings and its exact value: in 2D (1/r) d(r rho D dY/dr)/dr / (eps/R^2) = 4 rho D + 2 r d(rho D)/dr,
    with rho D(r) = k/cp of the local composition; in 3D the cos(theta) coefficient is 3 rho D of the base state (the
    variation of rho D contributes only to the 0 and 2 theta modes, at O(eps^2))."""
    rd = os.path.join(HERE, "runs", f"{case}_{N}_T")
    d3, NP, nv, RR, TH = grid(case, N, rd)
    from case import SPECIES

    _, _, _, rho_d, ih2, _ = species_state(RR, TH, d3)
    row = 3 + (3 if d3 else 2) + ih2

    def ld(v):
        return np.fromfile(os.path.join(HERE, "runs", f"{case}_{N}_{v}", "restart_data", "lustre_1.dat")).reshape(nv, NP, N, NX)[row].mean(axis=2)

    rate = (ld("T") - ld("F")) / SPECIES["dt"]
    eps, R = SPECIES["eps"], SPECIES["R"]
    c = (2.0 * (rate * np.cos(TH)).mean(axis=0) if d3 else rate[0]) / (eps / R**2)
    kb = max(3, N // 16)
    r = RR[0, kb : N - kb]
    if d3:
        cex = np.full_like(r, 3.0 * rho_d(0.0))
    else:
        h = 1.0e-3 * R
        rd_r = np.array([rho_d(eps * (x / R) ** 2) for x in r])
        drd = np.array([(rho_d(eps * ((x + h) / R) ** 2) - rho_d(eps * ((x - h) / R) ** 2)) / (2 * h) for x in r])
        cex = 4.0 * rd_r + 2.0 * r * drd
    return c[kb : N - kb], cex


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--case", default="A3", choices=[*RATES, "S2", "S3"])
    p.add_argument("-N", type=int, nargs="+", default=[32, 64, 128])
    p.add_argument("--no-run", action="store_true", help="only evaluate existing runs")
    a, flags = p.parse_known_args()
    flags = [f for f in flags if f != "--"]
    for N in a.N:
        if not a.no_run:
            for v in ("T", "F"):
                run_one(a.case, N, v, flags)
        if a.case.startswith("S"):
            c, cex = species_coefficient(a.case, N)
            print(f"{a.case} N={N}: rho Y_H2 rate, max |c - c_exact| / mean c_exact = {np.abs(c - cex).max() / cex.mean():.3e} (c mean {c.mean():.5e}, exact {cex.mean():.5e})")
            continue
        res = errors(a.case, N)
        print(f"{a.case} N={N}: " + "; ".join(f"{c} err {e:.3e} (from pre-fix model {m:.3e})" for c, (e, m) in res.items()))
