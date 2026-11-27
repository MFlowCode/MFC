#!/usr/bin/env python3
"""
Water-air shock tube (1D_kapilashocktube) against the exact two-material Riemann solution.

    check_riemann.py CASE_DIR [CASE_DIR ...] [--case-args '-N 1000' ...] [--tol TOL]

Each CASE_DIR is a finished run (pre_process, simulation, post_process; the Silo output is written next to case.py).
The left state is pure fluid 1 and the right state pure fluid 2, each a stiffened gas p = (gamma - 1) rho e - gamma pi_inf.
The exact solution (Toro, Riemann Solvers, ch. 4, with p replaced by p + pi_inf in each material) is sampled at the
cell centres at the time of the last save. Prints, per run, the star state, the L1 errors of rho, u and p (normalized by
the jump across the tube), the material-interface position from the volume fraction, the star-plateau errors of the finest run;
with several runs, the observed order between successive grids. Exits 1 if the L1 density error of the finest run
exceeds --tol.
"""

import argparse
import json
import math
import subprocess
import sys

import numpy as np

from mfc.viz import assemble_silo, discover_timesteps


def material(case, i):
    """gamma and pi_inf of fluid i from MFC's fluid_pp gamma = 1/(gamma - 1) and pi_inf = gamma pi_inf/(gamma - 1)."""
    g = 1.0 + 1.0 / case[f"fluid_pp({i})%gamma"]
    return g, case.get(f"fluid_pp({i})%pi_inf", 0.0) * (g - 1.0) / g


class Side:
    def __init__(self, rho, u, p, g, pi):
        self.rho, self.u, self.p, self.g, self.pi = rho, u, p, g, pi
        self.c = math.sqrt(g * (p + pi) / rho)

    def f(self, ps):
        """Velocity change across the wave to the star pressure ps, and its derivative."""
        g, pp, pk = self.g, ps + self.pi, self.p + self.pi
        if ps > self.p:  # shock
            a, b = 2.0 / ((g + 1.0) * self.rho), (g - 1.0) / (g + 1.0) * pk
            q = math.sqrt(a / (pp + b))
            return (ps - self.p) * q, q * (1.0 - 0.5 * (ps - self.p) / (pp + b))
        r = (pp / pk) ** ((g - 1.0) / (2.0 * g))
        return 2.0 * self.c / (g - 1.0) * (r - 1.0), r / (self.rho * self.c) * pk / pp

    def star_rho(self, ps):
        g, pp, pk = self.g, ps + self.pi, self.p + self.pi
        if ps > self.p:
            k = (g - 1.0) / (g + 1.0)
            return self.rho * (pp / pk + k) / (k * pp / pk + 1.0)
        return self.rho * (pp / pk) ** (1.0 / g)


def star(L, R):
    floor = -min(L.pi, R.pi)  # p + pi_inf > 0 in both materials
    ps = 0.5 * (L.p + R.p)
    for _ in range(200):
        fl, dl = L.f(ps)
        fr, dr = R.f(ps)
        new = ps - (fl + fr + R.u - L.u) / (dl + dr)
        new = 0.5 * (ps + floor) if new <= floor else new
        done = abs(new - ps) <= 1e-14 * (abs(ps) + abs(floor))
        ps = new
        if done:
            break
    return ps, 0.5 * (L.u + R.u) + 0.5 * (R.f(ps)[0] - L.f(ps)[0])


def sample(L, R, ps, us, xi):
    """rho, u, p at similarity coordinate xi = (x - x0)/t."""
    out = np.empty((3, xi.size))
    for j, s in enumerate(xi):
        left = s <= us
        K, sgn = (L, -1.0) if left else (R, 1.0)
        g, pk = K.g, K.p + K.pi
        if ps > K.p:  # shock
            ws = K.u + sgn * K.c * math.sqrt((g + 1.0) / (2.0 * g) * (ps + K.pi) / pk + (g - 1.0) / (2.0 * g))
            out[:, j] = (K.star_rho(ps), us, ps) if sgn * (s - ws) < 0 else (K.rho, K.u, K.p)
            continue
        rs = K.star_rho(ps)
        cs = math.sqrt(g * (ps + K.pi) / rs)
        head, tail = K.u + sgn * K.c, us + sgn * cs
        if sgn * (s - head) >= 0:
            out[:, j] = (K.rho, K.u, K.p)
        elif sgn * (s - tail) <= 0:
            out[:, j] = (rs, us, ps)
        else:  # inside the fan
            u = 2.0 / (g + 1.0) * (-sgn * K.c + 0.5 * (g - 1.0) * K.u + s)
            c = sgn * (s - u)
            rho = K.rho * (c / K.c) ** (2.0 / (g - 1.0))
            out[:, j] = (rho, u, pk * (c / K.c) ** (2.0 * g / (g - 1.0)) - K.pi)
    return out


def analyze(d, case_args):
    case = json.loads(subprocess.run([sys.executable, "case.py", *case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)
    gl, pil = material(case, 1)
    gr, pir = material(case, 2)
    p1, p2 = "patch_icpp(1)%", "patch_icpp(2)%"
    L = Side(case[p1 + "alpha_rho(1)"], case[p1 + "vel(1)"], case[p1 + "pres"], gl, pil)
    R = Side(case[p2 + "alpha_rho(2)"], case[p2 + "vel(1)"], case[p2 + "pres"], gr, pir)
    x0 = case[p2 + "x_centroid"] - 0.5 * case[p2 + "length_x"]
    steps = discover_timesteps(d, "silo")
    t = steps[-1] * case["dt"]
    data = assemble_silo(d, steps[-1])
    x = np.asarray(data.x_cc).ravel()
    V = {k: np.asarray(v).ravel() for k, v in data.variables.items()}
    ps, us = star(L, R)
    ex = sample(L, R, ps, us, (x - x0) / t)
    dx = x[1] - x[0]
    num = (V["alpha_rho1"] + V["alpha_rho2"], V["vel1"], V["pres"])
    scale = (abs(L.rho - R.rho), max(abs(us - L.u), abs(us - R.u)), abs(L.p - R.p))
    l1 = [float(np.sum(np.abs(n - e)) * dx / (s * (x[-1] - x[0] + dx))) for n, e, s in zip(num, ex, scale)]
    a1 = V["alpha1"]
    i = int(np.argmax(a1 < 0.5))
    x_int = x[i - 1] + (0.5 - a1[i - 1]) * dx / (a1[i] - a1[i - 1])
    # plateau values: the middle fifth of each star region (left: rarefaction tail to contact, right: contact to shock)
    tail = us - math.sqrt(L.g * (ps + L.pi) / L.star_rho(ps)) if ps < L.p else us
    shock = R.u + R.c * math.sqrt((R.g + 1.0) / (2.0 * R.g) * (ps + R.pi) / (R.p + R.pi) + (R.g - 1.0) / (2.0 * R.g)) if ps > R.p else us
    plateau = []
    for a, b, rho_s in ((tail, us, L.star_rho(ps)), (us, shock, R.star_rho(ps))):
        xi = (x - x0) / t
        m = np.abs(xi - 0.5 * (a + b)) <= 0.1 * (b - a)
        plateau.append([float(num[0][m].mean() / rho_s - 1.0), float(num[1][m].mean() / us - 1.0), float(num[2][m].mean() / ps - 1.0)])
    return dict(N=x.size, t=t, ps=ps, us=us, rl=L.star_rho(ps), rr=R.star_rho(ps), l1=l1, x_int=x_int, x_int_exact=x0 + us * t, dx=dx, plateau=plateau)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dirs", nargs="+")
    ap.add_argument("--case-args", nargs="*", default=None, help="case.py arguments of each run (default: none).")
    ap.add_argument("--tol", type=float, default=0.0, help="Fail if the L1 density error of the finest run exceeds this (0: no check).")
    args = ap.parse_args()
    cargs = args.case_args or [""] * len(args.case_dirs)
    res = sorted((analyze(d, c) for d, c in zip(args.case_dirs, cargs)), key=lambda r: r["N"])
    r = res[0]
    print(f"exact star state: p* = {r['ps']:.6e} Pa, u* = {r['us']:.4f} m/s, rho*L = {r['rl']:.4f}, rho*R = {r['rr']:.4f} kg/m^3, t = {r['t']:.4e} s")
    print(f"{'N':>6} {'L1 rho':>10} {'L1 u':>10} {'L1 p':>10} {'interface - exact [cells]':>26}")
    for k, r in enumerate(res):
        order = ""
        if k:
            q = res[k - 1]
            order = "  order " + " ".join(f"{math.log(q['l1'][i] / r['l1'][i]) / math.log(r['N'] / q['N']):.2f}" for i in range(3))
        print(f"{r['N']:6d} {r['l1'][0]:10.3e} {r['l1'][1]:10.3e} {r['l1'][2]:10.3e} {(r['x_int'] - r['x_int_exact']) / r['dx']:26.2f}{order}")
    for side, e in zip(("water", "air"), res[-1]["plateau"]):
        print(f"N = {res[-1]['N']}, {side} star plateau (middle fifth): rho {e[0]:+.2e}, u {e[1]:+.2e}, p {e[2]:+.2e}")
    ok = args.tol <= 0 or res[-1]["l1"][0] <= args.tol
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
