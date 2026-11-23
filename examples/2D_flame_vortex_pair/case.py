#!/usr/bin/env python3
"""
Premixed H2/air flame (phi = 1, 1 atm, 300 K) hit by a counter-rotating Lamb-Oseen vortex pair, 2D, no immersed boundary.

    ./mfc.sh run examples/2D_flame_vortex_pair/case.py -n 1 -t pre_process simulation post_process -- --N 128 --tend 1e-4
    python3 examples/2D_flame_vortex_pair/analyze_fv.py examples/2D_flame_vortex_pair --case-args "--N 128 --tend 1e-4"

Square [0, L]^2, L = 3.2 mm: slip wall at x = 0 (reactants at rest), characteristic subsonic outflow at x = L (bc_x%end = -8, GRCBC relaxation
to 1 atm), periodic in y. The pair (circulation +-0.05 m^2/s, core radius 0.25 mm, centres 0.8 mm apart, y-periodic image sums,
cyclostrophic pressure) starts at x = 0.9 mm and moves toward the flame, which starts at x = 1.8 mm (Cantera profile, S_L = 2.33 m/s).
The initial condition is written to IC/ at the first run (hcid 274). Fixed dt on every grid. --gamma 0 gives the flat-flame control.
Production values of the validation study: N = 64 to 512, --tend 1e-4 to 1.5e-4.
"""

import argparse
import json
import os

import numpy as np
from fv_common import DT, HERE, P0, L, apply_set, flame_on, gas_unburnt, numerics, profile

parser = argparse.ArgumentParser(prog="2D_flame_vortex_pair")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--N", type=int, default=64, help="Cells in x (and in y unless --ny).")
parser.add_argument("--ny", type=int, default=0, help="Cells in y (0: N). An N x ny strip keeps dx (flat-flame control).")
parser.add_argument("--gamma", type=float, default=0.05, help="Circulation of each vortex [m^2/s]; 0: flat flame.")
parser.add_argument("--tend", type=float, default=1.0e-4, help="End time [s].")
parser.add_argument("--frames", type=int, default=10, help="Number of saves.")
parser.add_argument("--dt", type=float, default=DT, help="Fixed time step [s].")
parser.add_argument("--set", action="append", default=[], metavar="KEY=VAL", help="Override a case parameter.")
args, _ = parser.parse_known_args()

RC, SEP, XV = 0.25e-3, 0.8e-3, 0.9e-3  # core radius, pair separation, pair centre x
N, NY = args.N, args.ny or args.N
dx = L / N
IC = os.path.join(HERE, "IC", f"N{N}_y{NY}_g{args.gamma:g}")


def vortex_velocity(x, y, sign, yc):
    """Lamb-Oseen vortex (circulation sign*gamma) at (XV, yc), summed over y-periodic images of period NY dx."""
    u, v = np.zeros_like(x), np.zeros_like(x)
    for k in range(-3, 4):
        ddx, ddy = x - XV, y - (yc + k * NY * dx)
        r2 = ddx**2 + ddy**2 + 1e-30
        ut = sign * args.gamma / (2 * np.pi * r2) * (1 - np.exp(-r2 / RC**2))
        u += -ut * ddy
        v += ut * ddx
    return u, v


def vortex_pressure(r, rho):
    """Cyclostrophic pressure deficit p(r) - p(inf) of one Lamb-Oseen vortex."""
    rr = np.linspace(0, 40 * RC, 40001)
    ut = args.gamma / (2 * np.pi * np.maximum(rr, 1e-30)) * (1 - np.exp(-(rr**2) / RC**2))
    dp = rho * ut**2 / np.maximum(rr, 1e-30)
    dp[0] = 0.0
    P = np.concatenate([[0.0], np.cumsum(0.5 * (dp[1:] + dp[:-1]) * np.diff(rr))])
    return np.interp(r, rr, P - P[-1])


def write_ic():
    """Flame profile plus the vortex pair on the N x NY cell centres, as hcid 274 files (skipped if IC/ is current)."""
    if os.path.exists(os.path.join(IC, "prim.15.00.000000.dat")):
        return
    import cantera as ct

    f = profile()
    g = gas_unburnt()
    ns = g.n_species
    xc, yc = (np.arange(N) + 0.5) * dx, (np.arange(NY) + 0.5) * dx
    X, Y = np.meshgrid(xc, yc, indexing="ij")  # x-major files: outer x, inner y
    T1, u1, Y1 = flame_on(xc, f)
    u, v, dp = np.zeros_like(X), np.zeros_like(X), np.zeros_like(X)
    rho_u = g.density
    for s, c in zip((+1, -1), (0.5 * NY * dx - 0.5 * SEP, 0.5 * NY * dx + 0.5 * SEP)):  # pair moves +x
        a, b = vortex_velocity(X, Y, -s, c)
        u, v = u + a, v + b
        dp += vortex_pressure(np.hypot(X - XV, Y - c), rho_u)
    u += u1[:, None]
    p = P0 + dp
    T = np.repeat(T1[:, None], NY, 1)
    Ys = np.repeat(Y1[:, :, None], NY, 2)
    sa = ct.SolutionArray(g, T.size)
    sa.TPY = T.ravel(), p.ravel(), Ys.reshape(ns, -1).T
    rho = sa.density.reshape(T.shape)
    os.makedirs(IC, exist_ok=True)
    xy = np.column_stack([X.ravel(), Y.ravel()])
    for k, c in enumerate([rho, u, v, p, np.ones_like(rho)] + [Ys[i] for i in range(ns)], 1):
        np.savetxt(os.path.join(IC, f"prim.{k}.00.000000.dat"), np.column_stack([xy, c.ravel()]), fmt="%.15e")


write_ic()
g = gas_unburnt()
nsave = max(1, round(args.tend / args.frames / args.dt))
case = numerics(g, args.mfc)
case.update(
    {
        "x_domain%beg": 0.0,
        "x_domain%end": L,
        "y_domain%beg": 0.0,
        "y_domain%end": NY * dx,
        "m": N - 1,
        "n": NY - 1,
        "p": 0,
        "dt": args.dt,
        "t_step_start": 0,
        "t_step_stop": nsave * args.frames,
        "t_step_save": nsave,
        "t_step_print": nsave,
        "bc_y%beg": -1,
        "bc_y%end": -1,
        "num_patches": 1,
        "patch_icpp(1)%geometry": 3,
        "patch_icpp(1)%hcid": 274,
        "patch_icpp(1)%x_centroid": 0.5 * L,
        "patch_icpp(1)%y_centroid": 0.5 * NY * dx,
        "patch_icpp(1)%length_x": L,
        "patch_icpp(1)%length_y": NY * dx,
        "patch_icpp(1)%vel(1)": 0.0,
        "patch_icpp(1)%vel(2)": 0.0,
        "patch_icpp(1)%pres": P0,
        "patch_icpp(1)%alpha_rho(1)": 1.0,
        "patch_icpp(1)%alpha(1)": 1.0,
        "files_dir": IC,
        "file_extension": "000000",
    }
)
case.update({f"chem_wrt_Y({i + 1})": "T" for i in range(g.n_species)})
apply_set(case, args.set)

if __name__ == "__main__":
    print(json.dumps(case))
