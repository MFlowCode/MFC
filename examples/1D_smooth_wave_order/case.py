#!/usr/bin/env python3
"""
Smooth 1D Euler problems for the order of accuracy of the reconstruction, started from exact cell averages.

--problem advect: a density wave rho = 1 + 0.2 sin(2 pi x) carried at u = 1 through uniform p = 1 (gamma 1.4) for one
period; all primitives are linear in the conserved variables, so the reconstruction sees exact averages.
--problem gamma3: the isentropic gamma = 3 wave rho0 = 1 + 0.2 sin(2 pi x), u0 = 0, p = rho^3/3, for which u + c and u - c
each obey Burgers' equation (exact solution by characteristics, smooth until t = 1/(0.4 pi) = 0.80); run to --t-end 0.2.
--ic avg sets the conserved variables to their exact cell averages (10-point Gauss quadrature, written to IC/ next to this
file and read by hcid 170), --ic point to the point values at the cell centres. check_smooth.py measures the density error
against the exact cell averages. Each run needs its own copy of this directory (the IC files depend on -N).
"""

import argparse
import json
import math
import os

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=64)
parser.add_argument("--problem", choices=["advect", "gamma3"], default="gamma3")
parser.add_argument("--ic", choices=["avg", "point"], default="avg")
parser.add_argument("--order", type=int, default=5)
parser.add_argument("--teno", action="store_true")
parser.add_argument("--cfl", type=float, default=0.02, help="CFL number on max |u| + c")
parser.add_argument("--t-end", type=float, default=None)
args = parser.parse_args()

a, N = 0.2, args.N
dx = 1.0 / N
gamma, u0, T, smax = (1.4, 1.0, 1.0, 1.0 + math.sqrt(1.4 / 0.8)) if args.problem == "advect" else (3.0, 0.0, 0.2, 1.0 + 3 * a)
xc = (np.arange(N) + 0.5) * dx
if args.ic == "avg":  # 10-point Gauss averages of rho and E; at t = 0 either u is uniform or zero, so avg(rho u)/avg(rho) = u0
    xg, wg = np.polynomial.legendre.leggauss(10)
    X = xc[:, None] + 0.5 * dx * xg[None, :]
    rho = (wg * (1 + a * np.sin(2 * np.pi * X))).sum(-1) / 2
    pres = np.ones(N) if args.problem == "advect" else (wg * (1 + a * np.sin(2 * np.pi * X)) ** 3).sum(-1) / 2 / 3  # p from avg(E)
else:
    rho = 1 + a * np.sin(2 * np.pi * xc)
    pres = np.ones(N) if args.problem == "advect" else rho**3 / 3
ic = os.path.join(os.path.dirname(os.path.abspath(__file__)), "IC")
os.makedirs(ic, exist_ok=True)
for k, v in enumerate((rho, np.full(N, u0), pres, np.ones(N)), start=1):  # rho, u, p, alpha
    np.savetxt(os.path.join(ic, f"prim.{k}.00.000000.dat"), np.column_stack((xc, v)), fmt="%.17e")
T = args.t_end if args.t_end is not None else T
Nt = math.ceil(T / (args.cfl * dx / smax))
dt = T / Nt

print(
    json.dumps(
        {
            "run_time_info": "F",
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "m": args.N - 1,
            "n": 0,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": Nt,
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "recon_type": 1,
            "weno_order": args.order,
            "weno_eps": 1.0e-40,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "F" if (args.order == 1 or args.teno) else "T",
            "null_weights": "F",
            "mp_weno": "F",
            "teno": "T" if args.teno else "F",
            **({"teno_CT": 1e-6 if args.order == 5 else 1e-9} if args.teno else {}),
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "F",
            "cons_vars_wrt": "T",
            "parallel_io": "F",
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%hcid": 170,
            "files_dir": ic,
            "file_extension": "000000",
            "patch_icpp(1)%vel(1)": u0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
