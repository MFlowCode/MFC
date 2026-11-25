#!/usr/bin/env python3
"""
Water disk in air advected diagonally through a periodic box at uniform pressure and velocity (6-equation model).

The exact solution keeps p and u uniform and translates the volume fraction; after --periods flow-throughs it equals the
initial condition. check_advection.py measures the pressure and velocity disturbance and the interface error and width.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=100, help="Cells per direction (default 100).")
parser.add_argument("--periods", type=float, default=1.0, help="Flow-through times along the diagonal (default 1).")
parser.add_argument("--cfl", type=float, default=0.4, help="Acoustic CFL number (default 0.4).")
parser.add_argument("--saves", type=int, default=20, help="Number of saves (default 20).")
parser.add_argument("--eps", type=float, default=1.0e-8, help="Volume fraction of the other fluid in each pure region (default 1e-8).")
args = parser.parse_args()

L, U, p0, eps = 1.0, 100.0, 1.0e05, args.eps
g_w, pi_w, rho_w = 2.35, 1.0e09, 1000.0  # stiffened-gas water
g_a, rho_a = 1.4, 1.0  # air
c_max = math.sqrt(g_w * (p0 + pi_w) / rho_w) + math.sqrt(2.0) * U
dx = L / args.N
Nt = math.ceil(args.periods * L / U / (args.cfl * dx / c_max))
Nt = args.saves * math.ceil(Nt / args.saves)
dt = args.periods * L / U / Nt


def state(alpha_w):
    return {"vel(1)": U, "vel(2)": U, "pres": p0, "alpha_rho(1)": alpha_w * rho_w, "alpha_rho(2)": (1.0 - alpha_w) * rho_a, "alpha(1)": alpha_w, "alpha(2)": 1.0 - alpha_w}


case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "y_domain%beg": 0.0,
    "y_domain%end": L,
    "m": args.N - 1,
    "n": args.N - 1,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": Nt // args.saves,
    "num_patches": 2,
    "model_eqns": "6eq",
    "alt_soundspeed": "F",
    "num_fluids": 2,
    "mpp_lim": "T",
    "mixture_err": "T",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "weno_Re_flux": "F",
    "weno_avg": "F",
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    # Patch 1: air
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * L,
    "patch_icpp(1)%y_centroid": 0.5 * L,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%length_y": L,
    # Patch 2: water disk, off the cell faces and the box centre
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.5,
    "patch_icpp(2)%x_centroid": 0.4627 * L,
    "patch_icpp(2)%y_centroid": 0.5381 * L,
    "patch_icpp(2)%radius": 0.2 * L,
    "patch_icpp(2)%alter_patch(1)": "T",
    "fluid_pp(1)%gamma": 1.0 / (g_w - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": g_w * pi_w / (g_w - 1.0),
    "fluid_pp(2)%gamma": 1.0 / (g_a - 1.0),
    "fluid_pp(2)%eos": "ideal_gas",
}
for i, alpha_w in ((1, eps), (2, 1.0 - eps)):
    case.update({f"patch_icpp({i})%{k}": v for k, v in state(alpha_w).items()})

print(json.dumps(case))
