#!/usr/bin/env python3
"""
2D isentropic vortex (eps = 5). By default a vortex at rest on a stretched grid for 10 steps.

With --convect: the vortex translates at (u0, v0) through a periodic box [-5, 5]^2 on a uniform N x N grid from exact
conservative cell averages (hcid 283), for the order-of-accuracy check against the exact translated vortex (check_vortex.py).
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("--convect", action="store_true", help="Translating vortex in a periodic box (order-of-accuracy mode)")
parser.add_argument("-N", type=int, default=64, help="--convect: cells per direction (default: 64)")
parser.add_argument("--order", type=int, default=5, help="--convect: WENO order 1, 3, 5 or 7 (default: 5)")
parser.add_argument("--scheme", choices=["wenojs", "wenom", "wenoz", "teno"], default="wenom", help="--convect: WENO variant (default: wenom)")
parser.add_argument("--cfl", type=float, default=0.4, help="--convect: CFL number on |u| + c (default: 0.4)")
parser.add_argument("--t-end", type=float, default=2.0, help="--convect: end time (default: 2)")
parser.add_argument("--time-stepper", type=int, default=3, help="--convect: 1 = Euler, 2 = RK2, 3 = RK3 (default: 3)")
args = parser.parse_args()

# Parameters
epsilon = "5d0"
alpha = "1d0"
gamma = "1.4"

# Initial conditions
vel1_i = "0d0"
vel2_i = "0d0"
T_i = "1d0"
alpha_rho1_i = "1d0"
pres_i = "1d0"

# Perturbations
vel1 = f"{vel1_i} + (y - yc)*({epsilon}/(2d0*pi))*" + f"exp({alpha}*(1d0 - (x - xc)**2d0 - (y - yc)**2d0))"
vel2 = f"{vel2_i} - (x - xc)*({epsilon}/(2d0*pi))*" + f"exp({alpha}*(1d0 - (x - xc)**2d0 - (y - yc)**2d0))"
alpha_rho1 = (
    f"{alpha_rho1_i}*(1d0 - ({alpha_rho1_i}/{pres_i})*({epsilon}/(2d0*pi))*"
    + f"({epsilon}/(8d0*{alpha}*({gamma} + 1d0)*pi))*"
    + f"exp(2d0*{alpha}*(1d0 - (x - xc)**2d0"
    + "- (y - yc)**2d0))"
    + f")**{gamma}"
)
pres = (
    f"{pres_i}*(1d0 - ({alpha_rho1_i}/{pres_i})*({epsilon}/(2d0*pi))*"
    + f"({epsilon}/(8d0*{alpha}*({gamma} + 1d0)*pi))*"
    + f"exp(2d0*{alpha}*(1d0 - (x - xc)**2d0"
    + "- (y - yc)**2d0))"
    + f")**({gamma} + 1d0)"
)

# Numerical setup
Nx = 399
dx = 10.0 / (1.0 * (Nx + 1))

c = 1.4**0.5
C = 0.3
mydt = C * dx / c * 0.01
Nt = 1 / mydt

# Configuring case dictionary
case = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    "x_domain%beg": -3,
    "x_domain%end": 3,
    "y_domain%beg": -3,
    "y_domain%end": 3,
    "stretch_x": "T",
    "stretch_y": "T",
    "loops_x": 2,
    "loops_y": 2,
    "a_x": 1.03,
    "a_y": 1.03,
    "x_a": -1.5,
    "y_a": -1.5,
    "x_b": 1.5,
    "y_b": 1.5,
    "m": Nx,
    "n": Nx,
    "p": 0,
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": 10,
    "t_step_save": 1,
    # Simulation Algorithm Parameters
    "num_patches": 1,
    "model_eqns": "6eq",
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -4,
    "bc_x%end": -4,
    "bc_y%beg": -4,
    "bc_y%end": -4,
    # Formatted Database Files Structure Parameters
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "parallel_io": "F",
    "omega_wrt(3)": "T",
    "fd_order": 2,
    # Patch 1
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0,
    "patch_icpp(1)%y_centroid": 0,
    "patch_icpp(1)%length_x": 10.0,
    "patch_icpp(1)%length_y": 10.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 0.0,
    "patch_icpp(1)%alpha_rho(1)": 0.0,
    "patch_icpp(1)%hcid": 280,
    "patch_icpp(1)%alpha(1)": 1.0,
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1.0e00),
    "fluid_pp(1)%eos": "ideal_gas",
}


if args.convect:
    u0, v0, half = 1.0, 0.5, 5.0
    dx = 2.0 * half / args.N
    # |u| + c bound: the vortex adds up to 0.93 to the speed and lowers c to about 0.84
    dt = args.cfl * dx / (math.hypot(u0, v0) + 0.93 + 1.2)
    nt = math.ceil(args.t_end / dt)
    scheme = {"wenojs": {}, "wenom": {"mapped_weno": "T"}, "wenoz": {"wenoz": "T"}, "teno": {"teno": "T", "teno_CT": 1.0e-6 if args.order == 5 else 1.0e-9}}
    for k in ("stretch_x", "stretch_y", "loops_x", "loops_y", "a_x", "a_y", "x_a", "y_a", "x_b", "y_b", "omega_wrt(3)", "fd_order"):
        case.pop(k)
    case.update(
        {
            "run_time_info": "F",
            "x_domain%beg": -half,
            "x_domain%end": half,
            "y_domain%beg": -half,
            "y_domain%end": half,
            "m": args.N - 1,
            "n": args.N - 1,
            "dt": args.t_end / nt,
            "t_step_stop": nt,
            "t_step_save": nt,
            "model_eqns": "5eq",
            "time_stepper": args.time_stepper,
            "weno_order": args.order,
            "weno_eps": 1.0e-40 if args.scheme != "wenojs" else 1.0e-6,
            "mapped_weno": "F",
            **scheme[args.scheme],
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "prim_vars_wrt": "F",
            "patch_icpp(1)%length_x": 2.0 * half,
            "patch_icpp(1)%length_y": 2.0 * half,
            "patch_icpp(1)%vel(1)": u0,
            "patch_icpp(1)%vel(2)": v0,
            "patch_icpp(1)%hcid": 283,
        }
    )

print(json.dumps(case))
