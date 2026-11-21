#!/usr/bin/env python3
"""
A two-way circle under a constant body force: exact kinematics of a moving immersed boundary.

A circle (moving_ibm = 2, mass per unit depth 1000 kg/m) sits in air at rest and is accelerated by the body force
g_x = 1e5 m/s^2 (bf_x = T). The gas force on it is negligible against the body force, so the centroid must follow
x = g t^2 / 2 and v = g t. A wrong Runge-Kutta stage weight in the position update shows up as a position error with an
exact velocity; the check is sensitive to that and nothing else.

    ./mfc.sh run examples/2D_mibm_constant_acceleration/case.py -n 1 --gpu
    python3 examples/2D_mibm_constant_acceleration/check_kinematics.py <case directory>

Defaults are a CI-sized run (5 cells per radius, 200 fixed steps). Production-like values: R = 50 um, 10 cells per radius,
a run of 40 R/c, air at 1e5 Pa and 1.2 kg/m^3.
"""

import argparse
import json

parser = argparse.ArgumentParser(prog="2D_mibm_constant_acceleration")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--cells_per_R", type=int, default=5, help="Cells per radius.")
parser.add_argument("--steps", type=int, default=200, help="Number of time steps.")
parser.add_argument("--g", type=float, default=1.0e5, help="Body acceleration [m/s^2].")
parser.add_argument("--stride", type=int, default=10, help="ib_force_stride.")
parser.add_argument("--weno_order", type=int, default=1, help="Reconstruction order: 1 is well conditioned (default), 5 is WENO5.")
args = parser.parse_args()

gam, rho, p0 = 1.4, 1.2, 1.0e5
c = (gam * p0 / rho) ** 0.5
R = 50.0e-6
half = 6  # domain half width [R]
dx = R / args.cells_per_R
nc = 2 * half * args.cells_per_R
dt = 0.4 * dx / c  # acoustic CFL 0.4; the body speed stays far below c
t_end = args.steps * dt

case = {
    "run_time_info": "F",
    "x_domain%beg": -half * R,
    "x_domain%end": half * R,
    "y_domain%beg": -half * R,
    "y_domain%end": half * R,
    "m": nc - 1,
    "n": nc - 1,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": args.steps,
    "t_step_save": args.steps,
    "num_patches": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": args.weno_order,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "ib": "T",
    "num_ibs": 1,
    "fd_order": 2,
    "viscous": "F",
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "F",
    "parallel_io": "T",
    "ib_state_wrt": "T",
    "ib_force_wrt": "T",
    "ib_force_stride": args.stride,
    "bf_x": "T",
    "g_x": args.g,
    "k_x": 0.0,
    "w_x": 0.0,
    "p_x": 0.0,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": 2 * half * R,
    "patch_icpp(1)%length_y": 2 * half * R,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": 0.0,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%radius": R,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%mass": 1000.0,
    "patch_ib(1)%moving_ibm": 2,
    "patch_ib(1)%vel(1)": 0.0,
    "patch_ib(1)%vel(2)": 0.0,
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%eos": "ideal_gas",
}

# Above first order, WENO with the usual smoothness weights
if args.weno_order > 1:
    case.update({"weno_eps": 1.0e-16, "mapped_weno": "T", "mp_weno": "T"})

print(json.dumps(case))
