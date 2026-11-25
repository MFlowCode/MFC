#!/usr/bin/env python3
# This case file demonstrates the Laplace pressure jump of a water droplet in air. The laplace pressure jump
# in 2D is given by delta = sigma / r where delta is the pressure jump, sigma is the surface tension coefficient,
# and r is the radius of the droplet. The results of this simulation agree with theory to well within 1%
# relative error.

import argparse
import json
import math

parser = argparse.ArgumentParser(description="Laplace pressure jump of a 2D drop (quarter domain); check_laplace.py compares it with sigma/R.")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=50, help="Cells per direction (default 50); dt scales with 1/N.")
parser.add_argument("--tend", type=float, default=0.5, help="End time [s] (default 0.5).")
parser.add_argument("--saves", type=int, default=100, help="Number of saves (default 100).")
args = parser.parse_args()

l = 0.375

# Numerical setup
r0 = 0.15
x0 = 0
x1 = l
y0 = 0
y1 = l

Nx = args.N - 1
Ny = args.N - 1

eps = 1e-9

mydt = 5e-6 * 50 / args.N
Nt = round(args.tend / mydt)

# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": x0,
    "x_domain%end": x1,
    "y_domain%beg": y0,
    "y_domain%end": y1,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": max(1, Nt // args.saves),
    # Simulation Algorithm
    "model_eqns": "6eq",
    "alt_soundspeed": "F",
    "mixture_err": "T",
    "mpp_lim": "F",
    "time_stepper": "rk3",
    "weno_order": 5,
    "avg_state": "arithmetic",
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "weno_Re_flux": "F",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "bc_x%beg": -2,
    "bc_x%end": -3,
    "bc_y%beg": -2,
    "bc_y%end": -3,
    "num_patches": 2,
    "num_fluids": 2,
    "weno_avg": "T",
    "surface_tension": "T",
    # Database Structure Parameters
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    "sigma": 8,
    # Fluid Parameters (Water)
    "fluid_pp(1)%gamma": 1.0e00 / (2.1e00 - 1.0e00),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 2.1e00 * 1.0e06 / (2.1e00 - 1.0e00),
    # Fluid Parameters (Gas)
    "fluid_pp(2)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(2)%eos": "ideal_gas",
    # Air Patch
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0,
    "patch_icpp(1)%y_centroid": 0,
    "patch_icpp(1)%length_x": 2,
    "patch_icpp(1)%length_y": 2,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": 100000,
    "patch_icpp(1)%alpha_rho(1)": eps * 1000,
    "patch_icpp(1)%alpha_rho(2)": (1 - eps) * 1,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1 - eps,
    "patch_icpp(1)%cf_val": 0,
    # Water Patch
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.95,
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0,
    "patch_icpp(2)%y_centroid": 0,
    "patch_icpp(2)%radius": r0,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": 100000,
    "patch_icpp(2)%alpha_rho(1)": (1 - eps) * 1000,
    "patch_icpp(2)%alpha_rho(2)": eps * 1,
    "patch_icpp(2)%alpha(1)": (1 - eps),
    "patch_icpp(2)%alpha(2)": eps,
    "patch_icpp(2)%cf_val": 1,
}

print(json.dumps(data))
