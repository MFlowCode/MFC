#!/usr/bin/env python3
"""Viscous Rayleigh-Taylor instability, heavy fluid (3 kg/m^3) over light (1 kg/m^3); check_rt.py compares the early growth with linear theory."""

import argparse
import json
import math

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=200, help="Cells across one wavelength (default 200); the domain is 6 wavelengths tall.")
parser.add_argument("--tend", type=float, default=2.0, help="End time [s] (default 2).")
parser.add_argument("--saves", type=int, default=100, help="Number of saves (default 100).")
args = parser.parse_args()

lam = 0.2
h = 1.2
k = 2 * math.pi / lam
amp = 0.05 / k

# Numerical setup
x0 = -lam / 2
x1 = lam / 2
y0 = 0.0
y1 = h

Nx = args.N - 1
Ny = 6 * args.N - 1

eps = 1e-6

dx = lam / (Nx + 1)
c = math.sqrt(1.4 * 1e5 / 1)
cfl = 0.5
dt = cfl * dx / c

Nt = math.ceil(args.tend / dt)
Ns = math.ceil(Nt / args.saves)

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
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": Ns,
    # Simulation Algorithm
    "model_eqns": "6eq",
    "alt_soundspeed": "F",
    "mixture_err": "T",
    "mpp_lim": "T",
    "time_stepper": "rk3",
    "avg_state": "arithmetic",
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -16,
    "bc_y%end": -16,
    "num_patches": 1,
    "num_fluids": 2,
    "viscous": "T",
    # Database Structure Parameters
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    # Fluid Parameters (Heavy Gas)
    "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(1)%eos": "ideal_gas",
    "fluid_pp(1)%Re(1)": 1 / 0.0219,
    # Fluid Parameters (Light Gas)
    "fluid_pp(2)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(2)%eos": "ideal_gas",
    "fluid_pp(2)%Re(1)": 1 / 0.0073,
    # Body Forces
    "bf_y": "T",
    "k_y": 0.0,
    "w_y": 0.0,
    "p_y": 0.0,
    "g_y": -9.81,
    # Water Patch
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%hcid": 204,
    "patch_icpp(1)%x_centroid": 0,
    "patch_icpp(1)%y_centroid": h / 2,
    "patch_icpp(1)%length_x": lam,
    "patch_icpp(1)%length_y": h,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 1e5,
    "patch_icpp(1)%alpha_rho(1)": (1 - eps),
    "patch_icpp(1)%alpha_rho(2)": eps * 1,
    "patch_icpp(1)%alpha(1)": 1 - eps,
    "patch_icpp(1)%alpha(2)": eps,
}

print(json.dumps(data))
