#!/usr/bin/env python3
"""Circularly polarized Alfven wave, one period (t = 1) in a periodic box; check_alfven.py compares it with the initial state.

Case of section 6.3.1 of Toth, "The div B = 0 constraint in shock-capturing magnetohydrodynamics codes" (2000).
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mfc", type=json.loads, default={})
parser.add_argument("-N", type=int, default=1024, help="Cells")
parser.add_argument("--cfl", type=float, default=0.0, help="Fast-wave CFL number (0: the fixed dt = 1e-5 of the example)")
parser.add_argument("--weno", type=int, default=3, help="weno_order")
args = parser.parse_args()

# One period in steps; c_f = sqrt(gamma p / rho + |B|^2 / rho) = 1.085 bounds the fast speed
nt = 100000 if args.cfl <= 0 else math.ceil(1.085 * args.N / args.cfl)

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "m": args.N - 1,
            "n": 0,
            "p": 0,
            "dt": 1.0 / nt,
            "t_step_start": 0,
            "t_step_stop": nt,
            "t_step_save": nt // 10 if args.cfl <= 0 else nt,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": "5eq",
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": "rk3",
            "weno_order": args.weno,
            "weno_eps": 1e-7,
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": "hlld",
            "wave_speeds": "direct",
            "avg_state": "arithmetic",
            "bc_x%beg": -1,
            "bc_x%end": -1,
            # Formatted Database Files Structure Parameters
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "rho_wrt": "T",
            "parallel_io": "F",
            # MHD
            "mhd": "T",
            "Bx0": 1.0,
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%hcid": 150,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 0.1,
            "patch_icpp(1)%By": 0.0,
            "patch_icpp(1)%Bz": 0.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (5.0 / 3.0 - 1.0e00),
            "fluid_pp(1)%eos": "ideal_gas",
        }
    )
)
