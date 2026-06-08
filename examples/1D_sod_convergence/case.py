#!/usr/bin/env python3
"""
1D Sod shock tube convergence case.

Standard Sod problem: rho_L=1, u_L=0, p_L=1; rho_R=0.125, u_R=0, p_R=0.1.
Discontinuity at x=0.5, gamma=1.4, T_end=0.2 (shock, contact, rarefaction).
Used for L1 self-convergence study; outflow BCs (-3) at both ends.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="1D Sod shock tube convergence case")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=128, help="Grid points (default: 128)")
parser.add_argument("--order", type=int, default=5, help="WENO order: 1, 3, 5, or 7 (default: 5)")
parser.add_argument("--muscl", action="store_true", help="Use MUSCL-2 instead of WENO")
parser.add_argument("--muscl-lim", type=int, default=1, help="MUSCL limiter: 1=minmod 2=MC 4=VanLeer 5=SUPERBEE (default: 1)")
parser.add_argument("--teno", action="store_true", help="Use TENO instead of WENO")
parser.add_argument("--teno-ct", type=float, default=1e-6, help="TENO CT threshold (default: 1e-6)")
parser.add_argument("--cfl", type=float, default=0.3, help="CFL number (default: 0.3)")
args = parser.parse_args()

gamma = 1.4
N = args.N
m = N - 1
L = 1.0
dx = L / N

c_max = math.sqrt(gamma) + 1.0  # conservative: sound speed + max velocity
dt = args.cfl * dx / c_max
T_end = 0.2
Nt = max(4, math.ceil(T_end / dt))
dt = T_end / Nt

if args.muscl:
    scheme_params = {
        "recon_type": 2,
        "muscl_order": 2,
        "muscl_lim": args.muscl_lim,
    }
else:
    scheme_params = {
        "recon_type": 1,
        "weno_order": args.order,
        "weno_eps": 1.0e-40,
        "weno_Re_flux": "F",
        "weno_avg": "F",
        "mapped_weno": "F" if (args.order == 1 or args.teno) else "T",
        "null_weights": "F",
        "mp_weno": "F",
        "teno": "T" if args.teno else "F",
        **({"teno_CT": args.teno_ct} if args.teno else {}),
    }

print(
    json.dumps(
        {
            "run_time_info": "F",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "m": m,
            "n": 0,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": Nt,
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.25,
            "patch_icpp(1)%length_x": 0.5,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": 0.1,
            "patch_icpp(2)%alpha_rho(1)": 0.125,
            "patch_icpp(2)%alpha(1)": 1.0,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            **scheme_params,
        }
    )
)
