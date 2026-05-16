#!/usr/bin/env python3
"""
1D periodic advection convergence case.

Two identical fluids (same gamma, same density = 1) with a sine-wave volume
fraction.  Since both EOS are identical, alpha_1 advects passively at u = 1
with no acoustic coupling.  After exactly one period (T = L/u = 1), the
exact solution equals the IC, so L2(q_cons_vf1(T) - q_cons_vf1(0)) is
purely the scheme's accumulated spatial truncation error.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="1D advection convergence case")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=64, help="Grid points (default: 64)")
parser.add_argument("--order", type=int, default=5, help="WENO order: 1, 3, or 5")
parser.add_argument("--muscl", action="store_true", help="Use MUSCL-2 instead of WENO")
parser.add_argument("--cfl", type=float, default=0.4, help="CFL number (default: 0.4)")
parser.add_argument("--mp-weno", action="store_true", help="Enable MP-WENO limiter")
parser.add_argument("--muscl-lim", type=int, default=1, help="MUSCL limiter: 1=minmod 2=MC 3=VanAlbada 4=VanLeer 5=Superbee")
args = parser.parse_args()

gamma = 1.4
N = args.N
m = N - 1
L = 1.0
dx = L / N

# Max wave speed: acoustic speed + convective speed
# c_sound = sqrt(gamma * p / rho) = sqrt(gamma) ≈ 1.183 (for p=1, rho=1)
c_max = math.sqrt(gamma) + 1.0
dt = args.cfl * dx / c_max
T_end = 1.0  # exactly one period: u=1, L=1
Nt = max(4, math.ceil(T_end / dt))
dt = T_end / Nt  # snap to land exactly on T_end

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
        "weno_eps": 1.0e-16,
        "weno_Re_flux": "F",
        "weno_avg": "F",
        "mapped_weno": "F" if args.order == 1 else "T",
        "null_weights": "F",
        "mp_weno": "T" if args.mp_weno else "F",
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
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%vel(1)": 1.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": "0.5 + 0.2 * sin(2.0 * pi * x / lx)",
            "patch_icpp(1)%alpha_rho(2)": "0.5 - 0.2 * sin(2.0 * pi * x / lx)",
            "patch_icpp(1)%alpha(1)": "0.5 + 0.2 * sin(2.0 * pi * x / lx)",
            "patch_icpp(1)%alpha(2)": "0.5 - 0.2 * sin(2.0 * pi * x / lx)",
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(2)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
            **scheme_params,
        }
    )
)
