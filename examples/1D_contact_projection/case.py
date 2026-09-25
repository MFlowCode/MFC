#!/usr/bin/env python3
"""
1D air/water contact (833:1) at uniform pressure and velocity, periodic.

The exact solution is the initial profile translating unchanged, so any change in
the phase densities, pressure or velocity is error. --acfl sets the fixed step as
a multiple of the water acoustic limit, which the projection exists to exceed.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="1D air/water contact, all-Mach pressure projection")
parser.add_argument("--velocity", type=float, default=5.0, help="advection velocity [m/s] (default: %(default)s)")
parser.add_argument("--acfl", type=float, default=30.0, help="dt as a multiple of the water acoustic limit (default: %(default)s)")
parser.add_argument("--steps", type=int, default=0, help="override the step count; 0 => one traversal (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="run the explicit HLLC solver instead, as a control")
args, _ = parser.parse_known_args()

gamma_w, p_inf_w, rho_w = 4.4, 6.0e8, 1000.0
gamma_a, rho_a = 1.4, 1.2
p0 = 1.0e5

N, L = 200, 1.0
dx = L / N
dt = args.acfl * dx / math.sqrt(gamma_w * (p0 + p_inf_w) / rho_w)
Nt = args.steps if args.steps > 0 else max(int(L / max(abs(args.velocity), 1.0) / dt), 1)


def patch(i, x0, lx, water):
    return {
        f"patch_icpp({i})%geometry": 1,
        f"patch_icpp({i})%x_centroid": x0,
        f"patch_icpp({i})%length_x": lx,
        f"patch_icpp({i})%vel(1)": args.velocity,
        f"patch_icpp({i})%pres": p0,
        f"patch_icpp({i})%alpha_rho(1)": rho_w if water else 0.0,
        f"patch_icpp({i})%alpha_rho(2)": 0.0 if water else rho_a,
        f"patch_icpp({i})%alpha(1)": 1.0 if water else 0.0,
        f"patch_icpp({i})%alpha(2)": 0.0 if water else 1.0,
    }


case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "m": N - 1,
    "n": 0,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": Nt,
    "num_patches": 3,
    "model_eqns": 2,
    "num_fluids": 2,
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "F",
    **patch(1, 0.5 * L, 0.5 * L, True),
    **patch(2, 0.125 * L, 0.25 * L, False),
    **patch(3, 0.875 * L, 0.25 * L, False),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%gamma": 1.0 / (gamma_w - 1.0),
    "fluid_pp(1)%pi_inf": gamma_w * p_inf_w / (gamma_w - 1.0),
    "fluid_pp(2)%eos": "ideal_gas",
    "fluid_pp(2)%gamma": 1.0 / (gamma_a - 1.0),
}

if not args.explicit:
    case["proj_method"] = "T"

print(json.dumps(case))
