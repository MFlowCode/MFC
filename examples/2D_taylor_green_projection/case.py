#!/usr/bin/env python3
"""
Decaying 2D Taylor-Green vortex in water with the all-Mach pressure projection.

u = U sin(x) cos(y) exp(-2 nu t), v = -U cos(x) sin(y) exp(-2 nu t) on [-pi, pi]^2,
periodic, so the kinetic energy decays as exp(-4 nu t); see analyze.py. The flow is at
Mach ~1e-3, and the projection steps at the advective limit, a few hundred times the
water acoustic one; --explicit runs HLLC at the acoustic limit as a control.
"""

import argparse
import json
import math
import sys

parser = argparse.ArgumentParser(description="2D Taylor-Green vortex, all-Mach pressure projection")
parser.add_argument("--N", type=int, default=64, help="cells per direction (default: %(default)s)")
parser.add_argument("--Re", type=float, default=100.0, help="Reynolds number U L / nu (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.25, help="advective CFL; acoustic under --explicit (default: %(default)s)")
parser.add_argument("--tend", type=float, default=10.0, help="final time [s] (default: %(default)s)")
parser.add_argument("--saves", type=int, default=10, help="number of output snapshots (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead, as a control")
args, _ = parser.parse_known_args()

U, rho = 1.0, 1000.0
gamma_w, p_inf_w, p0 = 4.4, 6.0e8, 1.0e5
mu = rho * U / args.Re
c_w = math.sqrt(gamma_w * (p0 + p_inf_w) / rho)

dx = 2 * math.pi / args.N
dt = args.cfl * dx / (c_w + U if args.explicit else U)
Nt = int(math.ceil(args.tend / dt))
dt = args.tend / Nt
print(f"dt = {dt:.3e} s = {dt * c_w / dx:.3g} x water acoustic limit, {Nt} steps", file=sys.stderr)

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": -math.pi,
            "x_domain%end": math.pi,
            "y_domain%beg": -math.pi,
            "y_domain%end": math.pi,
            "m": args.N - 1,
            "n": args.N - 1,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": max(Nt // args.saves, 1),
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "viscous": "T",
            "proj_method": "F" if args.explicit else "T",
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "patch_icpp(1)%geometry": 20,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 2 * math.pi,
            "patch_icpp(1)%length_y": 2 * math.pi,
            "patch_icpp(1)%vel(1)": U,
            "patch_icpp(1)%vel(2)": 1.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho,
            "patch_icpp(1)%alpha(1)": 1.0,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma_w - 1.0),
            "fluid_pp(1)%pi_inf": gamma_w * p_inf_w / (gamma_w - 1.0),
            "fluid_pp(1)%Re(1)": 1.0 / mu,
        }
    )
)
