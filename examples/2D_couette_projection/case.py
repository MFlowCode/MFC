#!/usr/bin/env python3
"""
Startup of plane Couette flow at Mach 0.01 with the all-Mach pressure projection.

Fluid at rest between a fixed lower wall and an upper wall that starts moving at U (bc_y%ve1), periodic in x. Viscosity carries
the wall's momentum across the gap toward the linear profile u = U y/H (README.md compares the exact transient). rho = U = H =
1 and Re = U H/nu = 100. A stiffened gas with c = U/Mach makes the flow nearly incompressible, so the projection steps at the
flow's and viscosity's pace while --explicit (HLLC, low_Mach = 2) must resolve the sound speed, 100 times faster.
"""

import argparse
import json
import math
import sys

parser = argparse.ArgumentParser(description="Plane Couette flow startup, all-Mach pressure projection")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC (low_Mach = 2) at the acoustic limit instead of the projection")
args, _ = parser.parse_known_args()

U, H, rho, Re, mach, gamma, p0 = 1.0, 1.0, 1.0, 100.0, 0.01, 1.4, 1.0
ny, nx, tstop, saves, cfl = 64, 32, 30.0, 30, 0.5  # cells across and along the gap, end time (0.3 H^2/nu), outputs, CFL
nu = U * H / Re
c = U / mach
pi_inf = rho * c**2 / gamma - p0  # stiffened gas with sound speed c at p0

# Fixed step: the advective and viscous (4/3 nu dt sum 1/dx^2) limits, or the acoustic one with --explicit
dy = H / ny
dt = cfl * min(dy / U, dy**2 / (2 * 4 / 3 * nu)) if not args.explicit else cfl * dy / (U + c)
nt = int(math.ceil(tstop / dt))
dt = tstop / nt
print(f"dt = {dt:.3e}, {nt} steps", file=sys.stderr)

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": nx * dy,
            "y_domain%beg": 0.0,
            "y_domain%end": H,
            "m": nx - 1,
            "n": ny - 1,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": nt,
            "t_step_save": nt // saves,
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -16,
            "bc_y%end": -16,
            "bc_y%ve1": U,
            "viscous": "T",
            **({"low_Mach": 2} if args.explicit else {"proj_method": "T"}),
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * nx * dy,
            "patch_icpp(1)%y_centroid": 0.5 * H,
            "patch_icpp(1)%length_x": nx * dy,
            "patch_icpp(1)%length_y": H,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho,
            "patch_icpp(1)%alpha(1)": 1.0,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": gamma * pi_inf / (gamma - 1.0),
            "fluid_pp(1)%Re(1)": 1.0 / (rho * nu),
        }
    )
)
