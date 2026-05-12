#!/usr/bin/env python3
"""
2D single-fluid Euler convergence case (smooth diagonal advection).

Density wave rho = 1 + 0.2*sin(2*pi*(x+y)) with uniform velocity (u, v) = (1, 1)
on the unit square with periodic BCs.  Pressure p = 1 throughout.  For this IC
the Euler equations reduce to pure advection of all variables along the
diagonal at speed sqrt(2); the wave phase x+y - 2t returns to x+y after
T = 0.5.  L2(rho(T) - rho(0)) measures accumulated scheme spatial truncation
error.

The point of this case is to exercise WENO7 / TENO7 in 2D without the
primitive-to-conserved covariance floor that dominates the 2D isentropic
vortex (run_convergence.py) at testable resolutions.  Here all primitives are
linear in the conserved variables, so there is no floor.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="2D Euler diagonal-advection convergence case")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=64, help="Grid points per dim (default: 64)")
parser.add_argument("--order", type=int, default=5, help="WENO order: 1, 3, 5, or 7")
parser.add_argument("--muscl", action="store_true", help="Use MUSCL-2 instead of WENO")
parser.add_argument("--teno", action="store_true", help="Use TENO instead of WENO")
parser.add_argument("--teno-ct", type=float, default=1e-6, help="TENO CT threshold (default: 1e-6)")
parser.add_argument("--cfl", type=float, default=0.4, help="CFL number (default: 0.4; use 0.005 for WENO7/TENO7)")
parser.add_argument("--no-mapped", action="store_true", help="Disable mapped WENO")
parser.add_argument("--muscl-lim", type=int, default=0, help="MUSCL limiter: 0=unlimited 1=minmod ... (default: 0)")
parser.add_argument("--time-stepper", type=int, default=3, help="Time stepper: 1=Euler 2=RK2 3=RK3 (default: 3)")
parser.add_argument("--t-end", type=float, default=None, help="Override total simulation time (default: 0.5 = one diagonal period)")
args = parser.parse_args()

gamma = 1.4
N = args.N
m = N - 1
L = 1.0
dx = L / N

# Max wave speed: |u| + c with c = sqrt(gamma * p / rho_min). rho_min = 1 - 0.2.
c_max = math.sqrt(gamma / 0.8) + 1.0
dt = args.cfl * dx / c_max
T_end = args.t_end if args.t_end is not None else 0.5
Nt = max(1, math.ceil(T_end / dt))
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
        "mapped_weno": "F" if (args.order == 1 or args.no_mapped or args.teno) else "T",
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
            "y_domain%beg": 0.0,
            "y_domain%end": L,
            "m": m,
            "n": m,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": Nt,
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": args.time_stepper,
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": L,
            "patch_icpp(1)%vel(1)": 1.0,
            "patch_icpp(1)%vel(2)": 1.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": "1.0 + 0.2 * sin(2.0 * pi * (x / lx + y / ly))",
            "patch_icpp(1)%alpha(1)": 1.0,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            **scheme_params,
        }
    )
)
