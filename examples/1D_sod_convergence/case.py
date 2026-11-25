#!/usr/bin/env python3
"""
1D Riemann problems (Sod by default) for convergence studies and the exact-solution check (check_riemann.py).

Sod: rho_L=1, u_L=0, p_L=1; rho_R=0.125, u_R=0, p_R=0.1; T_end=0.2. Lax (--problem lax): (0.445, 0.698, 3.528) | (0.5, 0, 0.571),
T_end=0.14. 123 (--problem 123, Toro test 2): (1, -2, 0.4) | (1, 2, 0.4), T_end=0.15.
Discontinuity at x=0.5, gamma=1.4, outflow BCs (-3) at both ends.
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
parser.add_argument("--problem", choices=["sod", "lax", "123"], default="sod", help="Riemann problem (default: sod)")
parser.add_argument("--riemann", choices=["hll", "hllc", "lax_friedrichs"], default="hllc", help="Riemann solver (default: hllc)")
parser.add_argument("--igr", type=int, choices=[3, 5], default=None, help="Use IGR of this order instead of WENO/MUSCL (on a quasi-1D 2D strip)")
args = parser.parse_args()

# (rho, u, p) left and right, end time, bound on the signal speed (about 1.25 x the exact maximum)
problems = {
    "sod": ((1.0, 0.0, 1.0), (0.125, 0.0, 0.1), 0.2, math.sqrt(1.4) + 1.0),
    "lax": ((0.445, 0.698, 3.528), (0.5, 0.0, 0.571), 0.14, 3.3),
    "123": ((1.0, -2.0, 0.4), (1.0, 2.0, 0.4), 0.15, 3.45),
}
left, right, T_end, c_max = problems[args.problem]

gamma = 1.4
N = args.N
m = N - 1
L = 1.0
dx = L / N

dt = args.cfl * dx / c_max
Nt = max(4, math.ceil(T_end / dt))
dt = T_end / Nt

if args.igr:
    # IGR runs in 2D and 3D only: a quasi-1D strip of 8 cells, periodic in y. Without elliptic smoothing of the initial jump, Sod
    # gives NaN within 10 steps.
    scheme_params = {"igr": "T", "igr_order": args.igr, "alf_factor": 10, "num_igr_iters": 10, "num_igr_warm_start_iters": 50, "igr_iter_solver": 1}
    scheme_params.update({"elliptic_smoothing": "T", "elliptic_smoothing_iters": 10})
    scheme_params.update({"n": 7, "y_domain%beg": 0.0, "y_domain%end": 8 * dx, "bc_y%beg": -1, "bc_y%end": -1})
    for k in (1, 2):
        scheme_params.update({f"patch_icpp({k})%geometry": 3, f"patch_icpp({k})%y_centroid": 4 * dx, f"patch_icpp({k})%length_y": 8 * dx, f"patch_icpp({k})%vel(2)": 0.0})
elif args.muscl:
    scheme_params = {
        "recon_type": "muscl",
        "muscl_order": "second_order",
        "muscl_lim": args.muscl_lim,
    }
else:
    scheme_params = {
        "recon_type": "weno",
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
            "model_eqns": "5eq",
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": "rk3",
            "riemann_solver": args.riemann,
            "wave_speeds": "direct",
            "avg_state": "arithmetic",
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.25,
            "patch_icpp(1)%length_x": 0.5,
            "patch_icpp(1)%vel(1)": left[1],
            "patch_icpp(1)%pres": left[2],
            "patch_icpp(1)%alpha_rho(1)": left[0],
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%vel(1)": right[1],
            "patch_icpp(2)%pres": right[2],
            "patch_icpp(2)%alpha_rho(1)": right[0],
            "patch_icpp(2)%alpha(1)": 1.0,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%eos": "ideal_gas",
            **scheme_params,
        }
    )
)
