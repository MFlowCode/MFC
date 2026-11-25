#!/usr/bin/env python3
"""
Sod shock tube with AMR following the shock and the contact (see docs/documentation/amr.md).

rho, u, p = (1, 0, 1) | (0.125, 0, 0.1) at x = 0.5, gamma = 1.4, t = 0.2, outflow ends. --levels refined 2:1 levels regrid
every few steps onto the cells whose normalized density jump exceeds --tag-eps; --levels 0 is a uniform grid of -N cells
and --fine a uniform grid at the finest AMR spacing, N * 2**levels cells. check_sod.py compares the result with the exact
Riemann solution.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="1D Sod shock tube with AMR")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=200, help="base-grid cells")
parser.add_argument("--levels", type=int, default=2, help="refined levels; 0 = uniform grid")
parser.add_argument("--fine", action="store_true", help="uniform grid at the finest AMR spacing instead")
parser.add_argument("--subcycle", action="store_true", help="subcycle the refined levels in time")
# level 2 nests a few base cells inside level 1, so a narrower buffer leaves the shock near the level-2 edge
parser.add_argument("--buf", type=int, default=8, help="amr_buf")
parser.add_argument("--regrid", type=int, default=0, help="amr_regrid_int; 0 picks 2 (subcycled) or 2 * 2**levels")
parser.add_argument("--tag-eps", type=float, default=0.05, help="amr_tag_eps")
parser.add_argument("--cfl", type=float, default=0.4, help="CFL number on the cells that set dt")
parser.add_argument("--t-end", type=float, default=0.2, help="end time")
args = parser.parse_args()

gamma = 1.4
r = 2**args.levels
amr = args.levels > 0 and not args.fine
N = args.N * r if args.fine else args.N
dx = 1.0 / N
# subcycled, dt is the level-0 step; lock-step, the finest level's
dt = args.cfl * (dx if args.subcycle or not amr else dx / r) / (math.sqrt(gamma) + 1.0)
Nt = math.ceil(args.t_end / dt)
dt = args.t_end / Nt

case = {
    "run_time_info": "F",
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "m": N - 1,
    "n": 0,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": Nt,
    "num_patches": 2,
    "model_eqns": 2,
    "num_fluids": 1,
    "mpp_lim": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "format": 2,
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
    "fluid_pp(1)%eos": "ideal_gas",
}

if amr:
    case.update(
        {
            "amr": "T",
            "amr_max_level": args.levels,
            "amr_subcycle": "T" if args.subcycle else "F",
            # the initial block; regrid replaces it before the first step
            "amr_block_beg(1)": N // 2 - 8,
            "amr_block_end(1)": N // 2 + 7,
            # the shock moves about 0.3 base cells per base step and 0.3 finest cells per lock-step step
            "amr_regrid_int": args.regrid or (2 if args.subcycle else 2 * r),
            "amr_tag_eps": args.tag_eps,
            "amr_buf": args.buf,
            "amr_max_grid_size": N,
        }
    )

print(json.dumps(case))
