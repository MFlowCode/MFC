#!/usr/bin/env python3
"""
Conservation with AMR: a square of dense, high-pressure gas carried across a periodic box while the blocks regrid around it
(see docs/documentation/amr.md).

rho = p = 2 inside the square, 1 outside, in a uniform flow (0.5, 0.25). The square's edges lie on base-grid faces, so every
level starts from the same state. The density jump advects; the pressure jump launches shocks and rarefactions that cross
block edges and the periodic boundaries. check_conservation.py prints the drift of the total mass, momentum and energy.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="2D AMR conservation in a periodic box")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=64, help="base-grid cells per side")
parser.add_argument("--levels", type=int, default=2, help="refined levels")
parser.add_argument("--subcycle", action="store_true", help="subcycle the refined levels in time")
parser.add_argument("--steps", type=int, default=0, help="steps; 0 runs to t = 0.4")
parser.add_argument("--saves", type=int, default=4, help="saves over the run")
args = parser.parse_args()

gamma, u, v = 1.4, 0.5, 0.25
N, r = args.N, 2**args.levels
dx = 1.0 / N
# CFL 0.3 on the cells that set dt (level 0 subcycled, the finest level lock-step); 2 bounds the fastest wave
dt = 0.3 * (dx if args.subcycle else dx / r) / 2.0
Nt = args.steps or math.ceil(0.4 / dt)
Nt = args.saves * math.ceil(Nt / args.saves)


# square over base cells [N/4, 3N/8) x [5N/16, 7N/16)
lo, w = (N // 4, 5 * N // 16), N // 8

print(
    json.dumps(
        {
            "run_time_info": "F",
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "m": N - 1,
            "n": N - 1,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": Nt // args.saves,
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
            **{f"bc_{d}%{e}": -1 for d in "xy" for e in ("beg", "end")},
            "format": 2,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": u,
            "patch_icpp(1)%vel(2)": v,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "num_patches": 2,
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": (lo[0] + w / 2) / N,
            "patch_icpp(2)%y_centroid": (lo[1] + w / 2) / N,
            "patch_icpp(2)%length_x": w / N,
            "patch_icpp(2)%length_y": w / N,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": u,
            "patch_icpp(2)%vel(2)": v,
            "patch_icpp(2)%pres": 2.0,
            "patch_icpp(2)%alpha_rho(1)": 2.0,
            "patch_icpp(2)%alpha(1)": 1.0,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%eos": "ideal_gas",
            "amr": "T",
            "amr_max_level": args.levels,
            "amr_subcycle": "T" if args.subcycle else "F",
            # the initial block; regrid replaces it before the first step
            **{f"amr_block_{e}({d + 1})": lo[d] + o for d in (0, 1) for e, o in (("beg", -4), ("end", w + 3))},
            "amr_regrid_int": 2 if args.subcycle else 2 * r,
            "amr_tag_eps": 0.02,
            "amr_buf": 4,
            "amr_max_grid_size": N,
            "amr_max_blocks": 256,
        }
    )
)
