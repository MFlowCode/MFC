#!/usr/bin/env python3
"""
3D Sedov-Taylor blast (see README.md).

Energy E = 1 released in a sphere of radius --r0 at the center of [-0.6, 0.6]^3 in gas at rest, rho0 = 1, p0 = 1e-3,
gamma = 1.4. --mode amr: one 2:1 level regridded onto the shock (lock-step, or --subcycle); --mode uniform at --N cells per
direction or, with --fine, 2N. MUSCL (minmod), HLLC, SSP-RK3, CFL-adaptive dt, saves every --t-save.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="3D Sedov-Taylor blast with AMR")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("--N", type=int, default=96, help="Base cells per direction (default: 96)")
parser.add_argument("--mode", choices=["uniform", "amr"], default="amr")
parser.add_argument("--fine", action="store_true", help="uniform: run at 2N cells")
parser.add_argument("--subcycle", action="store_true", help="Subcycle the refined level")
parser.add_argument("--r0", type=float, default=0.04, help="Deposit radius (default: 0.04)")
parser.add_argument("--t-end", type=float, default=0.16, help="End time (default: 0.16, R = 0.50)")
parser.add_argument("--t-save", type=float, default=0.02, help="Time between saves (default: 0.02)")
args = parser.parse_args()

gam, rho0, p0, E = 1.4, 1.0, 1.0e-3, 1.0
N = args.N
cells = 2 * N if args.fine else N
p_dep = (gam - 1.0) * E / (4.0 / 3.0 * math.pi * args.r0**3)

amr = {}
if args.mode == "amr":
    amr = {
        "amr": "T",
        **{f"amr_block_{e}({d})": v for d in (1, 2, 3) for e, v in (("beg", N // 2 - N // 16), ("end", N // 2 + N // 16 - 1))},
        # tag the absolute density jump: the shock (a jump of several rho0 per cell), not the rarefied interior
        "amr_regrid_int": 4,
        "amr_tag_abs": "T",
        "amr_tag_eps": 0.1,
        "amr_buf": 4,
        "amr_max_grid_size": N // 4,
        "amr_max_blocks": 4096,
        "amr_subcycle": "T" if args.subcycle else "F",
    }

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": -0.6,
            "x_domain%end": 0.6,
            "y_domain%beg": -0.6,
            "y_domain%end": 0.6,
            "z_domain%beg": -0.6,
            "z_domain%end": 0.6,
            "m": cells - 1,
            "n": cells - 1,
            "p": cells - 1,
            "cfl_adap_dt": "T",
            "cfl_target": 0.3,
            "n_start": 0,
            "t_stop": args.t_end,
            "t_save": args.t_save,
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 1,
            "mpp_lim": "F",
            "time_stepper": 3,
            "recon_type": 2,
            "muscl_order": 2,
            "muscl_lim": 1,
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "F",
            "parallel_io": "T",
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_x": 1.2,
            "patch_icpp(1)%length_y": 1.2,
            "patch_icpp(1)%length_z": 1.2,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(2)%geometry": 8,
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%z_centroid": 0.0,
            "patch_icpp(2)%radius": args.r0,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": p_dep,
            "patch_icpp(2)%alpha_rho(1)": rho0,
            "patch_icpp(2)%alpha(1)": 1.0,
            "fluid_pp(1)%eos": "ideal_gas",
            "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
            **amr,
        }
    )
)
