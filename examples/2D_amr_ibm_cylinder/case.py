#!/usr/bin/env python3
# 2D flow over a fixed cylinder with block-structured AMR (see docs/documentation/amr.md).
#
# A uniform subsonic stream meets a no-slip immersed cylinder. The blocks always contain the
# body, so its surface is resolved at the fine spacing, and dynamic regridding grows them over
# the waves the impulsive start sends out.
import json

N = 128
dx = 1.0 / N

# CFL about 0.15 on the fine cells (dx/2) at u + c ~ 1.3
dt = 0.1 * (dx / 2) / 1.3

# Initial refined block around the cylinder, in level-0 cell indices (0-based, inclusive)
lo, hi = int(0.3 * N), int(0.7 * N) - 1

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational domain
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "m": N - 1,
            "n": N - 1,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": 1000,
            "t_step_save": 50,
            # Simulation algorithm
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "fd_order": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # Output
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: uniform stream in +x
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0.1,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Immersed cylinder
            "ib": "T",
            "num_ibs": 1,
            "patch_ib(1)%geometry": 2,
            "patch_ib(1)%x_centroid": 0.5,
            "patch_ib(1)%y_centroid": 0.5,
            "patch_ib(1)%radius": 0.1,
            "patch_ib(1)%slip": "F",
            # Fluid (ideal gas, gamma = 1.4)
            "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            # --- AMR ---
            "amr": "T",
            "amr_block_beg(1)": lo,
            "amr_block_end(1)": hi,
            "amr_block_beg(2)": lo,
            "amr_block_end(2)": hi,
            "amr_regrid_int": 4,
            # the waves off the cylinder are weak: a low threshold keeps the blocks on them
            "amr_tag_eps": 1.0e-4,
            "amr_buf": 4,
            "amr_max_blocks": 64,
            "amr_max_grid_size": 64,
        }
    )
)
