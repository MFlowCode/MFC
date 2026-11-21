#!/usr/bin/env python3
# 2D circular blast with two refined levels (see docs/documentation/amr.md).
#
# A disc of high-pressure gas expands into a quiescent ambient. amr_max_level = 2 nests a second
# 2:1 level inside the first, so the shock and the contact are resolved at four times the base
# resolution. All levels advance at the same dt (lock-step), which the level-2 cells set.
import json

N = 128
dx = 1.0 / N

# CFL about 0.2 on the level-2 cells (dx/4) at the post-shock wave speed
dt = 0.1 * (dx / 4)

# Initial refined block over the disc, in level-0 cell indices (0-based, inclusive)
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
            "t_step_stop": 800,
            "t_step_save": 40,
            # Simulation algorithm
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
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # Output
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: quiescent ambient
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 0.1,
            "patch_icpp(1)%alpha_rho(1)": 0.125,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: high-pressure disc
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%x_centroid": 0.5,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%radius": 0.1,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 1.0,
            "patch_icpp(2)%alpha_rho(1)": 1.0,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluid (ideal gas, gamma = 1.4)
            "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            # --- AMR ---
            "amr": "T",
            "amr_block_beg(1)": lo,
            "amr_block_end(1)": hi,
            "amr_block_beg(2)": lo,
            "amr_block_end(2)": hi,
            "amr_regrid_int": 8,
            "amr_tag_eps": 0.05,
            # a wide tag buffer leaves level-1 blocks large enough to nest the second level
            "amr_buf": 8,
            "amr_max_level": 2,
            "amr_max_blocks": 256,
            "amr_max_grid_size": 32,
        }
    )
)
