#!/usr/bin/env python3
# 3D spherical blast with subcycled block-structured AMR (see docs/documentation/amr.md).
#
# A sphere of high-pressure gas expands into a quiescent ambient. The outgoing shock and the
# contact behind it are thin shells, so dynamic regridding (amr_regrid_int > 0) keeps a set of
# 2:1 blocks on them while the interior and the far field stay coarse. With amr_subcycle the
# refined level takes two steps per coarse step, so dt is set by the coarse cells.
import json

N = 64
dx = 1.0 / N

# Initial refined block over the sphere, in level-0 cell indices (0-based, inclusive)
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
            "z_domain%beg": 0.0,
            "z_domain%end": 1.0,
            "m": N - 1,
            "n": N - 1,
            "p": N - 1,
            # Coarse-cell CFL about 0.2 at the post-shock wave speed (subcycled: the fine level halves it)
            "dt": 0.1 * dx,
            "t_step_start": 0,
            "t_step_stop": 200,
            "t_step_save": 20,
            # Simulation algorithm
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 1,
            "mpp_lim": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            # Output
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: quiescent ambient
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%z_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%length_z": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            # a 5:1 pressure ratio: at 10:1 this blast goes unphysical near t = 0.3 at the fine resolution, with or without AMR
            "patch_icpp(1)%pres": 0.2,
            "patch_icpp(1)%alpha_rho(1)": 0.25,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: high-pressure sphere
            "patch_icpp(2)%geometry": 8,
            "patch_icpp(2)%x_centroid": 0.5,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%z_centroid": 0.5,
            "patch_icpp(2)%radius": 0.1,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
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
            "amr_block_beg(3)": lo,
            "amr_block_end(3)": hi,
            "amr_regrid_int": 4,
            "amr_tag_eps": 0.05,
            "amr_buf": 4,
            "amr_max_blocks": 1024,
            # two fine steps per coarse step
            "amr_subcycle": "T",
        }
    )
)
