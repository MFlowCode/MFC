#!/usr/bin/env python3
import json

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0e00,
            "x_domain%end": 1.0e00,
            "y_domain%beg": 0.0e00,
            "y_domain%end": 1.0e00,
            "m": 99,
            "n": 99,
            "p": 0,
            "dt": 5.0e-07,
            "t_step_start": 0,
            "t_step_stop": 1000,
            "t_step_save": 100,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 3,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "recon_type": 2,
            "muscl_order": 2,
            "muscl_lim": 2,
            "int_comp": "T",
            "null_weights": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: Base
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5e00,
            "patch_icpp(1)%y_centroid": 0.5e00,
            "patch_icpp(1)%length_x": 1.0e00,
            "patch_icpp(1)%length_y": 1.0e00,
            "patch_icpp(1)%vel(1)": 100.0e00,
            "patch_icpp(1)%vel(2)": 100.0e00,
            "patch_icpp(1)%pres": 1.0e05,
            "patch_icpp(1)%alpha_rho(1)": 1000.0e00,
            "patch_icpp(1)%alpha_rho(2)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0e-12,
            "patch_icpp(1)%alpha(2)": 1.0 - 1.0e-12,
            # Patch 2: Density to transport
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%smoothen": "T",
            "patch_icpp(2)%smooth_patch_id": 1,
            "patch_icpp(2)%smooth_coeff": 0.5e00,
            "patch_icpp(2)%x_centroid": 0.1e00,
            "patch_icpp(2)%y_centroid": 0.1e00,
            "patch_icpp(2)%radius": 0.1e00,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 100.0e00,
            "patch_icpp(2)%vel(2)": 100.0e00,
            "patch_icpp(2)%pres": 1.0e05,
            "patch_icpp(2)%alpha_rho(1)": 1.0,
            "patch_icpp(2)%alpha_rho(2)": 1.0,
            "patch_icpp(2)%alpha(1)": 0,
            "patch_icpp(2)%alpha(2)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (2.35e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 2.35e00 * 1.0e09 / (2.35e00 - 1.0e00),
            "fluid_pp(2)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(2)%pi_inf": 0.0e00,
        }
    )
)
