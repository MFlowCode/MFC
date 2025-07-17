#!/usr/bin/env python3
import json

eps = 1e-9

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0e00,
            "x_domain%end": 4.0e00,
            "stretch_x": "T",
            "a_x": 7,
            "x_a": -2,
            "x_b": 2,
            "y_domain%beg": 0.0e00,
            "y_domain%end": 4.0e00,
            "stretch_y": "T",
            "a_y": 7,
            "y_a": -2,
            "y_b": 2,
            "m": 199,
            "n": 199,
            "p": 0,
            "dt": 5.0e-06,
            "t_step_start": 0,
            "t_step_stop": 2000,
            "t_step_save": 200,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 3,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -2,
            "bc_x%end": -7,
            "bc_y%beg": -2,
            "bc_y%end": -7,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: Base
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%hcid": 200,
            "patch_icpp(1)%x_centroid": 4.0,
            "patch_icpp(1)%y_centroid": 4.0,
            "patch_icpp(1)%length_x": 8.0,
            "patch_icpp(1)%length_y": 8.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 10,
            "patch_icpp(1)%alpha_rho(1)": (1 - eps) * 1000,
            "patch_icpp(1)%alpha_rho(2)": eps * 1,
            "patch_icpp(1)%alpha(1)": 1 - eps,
            "patch_icpp(1)%alpha(2)": eps,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (2.35e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 2.35e00 * 1.0e09 / (2.35e00 - 1.0e00),
            "fluid_pp(2)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(2)%pi_inf": 0.0e00,
        }
    )
)
