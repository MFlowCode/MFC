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
            "m": 199,
            "n": 0,
            "p": 0,
            "dt": 5.0e-08,
            "t_step_start": 0,
            "t_step_stop": 15000,
            "t_step_save": 1000,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 3,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: Left state
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.25e00,
            "patch_icpp(1)%length_x": 0.5e00,
            "patch_icpp(1)%vel(1)": -100.0e00,
            "patch_icpp(1)%pres": 1.0e05,
            "patch_icpp(1)%alpha_rho(1)": 1000.0e00 * 0.99e00,
            "patch_icpp(1)%alpha_rho(2)": 10.0e00 * 0.01e00,
            "patch_icpp(1)%alpha(1)": 0.99e00,
            "patch_icpp(1)%alpha(2)": 0.01e00,
            # Patch 2: Right state
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.75e00,
            "patch_icpp(2)%length_x": 0.5e00,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 100.0e00,
            "patch_icpp(2)%pres": 1.0e05,
            "patch_icpp(2)%alpha_rho(1)": 1000.0e00 * 0.99e00,
            "patch_icpp(2)%alpha_rho(2)": 10.0e00 * 0.01e00,
            "patch_icpp(2)%alpha(1)": 0.99e00,
            "patch_icpp(2)%alpha(2)": 0.01e00,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 4.4e00 * 6.0e08 / (4.4e00 - 1.0e00),
            "fluid_pp(2)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(2)%pi_inf": 0.0e00,
        }
    )
)
