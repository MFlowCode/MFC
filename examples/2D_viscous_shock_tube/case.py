#!/usr/bin/env python3
import json
import math


mu = 1.0e-3
time_end = 1.0
time_save = time_end / 100.0
eps = 1.0e-6

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 0.5,
            "m": 1280,
            "n": 640,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.1,
            "n_start": 0,
            "t_stop": time_end,
            "t_save": time_save,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-10,
            "null_weights": "F",
            "mapped_weno": "T",
            "mp_weno": "F",
            "riemann_solver": 1,
            "wave_speeds": 1,
            "avg_state": 1,
            "bc_x%beg": -8,
            "bc_x%end": -15,
            "bc_y%beg": -16,
            "bc_y%end": -5,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Left
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * 0.5,
            "patch_icpp(1)%y_centroid": 0.5 * 0.5,
            "patch_icpp(1)%length_x": 0.5,
            "patch_icpp(1)%length_y": 0.5,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 120.0 / (7.0 / 5.0),
            "patch_icpp(1)%alpha_rho(1)": 120.0 * (1 - eps),
            "patch_icpp(1)%alpha(1)": 1.0 * (1 - eps),
            "patch_icpp(1)%alpha_rho(2)": 120.0 * eps,
            "patch_icpp(1)%alpha(2)": 1.0 * eps,
            # Right
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.5 + 0.5 * 0.5,
            "patch_icpp(2)%y_centroid": 0.5 * 0.5,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%length_y": 0.5,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 1.2 / (7.0 / 5.0),
            "patch_icpp(2)%alpha_rho(1)": 1.2 * eps,
            "patch_icpp(2)%alpha(1)": 1.0 * eps,
            "patch_icpp(2)%alpha_rho(2)": 1.2 * (1 - eps),
            "patch_icpp(2)%alpha(2)": 1.0 * (1 - eps),
            # Fluids Physical Parameters
            "viscous": "T",
            "fluid_pp(1)%gamma": 1.0e00 / ((7.0 / 5.0) - 1.0e00),
            "fluid_pp(1)%Re(1)": 1 / mu,
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(2)%gamma": 1.0e00 / ((7.0 / 5.0) - 1.0e00),
            "fluid_pp(2)%Re(1)": 1 / mu,
            "fluid_pp(2)%pi_inf": 0.0,
        }
    )
)
