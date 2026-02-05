#!/usr/bin/env python3
import json
import math

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
            "y_domain%end": 1.0,
            "m": 399,
            "n": 399,
            "p": 0,
            "dt": 0.0004,
            "t_step_start": 0,
            "t_step_stop": 375,
            "t_step_save": 5,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "mapped_weno": "T",
            "weno_eps": 1.0e-6,
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 1,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "rho_wrt": "T",
            "parallel_io": "T",
            # MHD
            "mhd": "T",
            "hyper_cleaning": "T",
            "hyper_cleaning_speed": 2.5,
            "hyper_cleaning_tau": 0.004,
            # Patch 1 - 2D MHD Rotor Problem
            # gamma = 1.4
            # Ambient medium (r > 0.1):
            #   rho = 1
            #   p = 1
            #   v = (0, 0, 0)
            #   B = (1, 0, 0)
            # Rotor (r <= 0.1):
            #   rho = 10
            #   v has angular velocity of 20
            "patch_icpp(1)%hcid": 252,
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%Bx": 5.0 / math.sqrt(math.pi * 4.0),
            "patch_icpp(1)%By": 0.0,
            "patch_icpp(1)%Bz": 0.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
