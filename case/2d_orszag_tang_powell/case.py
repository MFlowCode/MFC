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
            "x_domain%beg": 0,
            "x_domain%end": 1.,
            "y_domain%beg": 0,
            "y_domain%end": 1.,
            "m": 511,
            "n": 511,
            "p": 0,
            "dt": 0.0005,
            "t_step_start": 0,
            "t_step_stop": 2000,
            "t_step_save": 10,
            # "m": 99,
            # "n": 99,
            # "p": 0,
            # "dt": 0.0025,
            # "t_step_start": 0,
            # "t_step_stop": 400,
            # "t_step_save": 2,

            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 1,
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
            "powell": "T",
            "fd_order": 2,

            # Patch 1 - Analytical for v and B
            # gamma = 5/3
            #   rho = 25/(36π)
            #     p = 5/(12π)
            #     v = (-sin(2π y), sin(2π x), 0)
            #     B = (-sin(2π y)/sqrt(4π), sin(4π x)/sqrt(4π), 0)
            "patch_icpp(1)%hcid": 250,
            "patch_icpp(1)%geometry": 7,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.,
            "patch_icpp(1)%length_y": 1.,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 5.0/(12*math.pi),
            "patch_icpp(1)%Bx": 0.0,
            "patch_icpp(1)%By": 0.0,
            "patch_icpp(1)%Bz": 0.0,
            "patch_icpp(1)%alpha_rho(1)": 25.0/(36.0*math.pi),
            "patch_icpp(1)%alpha(1)": 1.0,

            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (5.0/3.0 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
