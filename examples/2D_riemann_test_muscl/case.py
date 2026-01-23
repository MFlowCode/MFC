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
            "m": 499,
            "n": 499,
            "p": 0,
            "dt": 8e-05,
            "t_step_start": 0,
            "t_step_stop": 1000,
            "t_step_save": 100,
            # Simulation Algorithm Parameters
            "num_patches": 4,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            # "mp_weno": "F",
            "recon_type": 2,
            # "weno_order": 5,
            # "weno_eps": 1e-16,
            "muscl_order": 2,
            "muscl_lim": 1,
            "int_comp": "T",
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
            "patch_icpp(1)%x_centroid": 0.4,
            "patch_icpp(1)%y_centroid": 0.4,
            "patch_icpp(1)%length_x": 0.8,
            "patch_icpp(1)%length_y": 0.8,
            "patch_icpp(1)%vel(1)": 4 / math.sqrt(11),
            "patch_icpp(1)%vel(2)": 4 / math.sqrt(11),
            "patch_icpp(1)%pres": 9 / 310,
            "patch_icpp(1)%alpha_rho(1)": 77 / 558,
            "patch_icpp(1)%alpha(1)": 1,
            # Patch 1: Base
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.4,
            "patch_icpp(2)%y_centroid": 0.9,
            "patch_icpp(2)%length_x": 0.8,
            "patch_icpp(2)%length_y": 0.2,
            "patch_icpp(2)%vel(1)": 4 / math.sqrt(11),
            "patch_icpp(2)%vel(2)": 0,
            "patch_icpp(2)%pres": 0.3,
            "patch_icpp(2)%alpha_rho(1)": 33 / 62,
            "patch_icpp(2)%alpha(1)": 1,
            # Patch 1: Base
            "patch_icpp(3)%geometry": 3,
            "patch_icpp(3)%x_centroid": 0.9,
            "patch_icpp(3)%y_centroid": 0.4,
            "patch_icpp(3)%length_x": 0.2,
            "patch_icpp(3)%length_y": 0.8,
            "patch_icpp(3)%vel(1)": 0,
            "patch_icpp(3)%vel(2)": 4 / math.sqrt(11),
            "patch_icpp(3)%pres": 0.3,
            "patch_icpp(3)%alpha_rho(1)": 33 / 62,
            "patch_icpp(3)%alpha(1)": 1,
            # Patch 1: Base
            "patch_icpp(4)%geometry": 3,
            "patch_icpp(4)%x_centroid": 0.9,
            "patch_icpp(4)%y_centroid": 0.9,
            "patch_icpp(4)%length_x": 0.2,
            "patch_icpp(4)%length_y": 0.2,
            "patch_icpp(4)%vel(1)": 0,
            "patch_icpp(4)%vel(2)": 0,
            "patch_icpp(4)%pres": 1.5,
            "patch_icpp(4)%alpha_rho(1)": 1.5,
            "patch_icpp(4)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0e00,
        }
    )
)
