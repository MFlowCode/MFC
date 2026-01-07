#!/usr/bin/env python3
import math
import json

eps = 1e-8
Nx = 699
Ny = 299

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "F",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 7.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 3.0,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.8,
            "n_start": 0,
            "t_stop": 4.0,
            "t_save": 0.04,
            # Simulation Algorithm Parameters
            "num_patches": 3,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "time_stepper": 3,
            "elliptic_smoothing": "T",
            "elliptic_smoothing_iters": 50,
            "igr": "T",
            "igr_order": 5,
            "igr_iter_solver": 1,
            "num_igr_iters": 3,
            "num_igr_warm_start_iters": 30,
            "alf_factor": 10,
            "bc_x%beg": -3,  # 11,
            "bc_x%end": -3,  # 12
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "schlieren_wrt": "T",
            "fd_order": 4,
            "schlieren_alpha(1)": 0.5,
            "schlieren_alpha(2)": 0.5,
            "parallel_io": "T",
            # Patch 1: Left state
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 3,
            "patch_icpp(1)%length_x": 1,
            "patch_icpp(1)%length_y": 6,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": (1 - eps) * 1.0,
            "patch_icpp(1)%alpha_rho(2)": eps,
            "patch_icpp(1)%alpha(1)": 1 - eps,
            "patch_icpp(1)%alpha(2)": eps,
            # Patch 2: Top right state
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 4,
            "patch_icpp(2)%y_centroid": 2.25,
            "patch_icpp(2)%length_x": 6,
            "patch_icpp(2)%length_y": 1.5,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 0.1,
            "patch_icpp(2)%alpha_rho(1)": (1 - eps) * 0.125,
            "patch_icpp(2)%alpha_rho(2)": eps,
            "patch_icpp(2)%alpha(1)": 1 - eps,
            "patch_icpp(2)%alpha(2)": eps,
            # Patch 3: Bottom right state
            "patch_icpp(3)%geometry": 3,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%x_centroid": 4,
            "patch_icpp(3)%y_centroid": 0.75,
            "patch_icpp(3)%length_x": 6,
            "patch_icpp(3)%length_y": 1.5,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%pres": 0.1,
            "patch_icpp(3)%alpha_rho(1)": eps,
            "patch_icpp(3)%alpha_rho(2)": (1 - eps) * 1.0,
            "patch_icpp(3)%alpha(1)": eps,  # 0.95
            "patch_icpp(3)%alpha(2)": 1 - eps,  # 0.05,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0 / (1.5 - 1.0),
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(2)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
        }
    )
)
