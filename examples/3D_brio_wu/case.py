#!/usr/bin/env python3
import json

# 1D Brio-Wu Problem in z; Zero Bx
# Case is for validation and test suites

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
            "z_domain%beg": 0.0,
            "z_domain%end": 1.0,
            "m": 15,
            "n": 15,
            "p": 49,
            "dt": 0.004,
            "t_step_start": 0,
            "t_step_stop": 50,
            "t_step_save": 1,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 1,
            "weno_eps": 1.0e-16,
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 1,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "rho_wrt": "T",
            "parallel_io": "T",
            # MHD
            "mhd": "T",
            # Patch 1 Bottom
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%z_centroid": 0.25,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%length_z": 0.5,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%Bx": 0.0,
            "patch_icpp(1)%By": 1.0,
            "patch_icpp(1)%Bz": 0.75,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2 Top
            "patch_icpp(2)%geometry": 9,
            "patch_icpp(2)%x_centroid": 0.5,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%z_centroid": 0.75,
            "patch_icpp(2)%length_x": 1.0,
            "patch_icpp(2)%length_y": 1.0,
            "patch_icpp(2)%length_z": 0.5,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": 0.1,
            "patch_icpp(2)%Bx": 0.0,
            "patch_icpp(2)%By": -1.0,
            "patch_icpp(2)%Bz": 0.75,
            "patch_icpp(2)%alpha_rho(1)": 0.125,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (2.0e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
