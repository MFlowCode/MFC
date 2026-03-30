#!/usr/bin/env python3
import json

# 3D hypoelastic Riemann problem for HLLD regression testing.
# Two-material pressure discontinuity at x=0.5:
#   Left:  liquid (fluid 1), p=1e6
#   Right: solid  (fluid 2, G=1e7), p=1e5
# Exercises HLLD star-state stress coupling from step 1.
# Very few steps to keep Debug-vs-Release error within 1e-7.

config = {
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": 0,
    "x_domain%end": 1.0,
    "y_domain%beg": 0,
    "y_domain%end": 1.0,
    "z_domain%beg": 0,
    "z_domain%end": 1.0,
    "m": 24,
    "n": 24,
    "p": 24,
    "cyl_coord": "F",
    "dt": 6.0e-6,
    "t_step_start": 0,
    "t_step_stop": 10,
    "t_step_save": 10,
    # Simulation Algorithm
    "num_patches": 2,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 2,
    "mpp_lim": "T",
    "mixture_err": "F",
    "time_stepper": 1,
    "weno_order": 1,
    "weno_eps": 1.0e-20,
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 4,
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
    "rho_wrt": "T",
    "parallel_io": "T",
    # Hypoelasticity
    "hypoelasticity": "T",
    "fd_order": 4,
    # Patch 1: Liquid fills entire domain (3D box)
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%y_centroid": 0.5,
    "patch_icpp(1)%z_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%length_y": 1.0,
    "patch_icpp(1)%length_z": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 10.0,
    "patch_icpp(1)%vel(3)": 5.0,
    "patch_icpp(1)%pres": 1e6,
    "patch_icpp(1)%tau_e(1)": 1e4,
    "patch_icpp(1)%alpha_rho(1)": 1000 * (1.0 - 1e-8),
    "patch_icpp(1)%alpha(1)": 1.0 - 1e-8,
    "patch_icpp(1)%alpha_rho(2)": 1000 * 1e-8,
    "patch_icpp(1)%alpha(2)": 1e-8,
    # Patch 2: Solid overwrites right half (3D box)
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%geometry": 9,
    "patch_icpp(2)%x_centroid": 0.75,
    "patch_icpp(2)%y_centroid": 0.5,
    "patch_icpp(2)%z_centroid": 0.5,
    "patch_icpp(2)%length_x": 0.5,
    "patch_icpp(2)%length_y": 1.0,
    "patch_icpp(2)%length_z": 1.0,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": -10.0,
    "patch_icpp(2)%vel(3)": -5.0,
    "patch_icpp(2)%pres": 1e5,
    "patch_icpp(2)%tau_e(1)": -1e4,
    "patch_icpp(2)%alpha_rho(1)": 1000 * 1e-8,
    "patch_icpp(2)%alpha(1)": 1e-8,
    "patch_icpp(2)%alpha_rho(2)": 1000 * (1.0 - 1e-8),
    "patch_icpp(2)%alpha(2)": 1.0 - 1e-8,
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 4.4e00 * 5.57e08 / (4.4e00 - 1.0e00),
    "fluid_pp(1)%G": 0.0,
    "fluid_pp(2)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
    "fluid_pp(2)%pi_inf": 4.4e00 * 5.57e08 / (4.4e00 - 1.0e00),
    "fluid_pp(2)%G": 1e7,
}

print(json.dumps(config, indent=4))
