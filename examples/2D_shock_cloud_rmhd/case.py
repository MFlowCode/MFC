#!/usr/bin/env python3
import json
import math

# An HLLC Riemann solver for relativistic flows – II. Magnetohydrodynamics
# A. Mignone and G. Bodo
# Section 4.3.2 Relativistic shock–cloud interaction

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
            "m": 1000,
            "n": 600,
            "p": 0,
            "dt": 0.0002,
            "t_step_start": 0,
            "t_step_stop": 5000,
            "t_step_save": 100,
            # Simulation Algorithm Parameters
            "num_patches": 3,
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
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -2,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "rho_wrt": "T",
            "parallel_io": "T",
            # MHD
            "mhd": "T",
            # "powell": "T",
            # "fd_order": 2,
            "relativity": "T",
            # Patch 1 - Right (pre-shock)
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.8,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 0.4,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": -math.sqrt(1.0 - 1e-2),
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 1e-3,
            "patch_icpp(1)%Bx": 0.0,
            "patch_icpp(1)%By": 0.0,
            "patch_icpp(1)%Bz": 0.5,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2 - Circular density clump (within right patch; pre-shock)
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%x_centroid": 0.8,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%radius": 0.15,
            "patch_icpp(2)%vel(1)": -math.sqrt(1.0 - 1e-2),
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": 1e-3,
            "patch_icpp(2)%Bx": 0.0,
            "patch_icpp(2)%By": 0.0,
            "patch_icpp(2)%Bz": 0.5,
            "patch_icpp(2)%alpha_rho(1)": 10.0,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Patch 3 - Left (post-shock)
            "patch_icpp(3)%geometry": 3,
            "patch_icpp(3)%x_centroid": 0.3,
            "patch_icpp(3)%y_centroid": 0.5,
            "patch_icpp(3)%length_x": 0.6,
            "patch_icpp(3)%length_y": 1.0,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%vel(3)": 0.0,
            "patch_icpp(3)%pres": 127.9483,
            "patch_icpp(3)%Bx": 0.0,
            "patch_icpp(3)%By": 0.0,
            "patch_icpp(3)%Bz": -2.12971,
            "patch_icpp(3)%alpha_rho(1)": 42.5942,
            "patch_icpp(3)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (4.0 / 3.0 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
