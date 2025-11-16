#!/usr/bin/env python3
import json
import math

# Case for smooth Alfven waves from section 6.3.1 of:
# The ∇·B=0 Constraint in Shock-Capturing Magnetohydrodynamics Codes
# Gábor Tóth

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "m": 1023,
            "n": 0,
            "p": 0,
            "dt": 1.0e-5,
            "t_step_start": 0,
            "t_step_stop": 100000,
            "t_step_save": 10000,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1e-7,
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 4,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "rho_wrt": "T",
            "parallel_io": "F",
            # MHD
            "mhd": "T",
            "Bx0": 1.0,
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%hcid": 150,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 0.1,
            "patch_icpp(1)%By": 0.0,
            "patch_icpp(1)%Bz": 0.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (5.0 / 3.0 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
