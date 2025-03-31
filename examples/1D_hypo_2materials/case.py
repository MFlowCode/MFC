#!/usr/bin/python
import math
import json

# Numerical setup
Nx = 399
dx = 1.0 / (1.0 * (Nx + 1))

Tend = 240e-06
Nt = 1000
mydt = Tend / (1.0 * Nt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0e00,
            "x_domain%end": 1.0e00,
            "m": 25,
            "n": 0,
            "p": 0,
            "dt": mydt,
            "t_step_start": 0,
            "t_step_stop": 50,
            "t_step_save": 50,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "F",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 1,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # Turning on Hypoelasticity
            "hypoelasticity": "T",
            "fd_order": 4,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1 L
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.35,
            "patch_icpp(1)%length_x": 0.7,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": 1e09,
            "patch_icpp(1)%alpha_rho(1)": 1000,
            "patch_icpp(1)%alpha_rho(2)": 0.0,
            "patch_icpp(1)%alpha(1)": 1,
            "patch_icpp(1)%alpha(2)": 0.0,
            "patch_icpp(1)%tau_e(1)": 0.0,
            # Patch 2 R
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.85,
            "patch_icpp(2)%length_x": 0.3,
            "patch_icpp(2)%vel(1)": 0,
            "patch_icpp(2)%pres": 1e05,
            "patch_icpp(2)%alpha_rho(1)": 0.0,
            "patch_icpp(2)%alpha_rho(2)": 50.0,
            "patch_icpp(2)%alpha(1)": 0,
            "patch_icpp(2)%alpha(2)": 1,
            "patch_icpp(2)%tau_e(1)": 0.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 4.4e00 * 6.0e08 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%G": 1.0e09,
            "fluid_pp(2)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%G": 0.0,
        }
    )
)
