#!/usr/bin/python
import math
import json

# Numerical setup
Nx = 399
dx = 1.0 / (1.0 * (Nx + 1))

Ny = 19
dy = 1.0 / (1.0 * (Ny + 1))

Tend = 6.4e-05
Nt = 200
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
            "y_domain%beg": -1.0e-01,
            "y_domain%end": 1.0e-01,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": mydt,
            "t_step_start": 0,
            "t_step_stop": int(Nt),
            "t_step_save": int(math.ceil(Nt / 1)),
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
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
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # Turning on Hypoelasticity
            "hypoelasticity": "T",
            "fd_order": 4,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1 L
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.25,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 0.5,
            "patch_icpp(1)%length_y": 0.2,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 100,
            "patch_icpp(1)%pres": 1e08,
            "patch_icpp(1)%alpha_rho(1)": 1000,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%tau_e(1)": 0.0,
            # Patch 2 R
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%length_y": 0.2,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": -100,
            "patch_icpp(2)%pres": 1e05,
            "patch_icpp(2)%alpha_rho(1)": 1000,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%tau_e(1)": 0.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 4.4e00 * 6.0e08 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%G": 10.0e9,
        }
    )
)
