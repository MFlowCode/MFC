#!/usr/bin/env python3
import math
import json

# Numerical setup
Nx = 1000
dx = 1.0 / (1.0 * (Nx + 1))

Tend, Nt = 1.8, 2000
mydt = Tend / (1.0 * Nt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 10.0,
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": mydt,
            "t_step_start": 0,
            "t_step_stop": int(Nt),
            "t_step_save": int(math.ceil(Nt / 10.0)),
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-40,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "rho_wrt": "T",
            "parallel_io": "T",
            # Background to cover whole domain with basic line patch
            # Patch 1 Left (0 < x < 1)
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%vel(1)": 2.629,
            "patch_icpp(1)%pres": 10.333,
            "patch_icpp(1)%alpha_rho(1)": 3.857,
            "patch_icpp(1)%alpha(1)": 1.0,
            # One analytic patch to take care of 1 < x < 10
            # Patch 2 Analytic
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 5.5,
            "patch_icpp(2)%length_x": 9.0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": 1.0,
            "patch_icpp(2)%alpha_rho(1)": 0,
            "patch_icpp(2)%hcid": 180,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
