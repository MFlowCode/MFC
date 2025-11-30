#!/usr/bin/env python3
import json
import math

# Case for smoothmagnetic vortex from:
# Implicit hybridized discontinuous Galerkin methods for
# compressible magnetohydrodynamics
# C. Ciuca, P. Fernandez, A. Christophe, N.C. Nguyen, J. Peraire

# A 2D magnetic vortex advects for a period of T=10 diagonally

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -5.0,
            "x_domain%end": 5.0,
            "y_domain%beg": -5.0,
            "y_domain%end": 5.0,
            "m": 1024,
            "n": 1024,
            "p": 0,
            "dt": 1.0e-4,
            "t_step_start": 0,
            "t_step_stop": 100000,
            "t_step_save": 1000,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-12,
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 4,
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
            "patch_icpp(1)%hcid": 252,
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 10.0,
            "patch_icpp(1)%length_y": 10.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 5.0 / (12 * math.pi),
            "patch_icpp(1)%Bx": 0.0,
            "patch_icpp(1)%By": 0.0,
            "patch_icpp(1)%Bz": 0.0,
            "patch_icpp(1)%alpha_rho(1)": 25.0 / (36.0 * math.pi),
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (5.0 / 3.0 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
