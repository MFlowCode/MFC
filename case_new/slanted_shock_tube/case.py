#!/usr/bin/env python3
import json
import math

# -----------------------------------------------------------------------------
# 2.1 Gaussian Divergence Pulse (1D via uniform y-slice)
# -----------------------------------------------------------------------------
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain:
            "x_domain%beg": -1.0,
            "x_domain%end": 2.0,
            "y_domain%beg":-1.0,
            # "y_domain%end": 2.0/256.0,
            "y_domain%end": 2.0,
            "m": 256*3,
            "n": 256*3,
            "dt": 0.08/math.sqrt(5.0)/170.0/2.0,
            "t_step_start": 0,
            "t_step_stop": 170*2,
            "t_step_save": 10*2,
            # Numerical Method
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-6,
            # "wenoz": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 1,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # I/O
            # "format": 1,
            "format": 2,
            "precision": 2,
            "prim_vars_wrt": "T",
            "rho_wrt": "T",
            "parallel_io": "T",
            # Physics
            "mhd": "T",
            # "hyper_cleaning": "T",
            # "hyper_cleaning_speed": 30.0,
            # "hyper_cleaning_tau": 0.001,
            # --- Patch 1: Gaussian Divergence Pulse IC ---
            "patch_icpp(1)%hcid": 262,
            "patch_icpp(1)%geometry": 7,   # Cartesian
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 4.0,
            "patch_icpp(1)%length_y": 4.0,
            "patch_icpp(1)%vel(1)": 1.0,
            "patch_icpp(1)%vel(2)": 1.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 6.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%Bx": 0.0,
            "patch_icpp(1)%By": 0.0,
            "patch_icpp(1)%Bz": 1/math.sqrt(4*math.pi),
            # "patch_icpp(1)%Bz": 1/(4*math.pi),
            # "patch_icpp(1)%Bz": 0.0,
            # Fluid EOS
            "fluid_pp(1)%gamma": 1.0e00 / (5.0 / 3.0 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        },
    )
)
