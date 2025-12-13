#!/usr/bin/env python3
import json

# -----------------------------------------------------------------------------
# 2.1 Gaussian Divergence Pulse (1D via uniform y-slice)
# -----------------------------------------------------------------------------
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain: [0,1]×[0,1] but n=1 → effectively 1D in x
            "x_domain%beg":-0.5,
            "x_domain%end": 0.5,
            "y_domain%beg":-0.5,
            "y_domain%end": 0.5,
            "m": 100,
            "n": 100,
            "p": 0,
            "dt": 0.005,
            "t_step_start": 0,
            "t_step_stop": 50,
            "t_step_save": 1,
            # Numerical Method
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 1,
            "weno_order": 1,
            "weno_eps": 1.0e-6,
            # "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 1,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            # I/O
            # "format": 1,
            "format": 2,
            "precision": 2,
            "prim_vars_wrt": "T",
            "rho_wrt": "T",
            "parallel_io": "T",
            # Physics
            "mhd": "T",
            "hyper_cleaning": "T",
            "hyper_cleaning_speed": 1.5,
            "hyper_cleaning_tau": 0.04,
            # --- Patch 1: Gaussian Divergence Pulse IC ---
            "patch_icpp(1)%geometry": 3,   # Cartesian
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%Bx": 0.0,
            "patch_icpp(1)%By": 0.0,
            "patch_icpp(1)%Bz": 0.0,
            # Fluid EOS
            "fluid_pp(1)%gamma": 1.0e00 / (5.0 / 3.0 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        },
    )
)
