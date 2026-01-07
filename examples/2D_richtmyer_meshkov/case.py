#!/usr/bin/env python3
import json
import math


mu = 1.0e-4
lambd = 1.0
time_end = 15.0
time_save = time_end / 20.0
eps = 1.0e-6

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 16.0 * lambd,
            "y_domain%beg": 0.0,
            "y_domain%end": lambd,
            "m": 4096,
            "n": 256,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.1,
            "n_start": 0,
            "t_stop": time_end,
            "t_save": time_save,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "recon_type": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-9,
            "null_weights": "F",
            "mapped_weno": "T",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -17,
            "bc_x%end": -17,
            "bc_y%beg": -15,
            "bc_y%end": -15,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Fluid #1 = Heavier Fluid
            # Fluid #2 = Lighter Fluid
            # Pre Shock
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%hcid": 208,
            "patch_icpp(1)%x_centroid": 0.5 * 0.7 * lambd,
            "patch_icpp(1)%y_centroid": 0.5 * lambd,
            "patch_icpp(1)%length_x": 0.7 * lambd,
            "patch_icpp(1)%length_y": lambd,
            "patch_icpp(1)%vel(1)": 1.24,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1.0 / 1.4,
            "patch_icpp(1)%alpha_rho(1)": 1.0 * eps,
            "patch_icpp(1)%alpha(1)": eps,
            "patch_icpp(1)%alpha_rho(2)": 1.0 * (1.0 - eps),
            "patch_icpp(1)%alpha(2)": (1.0 - eps),
            # Post Shock
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%hcid": 208,
            "patch_icpp(2)%x_centroid": 0.7 * lambd + 0.5 * 15.3 * lambd,
            "patch_icpp(2)%y_centroid": 0.5 * lambd,
            "patch_icpp(2)%length_x": 15.3 * lambd,
            "patch_icpp(2)%length_y": lambd,
            "patch_icpp(2)%vel(1)": 0.8787,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 1.6272 / 1.4,
            "patch_icpp(2)%alpha_rho(1)": 1.4112 * eps,
            "patch_icpp(2)%alpha(1)": eps,
            "patch_icpp(2)%alpha_rho(2)": 1.4112 * (1.0 - eps),
            "patch_icpp(2)%alpha(2)": (1.0 - eps),
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.093 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0e00,
            "fluid_pp(2)%gamma": 1.0e00 / (1.4 - 1.0e00),
            "fluid_pp(2)%pi_inf": 0.0e00,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1 / mu,
            "fluid_pp(2)%Re(1)": 1 / mu,
        }
    )
)
