#!/usr/bin/env python3
import json

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 0.001,
            "m": 199,
            "n": 0,
            "p": 0,
            "dt": 2e-9,
            "t_step_start": 0,
            "t_step_stop": 250,
            "t_step_save": 1,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "teno": "T",
            "teno_CT": 1e-8,
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -6,
            "bc_x%end": -6,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "rho_wrt": "T",
            "parallel_io": "T",
            # Patch 1 Liquid
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.0005,
            "patch_icpp(1)%length_x": 0.001,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": 1e05,
            "patch_icpp(1)%alpha_rho(1)": 1100,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Acoustic source
            "acoustic_source": "T",
            "num_source": 1,
            "acoustic(1)%support": 1,
            "acoustic(1)%loc(1)": 0.0001,
            "acoustic(1)%pulse": 2,
            "acoustic(1)%npulse": 1,
            "acoustic(1)%dir": 1.0,
            "acoustic(1)%mag": 1.0,
            "acoustic(1)%gauss_sigma_dist": 2e-8 * 1493,
            "acoustic(1)%delay": 1e-7,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 4.4e00 * 5.57e08 / (4.4e00 - 1.0e00),
        }
    )
)
