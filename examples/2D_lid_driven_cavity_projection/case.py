#!/usr/bin/env python3
"""
2D lid-driven cavity at Re = 5000 (lid speed 0.5, rho = 1, Re(1) = 1/mu = 1e4) with the all-Mach pressure projection: a single
ideal gas, no-slip walls, the top one moving. --explicit runs the explicit HLLC solver (low_Mach = 2) at the acoustic limit
instead.
"""

import argparse
import json

parser = argparse.ArgumentParser(description="2D lid-driven cavity, all-Mach pressure projection")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC (low_Mach = 2) at the acoustic limit instead of the projection")
args, _ = parser.parse_known_args()

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
            "y_domain%end": 1.0,
            "m": 199,
            "n": 199,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.8,
            "n_start": 0,
            "t_stop": 200.0,
            "t_save": 1.0,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": "5eq",
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            # "mixture_err": "T",
            "time_stepper": "rk3",
            "weno_order": 5,
            "weno_eps": 1e-16,
            "mapped_weno": "T",
            # "weno_Re_flux": "T",
            "mp_weno": "T",
            "weno_avg": "T",
            "riemann_solver": "hllc",
            "wave_speeds": "direct",
            "avg_state": "arithmetic",
            "bc_x%beg": -16,
            "bc_x%end": -16,
            "bc_y%beg": -16,
            "bc_y%end": -16,
            "bc_y%ve1": 0.5,
            "viscous": "T",
            **({"low_Mach": 2} if args.explicit else {"proj_method": "T", "proj_max_acfl": 200}),
            # Formatted Database Files Structure Parameters
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "omega_wrt(3)": "T",
            "fd_order": 4,
            "parallel_io": "T",
            # Patch 1: Base
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 5,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(1)%eos": "ideal_gas",
            "fluid_pp(1)%Re(1)": 1e4,
        }
    )
)
