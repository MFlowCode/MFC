#!/usr/bin/env python3
"""Gaussian acoustic pulse from a 1D planar source in water; check_acoustic.py compares it with linear acoustics."""

import argparse
import json

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mfc", type=json.loads, default={})
parser.add_argument("-N", type=int, default=200, help="Cells (dt and steps scale with it to keep the CFL number and end time)")
parser.add_argument("--stretch", type=float, default=0.0, help="a_x of a stretched grid clustered on 0.3-0.6 mm (0: uniform)")
args = parser.parse_args()

r = args.N / 200
# The stretched grid ends past x = 1 mm (about 1.7 mm at a_x = 4), so the patch then covers 0-10 mm
stretch = {"stretch_x": "T", "a_x": args.stretch, "x_a": 0.0003, "x_b": 0.0006, "loops_x": 1, "patch_icpp(1)%x_centroid": 0.005, "patch_icpp(1)%length_x": 0.01} if args.stretch > 0 else {}

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 0.001,
            "m": args.N - 1,
            "n": 0,
            "p": 0,
            "dt": 2e-9 / r,
            "t_step_start": 0,
            "t_step_stop": round(250 * r),
            "t_step_save": round(r),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": "5eq",
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": "rk3",
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "teno": "T",
            "teno_CT": 1e-8,
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": "hllc",
            "wave_speeds": "direct",
            "avg_state": "arithmetic",
            "bc_x%beg": -6,
            "bc_x%end": -6,
            # Formatted Database Files Structure Parameters
            "format": "silo",
            "precision": "double",
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
            "acoustic(1)%gauss_sigma_time": 2e-8,
            "acoustic(1)%delay": 1e-7,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%pi_inf": 4.4e00 * 5.57e08 / (4.4e00 - 1.0e00),
            **stretch,
        }
    )
)
