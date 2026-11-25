#!/usr/bin/env python3
import argparse
import json

# An upwind differencing scheme for the equations of ideal magnetohydrodynamics
# M. Brio and C. C. Wu
# Note: HLLD is not used in the paper

parser = argparse.ArgumentParser()
parser.add_argument("--mfc", type=json.loads, default={})
parser.add_argument("-N", type=int, default=801, help="Cells (dt scales with 801 / N; end time 0.2)")
parser.add_argument("--gamma", type=float, default=2.0, help="Ratio of specific heats (2 in Brio and Wu; 5/3 in Ryu and Jones)")
args = parser.parse_args()
r = args.N / 801

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "m": args.N - 1,
            "n": 0,
            "p": 0,
            "dt": 0.0002 / r,
            "t_step_start": 0,
            "t_step_stop": round(1000 * r),
            "t_step_save": round(1000 * r) // 10,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": "5eq",
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": "rk3",
            "weno_order": 5,
            "mapped_weno": "T",
            "weno_eps": 1.0e-16,
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": "hlld",
            "wave_speeds": "direct",
            "avg_state": "arithmetic",
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # Formatted Database Files Structure Parameters
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "rho_wrt": "T",
            "parallel_io": "T",
            # MHD
            "mhd": "T",
            "Bx0": 0.75,
            # Patch 1 Left
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.25,
            "patch_icpp(1)%length_x": 0.5,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%By": 1.0,
            "patch_icpp(1)%Bz": 0.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2 Right
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": 0.1,
            "patch_icpp(2)%By": -1.0,
            "patch_icpp(2)%Bz": 0.0,
            "patch_icpp(2)%alpha_rho(1)": 0.125,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (args.gamma - 1.0e00),
            "fluid_pp(1)%eos": "ideal_gas",
        }
    )
)
