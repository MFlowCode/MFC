#!/usr/bin/env python3
import json

# 1D air + JWL mixing shock tube exercising the Rocflu state-interpolated closure.
# Fluid 1 = TNT JWL products (eos=2); fluid 2 = ambient air (ideal gas, eos=1).
# A high-pressure products slug (0 <= x <= 0.3 m) drives a shock into ambient air,
# so the mixed band 0 < Y < 1 forms at the contact and rides the Rocflu closure.
jwl_A = 3.712e11
jwl_B = 3.231e9
jwl_R1 = 4.15
jwl_R2 = 0.95
jwl_omega = 0.30
jwl_rho0 = 1630.0
jwl_E0 = 1.0089e10
jwl_Cv = 613.5

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "m": 399,
            "n": 0,
            "p": 0,
            "dt": 5.0e-8,
            "t_step_start": 0,
            "t_step_stop": 600,
            "t_step_save": 150,
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "alt_soundspeed": "F",
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: ambient air background over the whole domain.
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": 101325.0,
            "patch_icpp(1)%alpha_rho(1)": 1.63e-5,
            "patch_icpp(1)%alpha_rho(2)": 1.22499998775,
            "patch_icpp(1)%alpha(1)": 1.0e-8,
            "patch_icpp(1)%alpha(2)": 0.99999999,
            # Patch 2: high-pressure JWL products slug (0 <= x <= 0.3 m).
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.15,
            "patch_icpp(2)%length_x": 0.3,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": 1.2e10,
            "patch_icpp(2)%alpha_rho(1)": 1629.9999837,
            "patch_icpp(2)%alpha_rho(2)": 1.225e-8,
            "patch_icpp(2)%alpha(1)": 0.99999999,
            "patch_icpp(2)%alpha(2)": 1.0e-8,
            # Fluid 1: TNT JWL products (Rocflu closure).
            "fluid_pp(1)%eos": 2,
            "fluid_pp(1)%gamma": 2.5,
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%cv": jwl_Cv,
            "fluid_pp(1)%jwl_A": jwl_A,
            "fluid_pp(1)%jwl_B": jwl_B,
            "fluid_pp(1)%jwl_R1": jwl_R1,
            "fluid_pp(1)%jwl_R2": jwl_R2,
            "fluid_pp(1)%jwl_omega": jwl_omega,
            "fluid_pp(1)%jwl_rho0": jwl_rho0,
            "fluid_pp(1)%jwl_E0": jwl_E0,
            "fluid_pp(1)%jwl_air_e0": 2.5575e5,
            "fluid_pp(1)%jwl_air_rho0": 1.225,
            # Fluid 2: ambient air (ideal gas).
            "fluid_pp(2)%eos": 1,
            "fluid_pp(2)%gamma": 2.5,
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%cv": 717.5,
        }
    )
)
