#!/usr/bin/env python3
import json

# PETN JWL products parameters. The case passes Q (J/kg) and MFC derives
# jwl_E0 = rho0*Q internally.
jwl_A = 6.17e11
jwl_B = 2.11e9
jwl_R1 = 4.40
jwl_R2 = 1.20
jwl_omega = 0.25
jwl_rho0 = 1770.0
jwl_Cv = 900.0
jwl_Q = 5.96e6

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "m": 399,
            "n": 0,
            "p": 0,
            "dt": 1.0e-8,
            "t_step_start": 0,
            "t_step_stop": 600,
            "t_step_save": 60,
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "recon_type": 2,
            "muscl_order": 2,
            "muscl_lim": 2,
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.25,
            "patch_icpp(1)%length_x": 0.5,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": 30.0e9,
            "patch_icpp(1)%alpha_rho(1)": jwl_rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": 10.0e9,
            "patch_icpp(2)%alpha_rho(1)": jwl_rho0,
            "patch_icpp(2)%alpha(1)": 1.0,
            "fluid_pp(1)%eos": 2,
            "fluid_pp(1)%gamma": 1.0 / 0.4,
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%cv": jwl_Cv,
            "fluid_pp(1)%jwl_A": jwl_A,
            "fluid_pp(1)%jwl_B": jwl_B,
            "fluid_pp(1)%jwl_R1": jwl_R1,
            "fluid_pp(1)%jwl_R2": jwl_R2,
            "fluid_pp(1)%jwl_omega": jwl_omega,
            "fluid_pp(1)%jwl_rho0": jwl_rho0,
            "fluid_pp(1)%jwl_Q": jwl_Q,
            "fluid_pp(1)%jwl_air_e0": 2.5575e5,
            "fluid_pp(1)%jwl_air_rho0": 1.225,
        }
    )
)
