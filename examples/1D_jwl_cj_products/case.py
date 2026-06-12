#!/usr/bin/env python3
"""1D six-equation JWL products with a moving detonation-product state.

The left state is not a static high-pressure block. It is initialized from a
simple Rankine-Hugoniot product estimate behind a right-running detonation:

    u_p = D * (1 - rho_0 / rho_products)
    p_products = p_0 + rho_0 * D * u_p

The numbers are TNT-like and intended as a solver/debugging benchmark, not a
calibrated explosive model.
"""

import json

eps = 1.0e-8

rho0_jwl = 1630.0
rho_air = 1.225
p0 = 101325.0

detonation_speed = 6900.0
compression_ratio = 1.80
rho_products = compression_ratio * rho0_jwl
u_products = detonation_speed * (1.0 - rho0_jwl / rho_products)
p_products = p0 + rho0_jwl * detonation_speed * u_products

jwl = {
    "A": 3.712e11,
    "B": 3.231e9,
    "R1": 4.15,
    "R2": 0.95,
    "omega": 0.30,
    "rho0": rho0_jwl,
    "E0": 1.0089e10,
    "air_e0": 2.5575e5,
    "air_rho0": rho_air,
    "air_gamma": 0.4,
}

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "m": 399,
            "n": 0,
            "p": 0,
            "t_step_start": 0,
            "t_step_stop": 400,
            "t_step_save": 40,
            "cfl_adap_dt": "T",
            "cfl_target": 0.35,
            "n_start": 0,
            "t_stop": 1.0e-6,
            "t_save": 2.0e-7,
            "num_patches": 2,
            "model_eqns": 3,
            "num_fluids": 2,
            "relax": "T",
            "relax_model": 6,
            "mpp_lim": "T",
            "mixture_err": "T",
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
            "cons_vars_wrt": "F",
            "rho_wrt": "T",
            "pres_wrt": "T",
            "E_wrt": "T",
            "c_wrt": "T",
            "parallel_io": "F",
            # Ambient air.
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": eps * rho0_jwl,
            "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_air,
            "patch_icpp(1)%alpha(1)": eps,
            "patch_icpp(1)%alpha(2)": 1.0 - eps,
            # Moving detonation products estimated from the jump conditions.
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.10,
            "patch_icpp(2)%length_x": 0.20,
            "patch_icpp(2)%vel(1)": u_products,
            "patch_icpp(2)%pres": p_products,
            "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_products,
            "patch_icpp(2)%alpha_rho(2)": eps * rho_air,
            "patch_icpp(2)%alpha(1)": 1.0 - eps,
            "patch_icpp(2)%alpha(2)": eps,
            # Fluid 1: fully reacted JWL products.
            "fluid_pp(1)%eos": 2,
            "fluid_pp(1)%gamma": 1.0 / 0.4,
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%cv": 613.5,
            "fluid_pp(1)%jwl_A": jwl["A"],
            "fluid_pp(1)%jwl_B": jwl["B"],
            "fluid_pp(1)%jwl_R1": jwl["R1"],
            "fluid_pp(1)%jwl_R2": jwl["R2"],
            "fluid_pp(1)%jwl_omega": jwl["omega"],
            "fluid_pp(1)%jwl_rho0": jwl["rho0"],
            "fluid_pp(1)%jwl_E0": jwl["E0"],
            "fluid_pp(1)%jwl_air_e0": jwl["air_e0"],
            "fluid_pp(1)%jwl_air_rho0": jwl["air_rho0"],
            "fluid_pp(1)%jwl_air_gamma": jwl["air_gamma"],
            # Fluid 2: ideal-gas air.
            "fluid_pp(2)%eos": 1,
            "fluid_pp(2)%gamma": 1.0 / 0.4,
            "fluid_pp(2)%pi_inf": 0.0,
        },
        indent=2,
    )
)
