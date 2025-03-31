#!/usr/bin/env python3
import json
import math

# Dynamic Viscosity
Mu1 = 0.0000184
# Mu2 = 0.01
rho1 = 1.19  # 0.2199
gam_a = 1.4
# Patch Design
D = 5

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -6 * D,
            "x_domain%end": 6 * D,
            "y_domain%beg": -3 * D,
            "y_domain%end": 3 * D,
            "m": 159,
            "n": 79,
            "p": 0,
            "dt": 1.0e-9,
            "t_step_start": 0,
            "t_step_stop": 3000,
            "t_step_save": 30,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "avg_state": 2,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "viscous": "F",
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "ib": "T",
            "num_ibs": 1,
            # Formatted Database Files Structure Parameters
            # Export primitive variables in double precision with parallel
            # I/O to minimize I/O computational time during large simulations
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch: Middle
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0,
            "patch_icpp(1)%y_centroid": 0,
            "patch_icpp(1)%length_x": 1000 * D,
            "patch_icpp(1)%length_y": 1000 * D,
            "patch_icpp(1)%vel(1)": 0.001,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": 100000,
            "patch_icpp(1)%alpha_rho(1)": (1.0) * rho1,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_ib(1)%geometry": 5,
            "patch_ib(1)%model_filepath": "Circle_IBM.stl",
            "patch_ib(1)%model_translate(1)": -0.05,
            "patch_ib(1)%model_translate(2)": -0.05,
            "patch_ib(1)%model_spc": 100,
            "patch_ib(1)%model_threshold": 0.95,
            "patch_ib(1)%slip": "F",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),
            "fluid_pp(1)%pi_inf": 0,
        }
    )
)
