#!/usr/bin/env python3
import json
import math

# Dynamic Viscosity
Mu1 = 0.0000184
rho1 = 1.19
gam_a = 1.4
# Patch Design
D = 0.1

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -1 * D,
            "x_domain%end": 1 * D,
            "y_domain%beg": -0.3 * D,
            "y_domain%end": 0.3 * D,
            "m": 599,
            "n": 199,
            "p": 0,
            "dt": 5.0e-7,
            "t_step_start": 0,
            "t_step_stop": 6000,
            "t_step_save": 20,
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
            "viscous": "T",
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
            "patch_icpp(1)%vel(1)": 20,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": 100000,
            "patch_icpp(1)%alpha_rho(1)": (1.0) * rho1,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_ib(1)%geometry": 5,
            "patch_ib(1)%model_filepath": "Character_IBM.stl",
            "patch_ib(1)%model_translate(1)": -0.0469828735,
            "patch_ib(1)%model_translate(2)": -0.0094892600,
            "patch_ib(1)%model_spc": 200,
            "patch_ib(1)%model_threshold": 0.99,
            "patch_ib(1)%slip": "F",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(1)%Re(1)": 10000,
        }
    )
)
