#!/usr/bin/env python3
import json
import math

# Dynamic Viscosity
Mu1 = 0.0000184
# Mu2 = 0.01
rho1 = 1.19  # 0.2199
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
            "x_domain%beg": -3 * D,
            "x_domain%end": 3 * D,
            "y_domain%beg": -1.5 * D,
            "y_domain%end": 1.5 * D,
            "m": 399,
            "n": 199,
            "p": 0,
            "dt": 1.0e-6,
            "t_step_start": 0,
            "t_step_stop": 1000,
            "t_step_save": 10,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": "5eq",
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "time_stepper": "rk3",
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "avg_state": "arithmetic",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": "hllc",
            "wave_speeds": "direct",
            "viscous": "T",
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "ib": "T",
            "num_ibs": 1,
            "fd_order": 2,
            # Formatted Database Files Structure Parameters
            # Export primitive variables in double precision with parallel
            # I/O to minimize I/O computational time during large simulations
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch: Middle
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0,
            "patch_icpp(1)%y_centroid": 0,
            "patch_icpp(1)%length_x": 1000 * D,
            "patch_icpp(1)%length_y": 1000 * D,
            "patch_icpp(1)%vel(1)": 527.2e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": 10918.2549,
            "patch_icpp(1)%alpha_rho(1)": (1.0) * rho1,
            "patch_icpp(1)%alpha(1)": 1.0,
            "num_stl_models": 1,
            "patch_ib(1)%geometry": 5,
            "patch_ib(1)%model_id": 1,
            "stl_models(1)%model_filepath": "Wedge2D_IBM.stl",
            "stl_models(1)%model_translate(1)": -0.0500000000,
            "stl_models(1)%model_translate(2)": -0.0373970250,
            "stl_models(1)%model_threshold": 0.95,
            "patch_ib(1)%slip": "F",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(1)%Re(1)": 7535533.2,
        }
    )
)
