#!/usr/bin/env python3
import math
import json

gam_a = 1.4
Mu = 10000  # Define the fluid's dynamic viscosity
# Dimension of the Vortex
l = 1.0
L = math.pi * l

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            # Periodic Domain from -pi*l to pi*l
            "x_domain%beg": -L,
            "x_domain%end": L,
            "y_domain%beg": -L,
            "y_domain%end": L,
            "m": 199,
            "n": 199,
            "p": 0,
            "dt": 1.0e-08,
            "t_step_start": 0,
            "t_step_stop": 10000,
            "t_step_save": 100,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "weno_Re_flux": "T",
            "weno_avg": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: Base
            # Uncomment the configuration of the Taylor Green Vortex in
            # 2D analytical patch under m_patch.f90
            "patch_icpp(1)%geometry": 20,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 2.0 * L,
            "patch_icpp(1)%length_y": 2.0 * L,
            "patch_icpp(1)%vel(1)": 1.0e-01,  # Define the characteristic velocity of the vortex
            "patch_icpp(1)%vel(2)": l,  # Define the characteristic length of the vortex
            "patch_icpp(1)%pres": 1.0e05,
            "patch_icpp(1)%alpha_rho(1)": 1.22,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
            # Shear viscosity of STD air
            "fluid_pp(1)%Re(1)": 1 / Mu,
        }
    )
)
