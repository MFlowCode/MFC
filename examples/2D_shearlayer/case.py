#!/usr/bin/env python2
import json

myv = 1.0
Mu1 = 0.01
Mu2 = 0.01

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -0.500000000000000e00,
            "x_domain%end": 0.500000000000000e00,
            "y_domain%beg": -0.250000000000000e00,
            "y_domain%end": 0.250000000000000e00,
            "m": 319,
            "n": 159,
            "p": 0,
            "dt": 10.000000000000000e-7,
            "t_step_start": 0,
            "t_step_stop": int(4e5),
            "t_step_save": int(1e4),
            # Simulation Algorithm Parameters
            "model_eqns": 2,
            "num_fluids": 2,
            "num_patches": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.00000000000000e-16,
            "weno_Re_flux": "T",
            "weno_avg": "T",
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -5,
            "bc_y%end": -5,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: initialize entire domain
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.00000000000000e00,
            "patch_icpp(2)%y_centroid": 0.00000000000000e00,
            "patch_icpp(2)%length_x": 2.00000000000000e00,
            "patch_icpp(2)%length_y": 2.00000000000000e00,
            "patch_icpp(2)%vel(1)": myv,
            "patch_icpp(2)%vel(2)": 0.00000000000000e00,
            "patch_icpp(2)%pres": 1.01325000000000e05,
            "patch_icpp(2)%alpha_rho(1)": 1000.00000000000e00,
            "patch_icpp(2)%alpha_rho(2)": 1000.0 * 1e-12,
            "patch_icpp(2)%alpha(1)": 1.00000000000000e00 - 1.0e-12,
            "patch_icpp(2)%alpha(2)": 1.0e-12,
            # Patch 2: overwrite lower half plane
            "patch_icpp(1)%geometry": 4,
            "patch_icpp(1)%x_centroid": 0.00000000000000e00,
            "patch_icpp(1)%y_centroid": 0.00000000000000e00,
            #'patch_icpp(1)%length_x'       : 1.00000000000000E+00,
            #'patch_icpp(1)%length_y'       : 0.50000000000000E+00,
            "patch_icpp(1)%normal(1)": 0.00624987793326e00,
            "patch_icpp(1)%normal(2)": -0.99998046932219e00,
            #'patch_icpp(1)%smooth_patch_id': 1,
            #'patch_icpp(1)%smooth_coeff'   : 1.00000000000000E+00,
            "patch_icpp(1)%vel(1)": -myv,
            "patch_icpp(1)%vel(2)": 0.00000000000000e00,
            "patch_icpp(1)%pres": 1.01325000000000e05,
            "patch_icpp(1)%alpha_rho(1)": 1000 * 1.0e-12,
            "patch_icpp(1)%alpha_rho(2)": 1000.000000000000e00,
            "patch_icpp(1)%alpha(1)": 1.00000000000000e-12,
            "patch_icpp(1)%alpha(2)": 1 - 1.00000000000000e-12,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 4.4e00 * 6.0e08 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%Re(1)": 1 / Mu1,
            "fluid_pp(2)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(2)%pi_inf": 4.4e00 * 6.0e08 / (4.4e00 - 1.0e00),
            "fluid_pp(2)%Re(1)": 1 / Mu2,
        }
    )
)
