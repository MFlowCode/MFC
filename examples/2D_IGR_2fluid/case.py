
#!/usr/bin/env python3
# This case file demonstrates the Laplace pressure jump of a water droplet in air. The laplace pressure jump
# in 2D is given by delta = sigma / r where delta is the pressure jump, sigma is the surface tension coefficient,
# and r is the radius of the droplet. The results of this simulation agree with theory to well within 1%
# relative error.

import math
import json

l = 1
eps = 1e-6

# Numerical setup
r0 = l / 4
x0 = -l
x1 = l
y0 = -l
y1 = l

Nx = 200
Ny = 200

dt = 1e-4

# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": x0,
    "x_domain%end": x1,
    "y_domain%beg": y0,
    "y_domain%end": y1,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "cyl_coord": "F",
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": 5000,
    "t_step_save": 50,
    # Simulation Algorithm
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mixture_err": "T",
    "mpp_lim": "F",
    "time_stepper": 3,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "num_patches": 2,
    "num_fluids": 2,
    "elliptic_smoothing": "T",
    "elliptic_smoothing_iters": 10,
    "igr": "T",
    "igr_order": 5,
    "igr_iter_solver": 1,
    "num_igr_iters": 3,
    "num_igr_warm_start_iters": 30,
    "alf_factor": 10,
    "viscous": "T",
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    # Fluid Parameters (Gas)
    "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0.0e00,
    "fluid_pp(1)%Re(1)": 1e5,
    "fluid_pp(2)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(2)%pi_inf": 0.0e00,
    "fluid_pp(2)%Re(1)": 1e5,
    # Ambient pressure
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0,
    "patch_icpp(1)%y_centroid": 0,
    "patch_icpp(1)%length_x": 2,
    "patch_icpp(1)%length_y": 2,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": 1,
    "patch_icpp(1)%alpha_rho(1)": 1 - eps,
    "patch_icpp(1)%alpha(1)": 1 - eps,
    "patch_icpp(1)%alpha_rho(2)": eps,
    "patch_icpp(1)%alpha(2)": eps,
    # High pressure
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0,
    "patch_icpp(2)%y_centroid": 0,
    "patch_icpp(2)%radius": r0,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": 9.518,
    "patch_icpp(2)%alpha_rho(1)": eps,
    "patch_icpp(2)%alpha(1)": eps,
    "patch_icpp(2)%alpha_rho(2)": (1 - eps) * 5,
    "patch_icpp(2)%alpha(2)": 1 - eps,
}

print(json.dumps(data))
