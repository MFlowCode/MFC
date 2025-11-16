#!/usr/bin/env python3
import math
import json

# Domain parameters
D = 2.5  # Jet diameter
Nd = 711  # Cells per jet diameter

x0 = 0  # x_beg coordinate
x1 = 19 * D  # x_end coordinate
y0 = -22 * D / 2  # y_beg coordinate
y1 = 22 * D / 2  # y_end coordinate
z0 = -22 * D / 2  # y_beg coordinate
z1 = 22 * D / 2  # y_end coordinate
Nx = int(Nd * (x1 - x0) / D) - 1
Ny = int(Nd * (y1 - y0) / D) - 1
Nz = int(Nd * (z1 - z0) / D) - 1

time_end = 5
igrIters = 5

dx = D / Nd
dt = dx / 36

Nt = int(time_end / dt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": x0,
            "x_domain%end": x1,
            "y_domain%beg": y0,
            "y_domain%end": y1,
            "z_domain%beg": z0,
            "z_domain%end": z1,
            "m": int(Nx),
            "n": int(Ny),
            "p": int(Nz),
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": 10,  # Nt,
            "t_step_save": 10,  # int(Nt/50),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "num_bc_patches": 0,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "time_stepper": 3,
            "igr": "T",
            "igr_order": 3,
            "igr_pres_lim": "T",
            "igr_iter_solver": 1,
            "num_igr_iters": igrIters,
            "num_igr_warm_start_iters": 10 * igrIters,
            "alf_factor": 10,
            "bc_x%beg": -17,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 1,
            "prim_vars_wrt": "T",
            "file_per_process": "T",
            "parallel_io": "T",
            "down_sample": "T",
            # Patch
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": (x1 + x0) / 2,
            "patch_icpp(1)%y_centroid": (y1 + y0) / 2,
            "patch_icpp(1)%z_centroid": (z1 + z0) / 2,
            "patch_icpp(1)%length_x": 2 * (x1 - x0),
            "patch_icpp(1)%length_y": 2 * (y1 - y0),
            "patch_icpp(1)%length_z": 2 * (z1 - z0),
            "patch_icpp(1)%hcid": 303,
            "patch_icpp(1)%vel(1)": 0.0e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%vel(3)": 0.0e00,
            "patch_icpp(1)%pres": 1.0e00,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Perturbation of velocity field
            "simplex_perturb": "T",
            "simplex_params%perturb_vel(1)": "T",
            "simplex_params%perturb_vel_freq(1)": 3,
            "simplex_params%perturb_vel_scale(1)": 0.02,
            "simplex_params%perturb_vel_offset(1,1)": 12.3,
            "simplex_params%perturb_vel_offset(1,2)": -11.3,
            "simplex_params%perturb_vel_offset(1,3)": 34.6,
            "simplex_params%perturb_vel(2)": "T",
            "simplex_params%perturb_vel_freq(2)": 2,
            "simplex_params%perturb_vel_scale(2)": 0.02,
            "simplex_params%perturb_vel_offset(2,1)": -70.3,
            "simplex_params%perturb_vel_offset(2,2)": 33.4,
            "simplex_params%perturb_vel_offset(2,3)": -34.6,
            "simplex_params%perturb_vel(3)": "T",
            "simplex_params%perturb_vel_freq(3)": 2,
            "simplex_params%perturb_vel_scale(3)": 0.02,
            "simplex_params%perturb_vel_offset(3,1)": 123.3,
            "simplex_params%perturb_vel_offset(3,2)": -654.3,
            "simplex_params%perturb_vel_offset(3,3)": -64.5,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 5e5,
        },
        indent=4,
    )
)
