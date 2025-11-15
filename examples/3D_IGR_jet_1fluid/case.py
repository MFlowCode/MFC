#!/usr/bin/env python3
import json

# Domain parameters
alfFactor = 10
igrIters = 3

NPx = 2
NPy = 2
NPz = 2

x0 = 0
x1 = 2 * NPx
y0 = -1.0 * NPy
y1 = 1.0 * NPy
z0 = -1.0 * NPz
z1 = 1.0 * NPz

N = 1383

Nx = N * NPx - 1
Ny = N * NPy - 1
Nz = N * NPz - 1

dx = (x1 - x0) / Nx
dt = dx / 20000

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "F",
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
            "t_step_save": 10,  # int(Nt / 20),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "num_bc_patches": 0,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "time_stepper": 3,
            "igr": "T",
            "igr_order": 5,
            "igr_pres_lim": "T",
            "igr_iter_solver": 1,
            "num_igr_iters": igrIters,
            "num_igr_warm_start_iters": igrIters,
            "alf_factor": 10,
            "bc_x%beg": -17,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "c_wrt": "F",
            "parallel_io": "T",
            "file_per_process": "T",
            "down_sample": "F",
            # Background
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": (x1 + x0) / 2,
            "patch_icpp(1)%y_centroid": (y1 + y0) / 2,
            "patch_icpp(1)%z_centroid": (z1 + z0) / 2,
            "patch_icpp(1)%length_x": (x1 - x0),
            "patch_icpp(1)%length_y": (y1 - y0),
            "patch_icpp(1)%length_z": (z1 - z0),
            "patch_icpp(1)%hcid": 302,
            "patch_icpp(1)%vel(1)": 1.0,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%vel(3)": 0.0e00,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0e00,
            "patch_icpp(1)%alpha(1)": 1.0e00,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 5e4,
        },
        indent=4,
    )
)
