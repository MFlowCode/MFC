#!/usr/bin/env python3
import math
import json

lam = 0.2
h = 1.2
k = 2 * math.pi / lam
amp = 0.05 / k

# Numerical setup
x0 = 0.0
x1 = lam / 2
y0 = 0.0
y1 = h
z0 = 0.0
z1 = lam / 2

Nx = 99
Ny = 1199
Nz = 99

eps = 1e-6

dx = lam / 2 / (Nx + 1)
c = math.sqrt(1.4 * 1e5 / 1)
cfl = 0.4
dt = cfl * dx / c

Nt = math.ceil(0.2 / dt)
Ns = math.ceil(Nt / 100)

# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": x0,
    "x_domain%end": x1,
    "y_domain%beg": y0,
    "y_domain%end": y1,
    "z_domain%beg": z0,
    "z_domain%end": z1,
    "m": Nx,
    "n": Ny,
    "p": Nz,
    "cyl_coord": "F",
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": Ns,
    # Simulation Algorithm
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mixture_err": "T",
    "mpp_lim": "T",
    "time_stepper": 3,
    "recon_type": 2,
    "muscl_order": 2,
    "muscl_lim": 4,
    "int_comp": "T",
    "avg_state": 2,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "bc_x%beg": -2,
    "bc_x%end": -3,
    "bc_y%beg": -16,
    "bc_y%end": -16,
    "bc_z%beg": -2,
    "bc_z%end": -3,
    "num_patches": 1,
    "num_fluids": 2,
    "viscous": "T",
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    # Fluid Parameters (Heavy Gas)
    "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0.0e00,
    "fluid_pp(1)%Re(1)": 1 / 0.0219,
    # Fluid Parameters (Light Gas)
    "fluid_pp(2)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(2)%pi_inf": 0.0e00,
    "fluid_pp(2)%Re(1)": 1 / 0.0073,
    # Body Forces
    "bf_y": "T",
    "k_y": 0.0,
    "w_y": 0.0,
    "p_y": 0.0,
    "g_y": -98.1,
    # Water Patch
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%hcid": 300,
    "patch_icpp(1)%x_centroid": 0,
    "patch_icpp(1)%y_centroid": h / 2,
    "patch_icpp(1)%z_centroid": 0,
    "patch_icpp(1)%length_x": lam,
    "patch_icpp(1)%length_y": h,
    "patch_icpp(1)%length_z": h,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": 1e5,
    "patch_icpp(1)%alpha_rho(1)": (1 - eps),
    "patch_icpp(1)%alpha_rho(2)": eps * 1,
    "patch_icpp(1)%alpha(1)": 1 - eps,
    "patch_icpp(1)%alpha(2)": eps,
}

print(json.dumps(data))
