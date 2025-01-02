#!/usr/bin/env python3
import math
import json

# SURROUNDING FLOW
# Nondimensional parameters
Re0 = 50.0  # Reynolds number
M0 = 0.2  # Mach number

# Fluid properties
gamma = 1.4
pi_inf = 0.0

# Free stream velocity & pressure
u0 = 1.0
pres0 = 1.0 / (gamma * M0**2)

# Domain size
Lx = 59.0
Ly = 59.0
Lz = 59.0

# Number of grid cells
Nx = 191
Ny = 191
Nz = 191

# Grid spacing
dx = Lx / float(Nx)
dy = Ly / float(Ny)
dz = Lz / float(Nz)

# Time advancement
cfl = 0.1
T = 100.0
dt = cfl * dx / u0
Ntfinal = int(T / dt)
Ntstart = 0
Nfiles = 100
t_save = int(math.ceil((Ntfinal - 0) / float(Nfiles)))
Nt = t_save * Nfiles
t_step_start = Ntstart
t_step_stop = int(Nt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": Lx,
            "y_domain%beg": -Ly / 2.0,
            "y_domain%end": Ly / 2.0,
            "z_domain%beg": 0.0,
            "z_domain%end": Lz,
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "dt": dt,
            "t_step_start": t_step_start,
            "t_step_stop": t_step_stop,
            "t_step_save": t_save,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "T",
            "weno_avg": "T",
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -6,
            "bc_y%end": -6,
            "bc_z%beg": -1,
            "bc_z%end": -1,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "cons_vars_wrt": "T",
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "fd_order": 1,
            "omega_wrt(1)": "T",
            "omega_wrt(2)": "T",
            "omega_wrt(3)": "T",
            "qm_wrt": "T",
            # Patch 1
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": Lx / 2.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": Lz / 2.0,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%length_y": Ly,
            "patch_icpp(1)%length_z": Lz,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%vel(1)": u0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": pres0,
            # Mixing layer
            "mixlayer_vel_profile": "T",
            "mixlayer_vel_coef": 1.0,
            "mixlayer_domain": 1.0,
            "mixlayer_perturb": "T",
            # Fluids Physical Parameters
            # Surrounding liquid
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": Re0,
        }
    )
)
