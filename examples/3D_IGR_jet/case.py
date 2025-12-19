#!/usr/bin/env python3
import math
import json

# Domain parameters

D = 2.0  # Jet diameter
Nd = 100  # Cells per jet diameter
x0 = 0.0  # x_beg coordinate
x1 = 3 * D  # x_end coordinate
y0 = -3 * D / 2  # y_beg coordinate
y1 = 3 * D / 2  # y_end coordinate
z0 = -3 * D / 2  # y_beg coordinate
z1 = 3 * D / 2  # y_end coordinate

alfFactor = 10
igrIters = 5

tS = 20  # dimensionless time

gam = 1.4  # Fluid gamma
M = 10  # Mach number0

# Reference parameters
pR = 1.0
rhoR = 1.0
c = math.sqrt(1.4 * pR / rhoR)

# Ambient parameters
pA = 1 * pR
rhoA = 1 * rhoR
velA = 0.01 * M * c

# Jet parameters
pJ = 10 * pR
velJ = M * c
rhoJ = 1 * rhoR

Nx = int(Nd * (x1 - x0) / D) - 1
Ny = int(Nd * (y1 - y0) / D) - 1
Nz = int(Nd * (z1 - z0) / D) - 1

# time_end = tS * D / (M*c)
time_end = 1.0

dx = D / Nd
dt = dx / 18

Nt = int(time_end / dt)

eps = 1e-6

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
            "t_step_stop": Nt,
            "t_step_save": int(Nt / 20),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "num_bc_patches": 0,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
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
            "precision": 2,
            "prim_vars_wrt": "T",
            "c_wrt": "T",
            "parallel_io": "T",
            # Background
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": (x1 + x0) / 2,
            "patch_icpp(1)%y_centroid": (y1 + y0) / 2,
            "patch_icpp(1)%z_centroid": (z1 + z0) / 2,
            "patch_icpp(1)%length_x": (x1 - x0),
            "patch_icpp(1)%length_y": (y1 - y0),
            "patch_icpp(1)%length_z": (z1 - z0),
            "patch_icpp(1)%hcid": 302,
            "patch_icpp(1)%vel(1)": velA,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%vel(3)": 0.0e00,
            "patch_icpp(1)%pres": pA,
            "patch_icpp(1)%alpha_rho(1)": eps,
            "patch_icpp(1)%alpha(1)": eps,
            "patch_icpp(1)%alpha_rho(2)": rhoA,
            "patch_icpp(1)%alpha(2)": 1 - eps,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(2)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(2)%pi_inf": 0.0,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 5e4,
            "fluid_pp(2)%Re(1)": 5e4,
        },
        indent=4,
    )
)
