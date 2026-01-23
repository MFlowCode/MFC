#!/usr/bin/env python3
import math
import json

gam = 1.4
Ma0 = 1e-3
c0 = math.sqrt(gam)
u0 = Ma0 * c0

Nx = 256
Ny = 256
Nz = 0
dx = 1.0 / Nx

cfl = 0.5
T = 1.0 / Ma0
dt = cfl * dx / (u0 + c0)
Ntfinal = int(T / dt)
Ntstart = 0
Nfiles = 10
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
            "x_domain%end": 1,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
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
            # 'mapped_weno'                   : 'T',
            "wenoz": "T",
            "riemann_solver": 2,
            "low_Mach": 1,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "cons_vars_wrt": "T",
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "fd_order": 1,
            "omega_wrt(3)": "T",
            # Patch 1
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%hcid": 203,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1,
            "patch_icpp(1)%length_y": 1,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%vel(1)": u0,
            "patch_icpp(1)%vel(2)": Ma0,
            "patch_icpp(1)%pres": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
