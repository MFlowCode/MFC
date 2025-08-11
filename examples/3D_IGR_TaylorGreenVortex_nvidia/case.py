#!/usr/bin/env python3
import math
import json

N = 799
Nx = N
Ny = 2 * (N + 1) - 1
Nz = 2 * (N + 1) - 1

Re = 1600
L = 1
P0 = 101325
rho0 = 1
C0 = math.sqrt(1.4 * P0)
V0 = 0.1 * C0
mu = V0 * L / Re

cfl = 0.5
dx = 2 * math.pi * L / (Ny + 1)

dt = cfl * dx / (C0)

tC = L / V0
tEnd = 20 * tC

Nt = int(tEnd / dt)
Nt = 10


# Configuring case dictionary
print(
    json.dumps(
        {
            "rdma_mpi": "T",
            # Logistics
            "run_time_info": "F",
            # Computational Domain Parameters
            "x_domain%beg": -math.pi * L,
            "x_domain%end": math.pi * L,
            "y_domain%beg": -math.pi * L,
            "y_domain%end": math.pi * L,
            "z_domain%beg": -math.pi * L,
            "z_domain%end": math.pi * L,
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "cyl_coord": "F",
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": 10,  # Nt,
            "t_step_save": 10,  # int(Nt / 100),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "bc_z%beg": -1,
            "bc_z%end": -1,
            "igr": "T",
            "igr_order": 5,
            "igr_iter_solver": 1,
            "num_igr_iters": 3,
            "num_igr_warm_start_iters": 3,
            "alf_factor": 10,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "omega_wrt(1)": "T",
            "omega_wrt(2)": "T",
            "omega_wrt(3)": "T",
            "qm_wrt": "T",
            "fd_order": 4,
            "parallel_io": "T",
            # Patch 1: Background (AIR - 2)
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0,
            "patch_icpp(1)%y_centroid": 0,
            "patch_icpp(1)%z_centroid": 0,
            "patch_icpp(1)%length_x": 2 * math.pi * L,
            "patch_icpp(1)%length_y": 2 * math.pi * L,
            "patch_icpp(1)%length_z": 2 * math.pi * L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0,
            "patch_icpp(1)%pres": 0.0,
            "patch_icpp(1)%hcid": 380,
            "patch_icpp(1)%alpha_rho(1)": 1,
            "patch_icpp(1)%alpha(1)": 1,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1),
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(1)%Re(1)": 1 / mu,
            # NVIDIA UVM Options
            "nv_uvm_out_of_core": "T",
            "nv_uvm_igr_temps_on_gpu": 3,
            "nv_uvm_pref_gpu": "T",
        }
    )
)
