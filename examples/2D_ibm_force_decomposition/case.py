#!/usr/bin/env python3
"""The smallest case that shows MFC's immersed-boundary force depending on the MPI decomposition.

A static cylinder in uniform flow, Re 40, on a grid coarse enough to run on a laptop. The body sits at the
origin, which is also the domain center, so a decomposition with an even number of ranks in a direction puts
a subdomain edge straight through the body. Run it at one rank and at four and compare `restart_data/ib_state_200.dat`:
the cylinder is symmetric about y = 0 and the lift must be zero, and the drag must not care how the domain
was cut up.
"""

import json
import os

Re, Ma, d, gamma, U, rho = 40.0, 0.1, 1.0, 1.4, 1.0, 1.0
P = rho * U**2 / (gamma * Ma**2)
L = float(os.environ.get("L", 5.0))
dx = float(os.environ.get("DX", 0.1))
N = int(2 * L / dx)
NSTEP = int(os.environ.get("NSTEP", 200))

case = {
    "run_time_info": "T",
    "parallel_io": "T",
    "prim_vars_wrt": "T",
    "ib_state_wrt": "T",
    "format": "silo",
    "precision": "double",
    "x_domain%beg": -L,
    "x_domain%end": L,
    "y_domain%beg": -L,
    "y_domain%end": L,
    "m": N - 1,
    "n": N - 1,
    "p": 0,
    "cyl_coord": "F",
    "dt": 0.3 * dx / (U + U / Ma),
    "t_step_start": 0,
    "t_step_stop": NSTEP,
    "t_step_save": NSTEP,
    "num_patches": 1,
    "num_fluids": 1,
    "model_eqns": "5eq",
    "alt_soundspeed": "F",
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1.0e-10,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": "arithmetic",
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": "hllc",
    "low_Mach": 2,
    "wave_speeds": "direct",
    "viscous": "T",
    "fd_order": int(os.environ.get("FDORDER", 4)),
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": 2 * L,
    "patch_icpp(1)%length_y": 2 * L,
    "patch_icpp(1)%vel(1)": U,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": P,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%eos": "ideal_gas",
    "fluid_pp(1)%Re(1)": Re / (d * U),
    "bc_x%beg": -7,
    "bc_x%grcbc_in": "T",
    "bc_x%vel_in(1)": U,
    "bc_x%vel_in(2)": 0.0,
    "bc_x%pres_in": P,
    "bc_x%alpha_rho_in(1)": rho,
    "bc_x%alpha_in(1)": 1.0,
    "bc_x%end": -8,
    "bc_y%beg": -8,
    "bc_y%end": -8,
    "ib": "T",
    "num_ibs": 1,
    "ib_neighborhood_radius": 3,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": 0.0,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%radius": d / 2,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%moving_ibm": 0,
}
print(json.dumps(case, indent=4))
