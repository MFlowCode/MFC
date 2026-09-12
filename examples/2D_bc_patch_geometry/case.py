#!/usr/bin/env python3
"""A 2D boundary-condition patch, with the geometry the dimensionality actually supports.

`s_apply_boundary_patches` dispatches by dimensionality: geometry 1 (line segment) in 2D, geometry 2 (circle)
or 3 (rectangle) in 3D. A geometry belonging to the other case falls straight through the dispatch -- the
patch is never applied, and the face silently keeps whatever `bc_[xyz]` gave it.

Set GEOMETRY=3 to see the validator refuse what used to pass quietly:

    python3 case.py                 # geometry 1, valid in 2D
    GEOMETRY=3 python3 case.py      # a 3D geometry in a 2D case
"""

import json
import os

gamma, rho, U, Ma = 1.4, 1.0, 1.0, 0.6
cs = U / Ma
P = rho * cs**2 / gamma
D, L, H = 1.0, 8.0, 6.0
dx = D / 20
N, NY = int(L / dx), int(H / dx)

case = {
    "run_time_info": "T",
    "parallel_io": "T",
    "prim_vars_wrt": "T",
    "format": "silo",
    "precision": "double",
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "y_domain%beg": -0.5 * H,
    "y_domain%end": 0.5 * H,
    "m": N - 1,
    "n": NY - 1,
    "p": 0,
    "cyl_coord": "F",
    "dt": 0.2 * dx / (U + cs),
    "t_step_start": 0,
    "t_step_stop": 100,
    "t_step_save": 100,
    "num_patches": 2,
    "num_fluids": 1,
    "model_eqns": "5eq",
    "alt_soundspeed": "F",
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1.0e-10,
    "weno_avg": "T",
    "avg_state": "arithmetic",
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": "hllc",
    "low_Mach": 2,
    "wave_speeds": "direct",
    "viscous": "F",
    "fd_order": 4,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * L,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%length_y": H,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": P,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0,
    # The Dirichlet buffer is filled by pre_process from the initial condition *at the boundary face*, so the
    # inflow state has to be present there. A face initialised at rest stores rest and delivers nothing.
    "patch_icpp(2)%geometry": 3,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": 2 * dx,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%length_x": 4 * dx,
    "patch_icpp(2)%length_y": D,
    "patch_icpp(2)%vel(1)": U,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": P,
    "patch_icpp(2)%alpha_rho(1)": rho,
    "patch_icpp(2)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%eos": "ideal_gas",
    # a wall with a nozzle cut into it
    "bc_x%beg": -16,
    "num_bc_patches": 1,
    "patch_bc(1)%dir": 1,
    "patch_bc(1)%loc": -1,
    "patch_bc(1)%type": -17,
    "patch_bc(1)%geometry": int(os.environ.get("GEOMETRY", 1)),
    "patch_bc(1)%centroid(2)": 0.0,
    "patch_bc(1)%length(2)": D,
    "bc_x%end": -8,
    "bc_y%beg": -8,
    "bc_y%end": -8,
}
print(json.dumps(case, indent=4))
