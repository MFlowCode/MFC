#!/usr/bin/env python3
"""A case with probes, meant to be run twice in the same directory.

MFC appends to D/probe*_prim.dat when the file already exists. That is right when a run is being continued
and wrong when one is being started over: the second run's rows land on top of the first's, nothing in the
file marks the join, and the time column simply resets partway down. A reader sees one monotonic series.

Run it, run it again, and count the rows. See README.md.
"""

import json
import os

Re, Ma, gamma, U, rho = 40.0, 0.1, 1.4, 1.0, 1.0
P = rho * U**2 / (gamma * Ma**2)
L = 2.0
N = int(os.environ.get("N", 64))
NSTEP = int(os.environ.get("NSTEP", 20))
dx = 2 * L / N

case = {
    "run_time_info": "T",
    "parallel_io": "T",
    "prim_vars_wrt": "T",
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
    "fd_order": 4,
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
    "fluid_pp(1)%Re(1)": Re,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "probe_wrt": "T",
    "num_probes": 2,
    "probe(1)%x": 0.0,
    "probe(1)%y": 0.0,
    "probe(2)%x": 0.5,
    "probe(2)%y": 0.5,
}
print(json.dumps(case, indent=4))
