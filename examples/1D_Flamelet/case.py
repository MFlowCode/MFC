#!/usr/bin/env python3
import json
import argparse
import math
import os
import cantera as ct

current_dir = os.path.dirname(os.path.abspath(__file__))
ctfile = "sandiego.yaml"
sol_L = ct.Solution(ctfile)
sol_L.TPX = 300, 8000, "O2:2,N2:2,H2O:5"
L = 0.016
Nx = 1199
dx = L / Nx
dt = 1e-8
Tend = 0.60e-3
NT = int(Tend / dt)
SAVE_COUNT = 1000
NS = 1000
case = {
    "run_time_info": "T",
    "x_domain%beg": -L / 2,
    "x_domain%end": +L / 2,
    "m": Nx,
    "n": 0,
    "p": 0,
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": 100,
    "parallel_io": "F",
    "model_eqns": 2,
    "num_fluids": 1,
    "num_patches": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "weno_avg": "F",
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 2,
    "avg_state": 1,
    "bc_x%beg": -8,
    "bc_x%end": -8,
    "viscous": "F",
    "files_dir": os.path.join(current_dir, "IC"),
    "file_extension": "000000",
    "chemistry": "T" if not args.chemistry else "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "T",
    "chem_params%transport_model": 2,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%hcid": 170,
    "patch_icpp(1)%x_centroid": 0,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%vel(1)": "0",
    "patch_icpp(1)%pres": 1.01325e5,
    "patch_icpp(1)%alpha(1)": 1,
    "patch_icpp(1)%alpha_rho(1)": "1",
    "fluid_pp(1)%gamma": 1.0e00 / (1.5e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0,
    "cantera_file": ctfile,
}
if __name__ == "__main__":
    print(json.dumps(case))
