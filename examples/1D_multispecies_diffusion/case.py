#!/usr/bin/env python3
# References:
# + https://doi.org/10.1016/j.compfluid.2013.10.014: 4.4. Multicomponent diffusion test case

import json
import argparse
import math
import cantera as ct

ctfile = "gri30.yaml"
sol_L = ct.Solution(ctfile)
sol_L.TPX = 300, 8000, "O2:2,N2:2,H2O:5"

L = 0.05
Nx = 100
dx = L / Nx
dt = 0.3e-6
Tend = 0.05

NT = int(Tend / dt)
SAVE_COUNT = 2000
NS = 2000
case = {
    "run_time_info": "T",
    "x_domain%beg": 0,
    "x_domain%end": +L,
    "m": Nx,
    "n": 0,
    "p": 0,
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": NS,
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
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "viscous": "F",
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "F",
    "chem_params%transport_model": 2,  # Unity-Lewis
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "chem_wrt_T": "T",
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%hcid": 182,
    "patch_icpp(1)%x_centroid": L / 2,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%vel(1)": "0",
    "patch_icpp(1)%pres": 1.01325e5,
    "patch_icpp(1)%alpha(1)": 1,
    "patch_icpp(1)%alpha_rho(1)": 1,
    "fluid_pp(1)%gamma": 1.0e00 / (1.9326e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0,
    "cantera_file": ctfile,
}

for i in range(len(sol_L.Y)):
    case[f"chem_wrt_Y({i + 1})"] = "T"
    case[f"patch_icpp(1)%Y({i+1})"] = 0.0
if __name__ == "__main__":
    print(json.dumps(case))
