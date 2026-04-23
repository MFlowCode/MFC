#!/usr/bin/env python3
import json
import math

import cantera as ct

Lx = 0.05
Ly = 0.05

ctfile = "h2o2.yaml"
sol_L = ct.Solution(ctfile)
sol_L.TPY = 1125, ct.one_atm, "O2:0.21,N2:0.79"
# Configuring case dictionary
case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": Ly,
    "m": 699,
    "n": 699,
    "p": 0,
    "dt": 4.0e-08,
    "t_step_start": 0,
    "t_step_stop": 75000,
    "t_step_save": 4500,
    "num_patches": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "mp_weno": "F",
    "weno_order": 5,
    "weno_eps": 1e-16,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -7,
    "bc_x%end": -3,
    "bc_y%beg": -16,
    "bc_y%end": -3,
    "bc_y%isothermal_in": "T",
    "bc_y%Twall_in": 600.0,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "F",
    "cantera_file": ctfile,
    "chem_wrt_T": "T",
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%hcid": 291,
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%y_centroid": Ly / 2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0,
    "patch_icpp(1)%vel(2)": 0,
    "patch_icpp(1)%pres": 101325,
    "patch_icpp(1)%alpha_rho(1)": 1.00,
    "patch_icpp(1)%alpha(1)": 1,
    "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0.0e00,
    "viscous": "T",
    "fluid_pp(1)%Re(1)": 100000,
}
for i in range(len(sol_L.Y)):
    case[f"chem_wrt_Y({i + 1})"] = "T"
    case[f"patch_icpp(1)%Y({i + 1})"] = sol_L.Y[i]

if __name__ == "__main__":
    print(json.dumps(case))
