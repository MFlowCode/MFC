#!/usr/bin/env python3
# References:
# + https://doi.org/10.1016/j.compfluid.2013.10.014: 4.3. Multi-component inert shock tube

import json
import argparse
import math

import cantera as ct

parser = argparse.ArgumentParser(prog="nD_inert_shocktube", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--no-chem", dest="chemistry", default=True, action="store_false", help="Disable chemstry.")

args = parser.parse_args()

ctfile = "h2o2.yaml"
sol_L = ct.Solution(ctfile)
sol_L.TPX = 300, 101325, "O:1"
u_r = -487.34
L = 0.12
Nx = 400
dx = L / Nx
dt = dx / abs(u_r) * 0.02
Tend = 230e-6

NT = int(Tend / dt)
SAVE_COUNT = 100
NS = NT // SAVE_COUNT

case = {
    # Logistics ================================================================
    "run_time_info": "T",
    # ==========================================================================
    # Computational Domain Parameters ==========================================
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
    # ==========================================================================
    # Simulation Algorithm Parameters ==========================================
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
    "bc_x%beg": -2,
    "bc_x%end": -3,
    "viscous": "F",
    # ==========================================================================
    # Chemistry ================================================================
    "chemistry": "T" if not args.chemistry else "T",
    "chem_params%diffusion": "F",
    "chem_params%reactions": "T",
    # ==========================================================================
    # Formatted Database Files Structure Parameters ============================
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "F",
    # ==========================================================================
    # ==========================================================================
    "patch_icpp(1)%geometry": 15,
    "patch_icpp(1)%hcid": 170,
    "patch_icpp(1)%x_centroid": L / 2,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%vel(1)": "8*exp(-400*((x-0.002/2)/0.002)**2)",
    "patch_icpp(1)%pres": 1.01315e5,
    "patch_icpp(1)%alpha(1)": 1,
    # 'patch_icpp(1)%Y(1)'           : '0.4+x/0.002*0.2',
    # 'patch_icpp(1)%Y(4)'           : '0.6-x/0.002*0.2',
    # 'patch_icpp(1)%Y(5)'           : '0.6-x/0.002*0.2',
    "patch_icpp(1)%alpha_rho(1)": "1",  #'(1.01325*10**5-1100*exp(-400*((x-0.002/2)/0.002)**2))/(8.314*1000*( (0.4+x/0.002*0.2)/2.016+(0.6-x/0.002*0.2)/32)*(300+30*exp(-400*((x-0.002/2)/0.002)**2)))',
    # ==========================================================================
    # Fluids Physical Parameters ===============================================
    "fluid_pp(1)%gamma": 1.0e00 / (1.32e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0,
    # ==========================================================================
    # Chemistry ================================================================
    "cantera_file": ctfile,
    # ==========================================================================
}


if __name__ == "__main__":
    print(json.dumps(case))
