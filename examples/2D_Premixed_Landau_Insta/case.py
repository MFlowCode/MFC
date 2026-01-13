#!/usr/bin/env python3
import json
import argparse
import math
import os
import cantera as ct

parser = argparse.ArgumentParser(
    prog="nD_inert_shocktube",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default='{}', metavar="DICT",
                    help="MFC's toolchain's internal state.")
parser.add_argument("--no-chem", dest='chemistry', default=True, action="store_false",
                    help="Disable chemstry.")

args = parser.parse_args()
current_dir = os.path.dirname(os.path.abspath(__file__))
ctfile    = 'h2o2.yaml'
sol_L     = ct.Solution(ctfile)
sol_L.TPX =  300,  101325, 'H:1'

Lx    = 0.0155147
Ly    = 0.00775735
Nx   = 1279
Ny   = 649
dx   = Lx / Nx
dy   = Ly / Ny
dt   = 0.5e-8
Tend = 4.0e-3

NT         = int(Tend / dt)
SAVE_COUNT = 1
NS         = 4


# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": 0,
    "x_domain%end":Lx,
    "y_domain%beg": 0,
    "y_domain%end": Ly,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "cyl_coord": "F",
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": int(NT/10),
    "t_step_print": 1000,
    # Simulation Algorithm
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mixture_err": "F",
    "mpp_lim": "F",
    "time_stepper": 3,
    "avg_state": 1,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "weno_Re_flux": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "bc_x%beg": -7,
    "bc_x%end": -8,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "num_patches": 1,
    "num_fluids": 1,
    "viscous": "T",
    'chemistry'                    : 'T' if not args.chemistry else 'T',
    'chem_params%diffusion'        : 'T',
    'chem_params%reactions'        : 'T',
        'files_dir': os.path.join(current_dir, 'IC'),
    'file_extention'                : '000000',
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
 "chem_wrt_T"                  : "T",
    # Fluid Parameters (Heavy Gas)
    "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0.0e00,
    "fluid_pp(1)%Re(1)": 100000,
    # Fluid Parameters (Light Gas)

    # Body Forces
 
    # Water Patch
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%hcid": 272,
    "patch_icpp(1)%x_centroid": Lx/2,
    "patch_icpp(1)%y_centroid": Ly/2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 1e5,
    "patch_icpp(1)%alpha_rho(1)": 1,
 #   "patch_icpp(1)%alpha_rho(2)": eps * 1,
     "patch_icpp(1)%alpha(1)": 1 ,
  #  "patch_icpp(1)%alpha(2)": eps,
     'cantera_file'                 : ctfile,
}

print(json.dumps(data))

