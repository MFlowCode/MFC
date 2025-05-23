#!/usr/bin/env python3
import math
import json
import json, argparse
import cantera as ct

parser = argparse.ArgumentParser(
    prog="2D_detonation",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default='{}', metavar="DICT",
                    help="MFC's toolchain's internal state.")
parser.add_argument("--no-chem", dest='chemistry', default=True, action="store_false",
                                 help="Disable chemistry.")
parser.add_argument("--scale",   type=float,  default=1,    help="Scale.")

args = parser.parse_args()
ctfile    = 'h2o2.yaml'

sol_L     = ct.Solution(ctfile)
sol_L.DPX = 0.18075, 35594, 'H2:2,O2:1,AR:7'
Nx = 800
Ny = 200
Nz = 26

L = 0.24
dx =   L/Nx
dy = (L/8)/Ny
dz = (L/8)/Nz
dt = min(dx,dy)/abs(800)*0.05*0.1*0.2
Tend=830e-6

NT=int(Tend/dt)
SAVE_COUNT=100
NS=NT//SAVE_COUNT





# Configuring case dictionary

case =   {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0,
            "x_domain%end": L,
            "y_domain%beg": 0,
            "y_domain%end": L/4,
            "z_domain%beg": 0,
            "z_domain%end":  L,
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "cyl_coord": "F",
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": 1,
            "t_step_save": 1,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "bc_z%beg": -1,
            "bc_z%end": -1,
            "viscous": "F",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            # 'prim_vars_wrt'                :'T',
            "parallel_io": "T",
    "chemistry"                    : "T",
    'chem_params%diffusion'        : "F",
    'chem_params%reactions'        : "F",
    'cantera_file'                 : ctfile,

    'prim_vars_wrt'                : 'T',
    'chem_wrt_T'                   : 'T',
        'vel_wrt(1)'                   : 'T',
    'vel_wrt(2)'                   : 'T',
    'pres_wrt'                     : 'T',
    # ========================================
            # I will use 1 for WATER properties, and 2 for AIR properties
            # Patch 1: Background (AIR - 2)
            "patch_icpp(1)%geometry": 13,
            "patch_icpp(1)%hcid"           : 302,
            "patch_icpp(1)%x_centroid": L/2,
            "patch_icpp(1)%y_centroid": L/8,
            "patch_icpp(1)%z_centroid": L/2,
            "patch_icpp(1)%length_x":  L,
            "patch_icpp(1)%length_y": L/4,
            "patch_icpp(1)%length_z":  L,
            "patch_icpp(1)%vel(1)": 0,
            "patch_icpp(1)%vel(2)": 0,
            "patch_icpp(1)%vel(3)": 0,
            "patch_icpp(1)%pres": 1000,
            "patch_icpp(1)%alpha_rho(1)": 1,
            "patch_icpp(1)%alpha(1)": 1,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1),
            "fluid_pp(1)%pi_inf": 0,
        }


if args.chemistry:
      for i in range(len(sol_L.Y)):
        case[f'chem_wrt_Y({1})']    = 'T'
        case[f'patch_icpp(1)%Y({i+1})'] = sol_L.Y[i]
      #  case[f'patch_icpp(2)%Y({i+1})'] = sol_L.Y[i]
    #    #case[f'patch_icpp(3)%Y({i+1})'] = sol_L.Y[i]

if __name__ == '__main__':
    print(json.dumps(case))


