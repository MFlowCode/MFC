#!/usr/bin/env python3

# References:
# + https://doi.org/10.1016/j.compfluid.2013.10.014: 4.3. Multi-component inert shock tube

import json
import argparse
import math

import cantera as ct

parser = argparse.ArgumentParser(
    prog="nD_inert_shocktube",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default='{}', metavar="DICT",
                    help="MFC's toolchain's internal state.")
parser.add_argument("--no-chem", dest='chemistry', default=True, action="store_false",
                    help="Disable chemistry.")

args = parser.parse_args()

ctfile    = 'h2o2.yaml'
sol_L     = ct.Solution(ctfile)
sol_L.TPX =  300,  8000, 'O2:2,N2:2,H2O:5'

L    = 0.015
Nx   = 1280
dx   = L / Nx
dt   = 0.5e-8
Tend = 0.60e-3

NT         = int(Tend / dt)
SAVE_COUNT = 800
NS         = 800
case = {
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : 0,
    'x_domain%end'                 : +L,
    'm'                            : Nx,
    'n'                            : 0,
    'p'                            : 0,
    'dt'                           : float(dt),
    't_step_start'                 : 0,
    't_step_stop'                  : 1000,
    't_step_save'                  : 500,
    't_step_print'                 : 1,
    'parallel_io'                  : 'F',
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'model_eqns'                   : 2,
    'num_fluids'                   : 1,
    'num_patches'                  : 1,
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1E-16,
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'T',
    'mp_weno'                      : 'T',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 2,
    'avg_state'                    : 1,
    'bc_x%beg'                     :-7,
    'bc_x%end'                     :-8,
    'viscous'                      : 'T',
    # ==========================================================================

    # Chemistry ================================================================
    'chemistry'                    : 'T' if not args.chemistry else 'T',
    'chem_params%diffusion'        : 'T',
    'chem_params%reactions'        : 'T',
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                : 'T',
    # ==========================================================================

    # ==========================================================================
    'patch_icpp(1)%geometry'       : 15,
    'patch_icpp(1)%hcid'           : 102,
    'patch_icpp(1)%x_centroid'     : L/2,
    'patch_icpp(1)%length_x'       : L,
    'patch_icpp(1)%vel(1)'         : '0',
    'patch_icpp(1)%pres'           : 1.01325e5,
    'patch_icpp(1)%alpha(1)'       : 1,
    'patch_icpp(1)%alpha_rho(1)'   :  '1',
     # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.0E+00/(1.5E+00-1.0E+00),
    'fluid_pp(1)%pi_inf'           : 0,
    'fluid_pp(1)%Re(1)'            : 20000000,
    # ==========================================================================

    # Chemistry ================================================================
    'cantera_file'                 : ctfile,
    # ==========================================================================
}


if __name__ == '__main__':
    print(json.dumps(case))
