#!/usr/bin/env python3

# References:
# + https://doi.org/10.1016/j.compfluid.2013.10.014: 4.3. Multi-component inert shock tube

import json
import argparse
import math

import cantera as ct

parser = argparse.ArgumentParser(
    prog="1D_MultiCompo_Diffu",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default='{}', metavar="DICT",
                    help="MFC's toolchain's internal state.")
parser.add_argument("--no-chem", dest='chemistry', default=True, action="store_false",
                    help="Disable chemistry.")

args = parser.parse_args()

ctfile    = 'grigri.yaml'
sol_L     = ct.Solution(ctfile)
sol_L.TPX =  300,  8000, 'O2:2,N2:2,H2O:5'

L    = 0.05
Nx   = 100
dx   = L / Nx
dt   = 0.3e-6
Tend = 0.05

NT         = int(Tend / dt)
SAVE_COUNT = 5000
NS         = 5000
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
    't_step_stop'                  : NT,
    't_step_save'                  : NS,
    't_step_print'                 : NS,
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
    'bc_x%beg'                     :-1,
    'bc_x%end'                     :-1,
    'viscous'                      : 'F',
    # ==========================================================================

    # Chemistry ================================================================
    'chemistry'                    : 'T' if not args.chemistry else 'T',
    'chem_params%diffusion'        : 'T',
    'chem_params%reactions'        : 'F',
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                : 'T',
    "chem_wrt_T"                  : "T",
    # ==========================================================================

    # ==========================================================================
    'patch_icpp(1)%geometry'       : 1,
    'patch_icpp(1)%x_centroid'     : L/2,
    'patch_icpp(1)%length_x'       : L,
    'patch_icpp(1)%vel(1)'         : '0',
    'patch_icpp(1)%pres'           : 1.01325e5,
    'patch_icpp(1)%alpha(1)'       : 1,
    'patch_icpp(1)%Y(1)'           : '(0.195-0.142)*(1-0.5*exp(-(x-0.05d0/2.0d0)**2/(2.5d0*10.0d0**(-3.0d0))**2))+0.142',
    'patch_icpp(1)%Y(2)'           : '(0.0-0.1)*(1-0.5*exp(-(x-0.05d0/2.0d0)**2/(2.5d0*10.0d0**(-3.0d0))**2))+0.1',
    'patch_icpp(1)%Y(3)'          : '(0.214-0.0)*(1-0.5*exp(-(x-0.05d0/2.0d0)**2/(2.5d0*10.0d0**(-3.0d0))**2))+0.0',
    'patch_icpp(1)%Y(4)'          : '(0.591-0.758)*(1-0.5*exp(-(x-0.05d0/2.0d0)**2/(2.5d0*10.0d0**(-3.0d0))**2))+0.758',
    'patch_icpp(1)%alpha_rho(1)'   : '1.01325d0*10.0d0**(5.0d0)/(((320.0d0-1350.0d0)*(1.0d0-0.50d0*exp(-(x-0.05d0/2.0d0)**2/(2.5d0*10.0d0**(-3.0d0))**2))+1350.0d0)*8.3144626d0*1000.0d0*( ((0.195d0-0.142d0)*(1.0d0-0.5d0*exp(-(x-0.05d0/2.0d0)**2/(2.5d0*10.0d0**(-3.0d0))**2))+0.142d0)/31.998d0 +((0.0-0.1)*(1-0.5*exp(-(x-0.05d0/2.0d0)**2/(2.5d0*10.0d0**(-3.0d0))**2))+0.1)/18.01508d0+ ((0.214-0.0)*(1-0.5*exp(-(x-0.05d0/2.0d0)**2/(2.5d0*10.0d0**(-3.0d0))**2))+0.0)/16.04256 + ((0.591-0.758)*(1-0.5*exp(-(x-0.05d0/2.0d0)**2/(2.5d0*10.0d0**(-3.0d0))**2))+0.758)/28.0134))',
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.0E+00/(1.9326E+00-1.0E+00),
    'fluid_pp(1)%pi_inf'           : 0,
    # ==========================================================================

    # Chemistry ================================================================
    'cantera_file'                 : ctfile,
    # ==========================================================================
}

if args.chemistry:
    for i in range(4):
        case[f'chem_wrt_Y({i + 1})'] = 'T'
if __name__ == '__main__':
    print(json.dumps(case))
