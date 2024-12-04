#!/usr/bin/env python3

# References:
# + https://doi.org/10.1016/j.ijhydene.2023.03.190:  Verification of numerical method
# + https://doi.org/10.1016/j.compfluid.2013.10.014: 4.7. Multi-species reactive shock tube

import json, argparse
import cantera as ct

parser = argparse.ArgumentParser(
    prog="1D_reactive_shocktube",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default='{}', metavar="DICT",
                    help="MFC's toolchain's internal state.")
parser.add_argument("--no-chem", dest='chemistry', default=True, action="store_false",
                                 help="Disable chemistry.")
parser.add_argument("--scale",   type=float,  default=1,    help="Scale.")

args = parser.parse_args()

ctfile    = 'h2o2.yaml'
sol_L     = ct.Solution(ctfile)
sol_L.DPX =   0.072,  7173, 'H2:2,O2:1,AR:7'

sol_R     = ct.Solution(ctfile)
sol_R.DPX = 0.18075, 35594, 'H2:2,O2:1,AR:7'

u_l = 0
u_r = -487.34

L  = 0.12
NScale = 2
Nx = 1 * NScale * 400 * args.scale
Ny = 1 * NScale * 200 * args.scale
dx = L/Nx
dt = dx/abs(u_r)*0.05*0.1*0.2
Tend=230e-6/2

NT=int(Tend/dt)
SAVE_COUNT=100
NS=NT//SAVE_COUNT

case = {
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : 0,
    'x_domain%end'                 : L,
    'y_domain%beg'                 : 0,
    'y_domain%end'                 : L/2,
    'm'                            : Nx,
    'n'                            : Ny,
    'p'                            : 0,
    'dt'                           : float(dt),
    't_step_start'                 : 0,
    't_step_stop'                  : NT,
    't_step_save'                  : NS,
    't_step_print'                 : 10,
    'parallel_io'                  : 'F' if args.mfc.get("mpi", True) else 'F',

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
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     :-3,
    'bc_x%end'                     :-3,
    'bc_y%beg'                     :-1,
    'bc_y%end'                     :-1,

    # Chemistry ================================================================
    'chemistry'                    : 'F' if not args.chemistry else 'T',
    'chem_params%diffusion'        : 'F',
    'chem_params%reactions'        : 'T',
    'cantera_file'                 : ctfile,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    #'prim_vars_wrt'                : 'T',
    'chem_wrt_T'                   : 'F',
    # ==========================================================================
    'rho_wrt' : 'T',

    # ==========================================================================
    'patch_icpp(1)%geometry'       : 7,
    'patch_icpp(1)%x_centroid'     : L/2,
    'patch_icpp(1)%y_centroid'     : L/4,
    'patch_icpp(1)%length_x'       : L,
    'patch_icpp(1)%length_y'       : L/2,
    'patch_icpp(1)%hcid'           : 666,
    'patch_icpp(1)%vel(1)'         : 0,
    'patch_icpp(1)%vel(2)'         : 0,
    'patch_icpp(1)%pres'           : sol_L.P,
    'patch_icpp(1)%alpha(1)'       : 1,
    'patch_icpp(1)%alpha_rho(1)'   : sol_L.density,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.0E+00/(4.4E+00-1.0E+00),
    'fluid_pp(1)%pi_inf'           : 0,
    # ==========================================================================
}

if args.chemistry:
    for i in range(len(sol_L.Y)):
        #if sol_L.species_name(i) in ['H2', 'O2', 'OH', 'H2O']:
        #    case[f'chem_wrt_Y({i + 1})'] = 'T'
        #else:
        #    case[f'chem_wrt_Y({i + 1})'] = 'F'
        case[f'patch_icpp(1)%Y({i+1})'] = sol_L.Y[i]

if __name__ == '__main__':
    print(json.dumps(case))
