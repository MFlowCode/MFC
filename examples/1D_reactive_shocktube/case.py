#!/usr/bin/env python3

# References:
# + https://doi.org/10.1016/j.ijhydene.2023.03.190:  Verification of numerical method
# + https://doi.org/10.1016/j.compfluid.2013.10.014: 4.7. Multi-species reactive shock tube

import json
import cantera as ct

ctfile    = 'h2o2.yaml'
sol_L     = ct.Solution(ctfile)
sol_L.DPX =   0.072,  7173, 'H2:2,O2:1,AR:7'

sol_R     = ct.Solution(ctfile)
sol_R.DPX = 0.18075, 35594, 'H2:2,O2:1,AR:7'

u_l = 0
u_r = -487.34

L  = 0.12
Nx = 800
dx = L/Nx
dt = dx/abs(u_r)*0.01
Tend=230e-6

NT=int(Tend/dt)
SAVE_COUNT=100
NS=NT//SAVE_COUNT

chemistry = True

case = {
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : 0,
    'x_domain%end'                 : L,
    'm'                            : Nx,
    'n'                            : 0,
    'p'                            : 0,
    'dt'                           : float(dt),
    't_step_start'                 : 0,
    't_step_stop'                  : NT,
    't_step_save'                  : NS,
    't_step_print'                 : NS,
    'parallel_io'                  : 'F',

    # Simulation Algorithm Parameters ==========================================
    'model_eqns'                   : 2,
    'num_fluids'                   : 1,
    'num_patches'                  : 2,
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1E-16,
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'T',
    'mp_weno'                      : 'T',
    'riemann_solver'               : 1,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     :-2,
    'bc_x%end'                     :-3,

    # Chemistry ================================================================
    'chemistry'                    : 'F' if not chemistry else 'T',
    'chem_params%advection'        : 'T',
    'chem_params%diffusion'        : 'F',
    'chem_params%reactions'        : 'T',
    'cantera_file'                 : ctfile,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                : 'T',
    # ==========================================================================

    # ==========================================================================
    'patch_icpp(1)%geometry'       : 1,
    'patch_icpp(1)%x_centroid'     : L/4,
    'patch_icpp(1)%length_x'       : L/2,
    'patch_icpp(1)%vel(1)'         : u_l,
    'patch_icpp(1)%pres'           : sol_L.P,
    'patch_icpp(1)%alpha(1)'       : 1,
    'patch_icpp(1)%alpha_rho(1)'   : sol_L.density,
    # ==========================================================================

    # ==========================================================================
    'patch_icpp(2)%geometry'       : 1,
    'patch_icpp(2)%x_centroid'     : 3*L/4,
    'patch_icpp(2)%length_x'       : L/2,
    'patch_icpp(2)%vel(1)'         : u_r,
    'patch_icpp(2)%pres'           : sol_R.P,
    'patch_icpp(2)%alpha(1)'       : 1,
    'patch_icpp(2)%alpha_rho(1)'   : sol_R.density,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.0E+00/(4.4E+00-1.0E+00),
    'fluid_pp(1)%pi_inf'           : 0,
    # ==========================================================================
}

if chemistry:
    for i in range(len(sol_L.Y)):
        case[f'patch_icpp(1)%Y({i+1})'] = sol_L.Y[i]
        case[f'patch_icpp(2)%Y({i+1})'] = sol_R.Y[i]

if __name__ == '__main__':
    print(json.dumps(case))
