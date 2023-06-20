#!/usr/bin/env python3

#================================================================================
#
# In order to run this case uncomment the following lines in the subroutine
# s_1D_analytical() in src/pre_process/m_patches.f90
#
#    ! Shu-Osher Problem
#    if (x_cc(i) > -4d0) then
#       q_prim_vf(contxb)%sf(i, 0, 0) = 1 + 0.2*sin(5*x_cc(i))
#       q_prim_vf(momxb)%sf(i, 0, 0) = 0d0
#       q_prim_vf(E_idx)%sf(i, 0, 0) = 1d0
#    end if
#
#================================================================================

import math
import json

# Numerical setup
Nx = 1000
dx = 1./(1.*(Nx+1))

Tend, Nt = 1.8, 20000
mydt     = Tend/(1.*Nt)

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -5.,
    'x_domain%end'                 : 5. ,
    'm'                            : Nx,
    'n'                            : 0,
    'p'                            : 0,
    'dt'                           : mydt,
    't_step_start'                 : 0,
    't_step_stop'                  : int(Nt),
    't_step_save'                  : int(math.ceil(Nt/10.)),
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 1,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -3,
    'bc_x%end'                     : -3,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    # ==========================================================================
                                                            
    # Patch 1 L ================================================================
    'patch_icpp(1)%geometry'       : 15,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%length_x'       : 10.,
    'patch_icpp(1)%vel(1)'         : 2.629369,
    'patch_icpp(1)%pres'           : 10.3333,
    'patch_icpp(1)%alpha_rho(1)'   : 3.957143,
    'patch_icpp(1)%alpha(1)'       : 1.,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(1.4-1.E+00),
    'fluid_pp(1)%pi_inf'           : 0.0,
    # ==========================================================================
}))

# ==============================================================================
