#!/usr/bin/env python3

import math
import json

Ny = 359
Nx = 269

eps = 1e-6
dt = 1e-5
time_end = 0.0070

Nt = int(time_end/dt)

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 :  -6.,
    'x_domain%end'                 :  6.,
    'y_domain%beg'                 :  0.,
    'y_domain%end'                 :  9.,
    'm'                            : int(Nx),
    'n'                            : int(Ny),
    'p'                            : 0,
    'dt'                           : dt,
    't_step_start'                 : 0,
    't_step_stop'                  : Nt,
    't_step_save'                  : int(Nt/35.),
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 2,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 2,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'T',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'weno_Re_flux'                 : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -6,
    'bc_x%end'                     : -6,
    'bc_y%beg'                     : -6,
    'bc_y%end'                     : -5,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    # ==========================================================================
                                                                
    # Patch 1: Liquid ==========================================================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%y_centroid'     : 4.5,
    'patch_icpp(1)%length_x'       : 12,
    'patch_icpp(1)%length_y'       : 9,
    'patch_icpp(1)%vel(1)'         : 0.,
    'patch_icpp(1)%vel(2)'         : 0.E+00,
    'patch_icpp(1)%pres'           : 1.,
    'patch_icpp(1)%alpha_rho(1)'   : eps*1,
    'patch_icpp(1)%alpha_rho(2)'   : (1-eps)*1000,
    'patch_icpp(1)%alpha(1)'       : eps,
    'patch_icpp(1)%alpha(2)'       : 1-eps,
    # ==========================================================================

    # Patch 2: Air =============================================================
    'patch_icpp(2)%geometry'       : 2,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%x_centroid'     : 0,
    'patch_icpp(2)%y_centroid'     : 6.,
    'patch_icpp(2)%radius'         : 1,
    'patch_icpp(2)%vel(1)'         : 0,
    'patch_icpp(2)%vel(2)'         : 0.E+00,
    'patch_icpp(2)%pres'           : 8290,
    'patch_icpp(2)%alpha_rho(1)'   : (1-eps)*1270,
    'patch_icpp(2)%alpha_rho(2)'   : eps*1000,
    'patch_icpp(2)%alpha(1)'       : (1-eps),
    'patch_icpp(2)%alpha(2)'       : eps,
    # ==========================================================================

       # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(2E+00-1.E+00),
    'fluid_pp(1)%pi_inf'           : 0.,
    'fluid_pp(2)%gamma'            : 1.E+00/(7.15E+00-1.E+00),
    'fluid_pp(2)%pi_inf'           : 3.E+08*7.13E+00/(7.15E+00-1.E+00),
# ==============================================================================
}))

# ==============================================================================
