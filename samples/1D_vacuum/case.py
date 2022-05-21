#!/usr/bin/env python3

import json

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'case_dir'                     : '\'.\'',
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : 0.E+00,
    'x_domain%end'                 : 1.E+00,
    'm'                            : 199,
    'n'                            : 0,
    'p'                            : 0,
    'dt'                           : 5.E-08,
    't_step_start'                 : 0,
    't_step_stop'                  : 15000,
    't_step_save'                  : 1000,
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 2,
    'model_eqns'                   : 3,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 2,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'T',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_vars'                    : 2,
    'weno_order'                   : 3,
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
                                            
    # Patch 1: Left state ======================================================
    'patch_icpp(1)%geometry'       : 1,
    'patch_icpp(1)%x_centroid'     : 0.25E+00,
    'patch_icpp(1)%length_x'       : 0.5E+00,
    'patch_icpp(1)%vel(1)'         : -100.E+00,
    'patch_icpp(1)%pres'           : 1.E+05,
    'patch_icpp(1)%alpha_rho(1)'   : 1000.E+00*0.99E+00,
    'patch_icpp(1)%alpha_rho(2)'   : 10.E+00*0.01E+00,
    'patch_icpp(1)%alpha(1)'       : 0.99E+00,
    'patch_icpp(1)%alpha(2)'       : 0.01E+00,
    # ==========================================================================

    # Patch 2: Right state =====================================================
    'patch_icpp(2)%geometry'       : 1,
    'patch_icpp(2)%x_centroid'     : 0.75E+00,
    'patch_icpp(2)%length_x'       : 0.5E+00,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%vel(1)'         : 100.E+00,
    'patch_icpp(2)%pres'           : 1.E+05,
    'patch_icpp(2)%alpha_rho(1)'   : 1000.E+00*0.99E+00,
    'patch_icpp(2)%alpha_rho(2)'   : 10.E+00*0.01E+00,
    'patch_icpp(2)%alpha(1)'       : 0.99E+00,
    'patch_icpp(2)%alpha(2)'       : 0.01E+00,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(4.4E+00-1.E+00),
    'fluid_pp(1)%pi_inf'           : 4.4E+00*6.E+08/(4.4E+00-1.E+00),
    'fluid_pp(2)%gamma'            : 1.E+00/(1.4E+00-1.E+00),
    'fluid_pp(2)%pi_inf'           : 0.E+00,
    # ==========================================================================
}))

# ==============================================================================
