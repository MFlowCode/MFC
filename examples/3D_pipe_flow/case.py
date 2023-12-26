#!/usr/bin/env python3
import math
import json
eps = 1e-6
# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================
    'run_time_info'                : 'T',
    # ==========================================================

    # Computational Domain Parameters ==========================
    'x_domain%beg'                 : 0.,
    'x_domain%end'                 : 10.,
    'y_domain%beg'                 : 0.,
    'y_domain%end'                 : 1.,
    'z_domain%beg'                 : 0,
    'z_domain%end'                 : 2*math.pi,
    'cyl_coord'                    : 'T',
    'm'                            : 199,
    'n'                            : 99,
    'p'                            : 79,
    'dt'                           : 2.5E-04,
    't_step_start'                 : 0,
    't_step_stop'                  : 20000,
    't_step_save'                  : 1000,
    # ==========================================================

    # Simulation Algorithm Parameters ==========================
    'num_patches'                  : 1,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 2,
    'adv_alphan'                   : 'F',
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1e-16,
    'mapped_weno'                  : 'T',
    'weno_Re_flux'                 : 'T',
    'mp_weno'                      : 'T',
    'weno_avg'                     : 'T',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -16,
    'bc_x%vb1'                     : 0.1,
    'bc_x%end'                     : -8,
    'bc_y%beg'                     : -14,
    'bc_y%end'                     : -16,
    'bc_z%beg'                     : -1,
    'bc_z%end'                     : -1,
    # ==========================================================

    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    # ==========================================================
                                                            
    # Patch 1: Base ============================================
    'patch_icpp(1)%geometry'       : 9,
    'patch_icpp(1)%x_centroid'     : 5.,
    'patch_icpp(1)%y_centroid'     : 0.,
    'patch_icpp(1)%z_centroid'     : 0,
    'patch_icpp(1)%length_x'       : 10.,
    'patch_icpp(1)%length_y'       : 2.,
    'patch_icpp(1)%length_z'       : 2,
    'patch_icpp(1)%vel(1)'         : 0,
    'patch_icpp(1)%vel(2)'         : 0.,
    'patch_icpp(1)%vel(3)'         : 0,
    'patch_icpp(1)%pres'           : 5,
    'patch_icpp(1)%alpha_rho(1)'   : 0.5,
    'patch_icpp(1)%alpha(1)'       : 0.5,
    'patch_icpp(1)%alpha_rho(2)'   : 0.5,
    'patch_icpp(1)%alpha(2)'       : 0.5,
    # ==========================================================

    # Fluids Physical Parameters ===============================
    'fluid_pp(1)%gamma'            : 1./(1.4-1.),
    'fluid_pp(1)%pi_inf'           : 0.,
    'fluid_pp(1)%Re(1)'            : 1e4,
    'fluid_pp(2)%gamma'            : 1./(1.4-1.),
    'fluid_pp(2)%pi_inf'           : 0.E+00,
    'fluid_pp(2)%Re(1)'            : 1e4,
    # ==========================================================
}))
# ==============================================================================