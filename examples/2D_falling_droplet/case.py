#!/usr/bin/env python3
import math
import json

Ny = 1599.
Nx =399.
dx = 0.25/Nx #8.3e-6

eps = 1e-6

time_end = 0.1#50us
cfl = 0.8

c_l = math.sqrt(6.12*(101325 + 3.43e8)/1000)

dt = cfl * dx/c_l #5.3E-9
Nt = int(time_end/dt)#10000

print(json.dumps({
    # Logistics ================================================
    'run_time_info'                : 'F',
    # ==========================================================

    # Computational Domain Parameters ==========================
    'x_domain%beg'                 :  -0.125,
    'x_domain%end'                 :  0.125,
    'y_domain%beg'                 :  0,
    'y_domain%end'                 :  1,
    'm'                            : int(Nx),
    'n'                            : int(Ny),
    'p'                            : 0,
    'dt'                           : dt,
    't_step_start'                 : 0,
    't_step_stop'                  : 100000,
    't_step_save'                  : 1000, #math.ceil(Nt/100),
    # ==========================================================

    # Simulation Algorithm Parameters ==========================
    'num_patches'                  : 2,
    'ib'                           : 'T',
    'num_ibs'                      : 1,
    'model_eqns'                   : 3,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 2,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'T',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'F',
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -2,#11,
    'bc_x%end'                     : -2,#12
    'bc_y%beg'                     : -15,
    'bc_y%end'                     : -15,
    # ==========================================================

    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    # ==========================================================

    # Patch 1: Background  ============================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : 0,
    'patch_icpp(1)%y_centroid'     : 0.5,
    'patch_icpp(1)%length_x'       : 0.25,
    'patch_icpp(1)%length_y'       : 1,
    'patch_icpp(1)%vel(1)'         : 0.,
    'patch_icpp(1)%vel(2)'         : 0.,
    'patch_icpp(1)%pres'           : 101325.,
    'patch_icpp(1)%alpha_rho(1)'   : eps*1000,
    'patch_icpp(1)%alpha_rho(2)'   : (1-eps)*1,
    'patch_icpp(1)%alpha(1)'       : eps,
    'patch_icpp(1)%alpha(2)'       : 1-eps,
    'patch_icpp(1)%cf_val'         : 0,
    # ==========================================================

    # Patch 3: Droplet  ======================================
    'patch_icpp(2)%geometry'       : 2,
    'patch_icpp(2)%smoothen'       : 'T',
    'patch_icpp(2)%smooth_patch_id' : 1,
    'patch_icpp(2)%smooth_coeff'   : 0.5,
    'patch_icpp(2)%x_centroid'     : 0.,
    'patch_icpp(2)%y_centroid'     : 0.75,
    'patch_icpp(2)%radius'          : 0.05,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%vel(1)'         : 0.,
    'patch_icpp(2)%vel(2)'         : 0.,
    'patch_icpp(2)%pres'           : 101325.,
    'patch_icpp(2)%alpha_rho(1)'   : (1-eps)*1000,
    'patch_icpp(2)%alpha_rho(2)'   : eps*1,
    'patch_icpp(2)%alpha(1)'       : 1 - eps,# 0.95
    'patch_icpp(2)%alpha(2)'       : eps,#0.05,
    'patch_icpp(2)%cf_val'         : 1,
    # ==========================================================

    'bf_y'                         : 'T',
    'k_y'                          : 0,
    'w_y'                          : 0,
    'p_y'                          : 0,
    'g_y'                          : 196,

    'sigma'                        : 8.0,

    # IB_Patch Cylinder ========================================
    'patch_ib(1)%geometry'         : 2,
    'patch_ib(1)%x_centroid'       : 0,
    'patch_ib(1)%y_centroid'       : 0.25,
    'patch_ib(1)%radius'           : 0.05,
    'patch_ib(1)%slip'             : 'T',
    # ==========================================================

    # Fluids Physical Parameters ===============================
    'fluid_pp(1)%gamma'            : 1.E+00/(6.12-1.E+00),
    'fluid_pp(1)%pi_inf'           : 3.43e8*6.12/(6.12-1.E+00),
    'fluid_pp(2)%gamma'            : 1.E+00/(1.4-1.E+00),
    'fluid_pp(2)%pi_inf'           : 0.E+00,
    # ==========================================================
}))
