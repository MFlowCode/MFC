#!/usr/bin/env python3

import json
import math

# Dynamic Viscosity
Mu1 = 0.0000184
#Mu2 = 0.01
rho1 = 0.2199
gam_a = 1.4


# Patch Design
D = 0.01

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -6*D,
    'x_domain%end'                 : 6*D,
    'y_domain%beg'                 : -3*D,
    'y_domain%end'                 : 3*D,
    'm'                            : 399,
    'n'                            : 199,
    'p'                            : 0,
    'dt'                           : 2.0E-07,
    't_step_start'                 : 0,
    't_step_stop'                  : 3000,
    't_step_save'                  : 30,
    # ==========================================================================
    
    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 1,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'F',
    'weno_avg'                     : 'F',
    'avg_state'                    : 2,
    'mapped_weno'                  : 'F',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'viscous'                      : 'T',
    'bc_x%beg'                     : -3,
    'bc_x%end'                     : -3,
    'bc_y%beg'                     : -3,
    'bc_y%end'                     : -3,
    'ib'		                   : 'T',
    'num_ibs'                      :  1,
    'perturb_flow'                 : 'T',
    'perturb_flow_fluid'           :  1,
    'perturb_flow_mag'             : 0.0001,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    # Export primitive variables in double precision with parallel
    # I/O to minimize I/O computational time during large simulations
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    #'omega_wrt(1)'                :'T',  
    #'omega_wrt(2)'                :'T',  
    'omega_wrt(3)'                 :'T',  
    'fd_order'                     : 2, 
    # ==========================================================================

    # Patch: Middle ============================================================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : 0, 
    'patch_icpp(1)%y_centroid'     : 0,
    'patch_icpp(1)%length_x'       : 1000*D,
    'patch_icpp(1)%length_y'       : 1000*D,
    'patch_icpp(1)%vel(1)'         : 527.2E+00,
    'patch_icpp(1)%vel(2)'         : 0.0E+00,
    'patch_icpp(1)%pres'           : 10918.2549,
    'patch_icpp(1)%alpha_rho(1)'   : (1.0)*rho1,
    'patch_icpp(1)%alpha(1)'       : 1.0,
    # ==========================================================================
    'patch_ib(1)%geometry'              : 5,
    'patch_ib(1)%model%filepath'        : 'RootToMFC/MFC/examples/2D_IBM_STL_wedge/20deg_wedge.stl',
    'patch_ib(1)%model%translate(1)'    : -0.05 + 0.0424444350,
    'patch_ib(1)%model%translate(2)'    : -0.05 + 0.0445000000,
    'patch_ib(1)%model%spc'             : 200,
    'patch_ib(1)%model%threshold'       : 0.9,
    'patch_ib(1)%slip'                  : 'F',
    # # ========================================================================
   
    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(gam_a-1.E+00),  # 2.50(Not 1.40)
    'fluid_pp(1)%pi_inf'           : 0,
    'fluid_pp(1)%Re(1)'            : 7535533.2,
    # ==========================================================================
}))

