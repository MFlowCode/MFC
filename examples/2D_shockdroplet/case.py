#!/usr/bin/env python3
import math
import json

Ma = 2.4
ps = 664016.5 
rho_post_a=3.757918216
rho_a=1.17
rho_w=1000
gam_a = 1.4
gam_w=6.12
pi_w=3.43E8
vel = 575.4980523
rho = 1
c_l = math.sqrt( 1.4*ps/rho )

D = 0.022
Ny = 740.
Nx =5000.
dx = 0.25/Nx #8.3e-6

time_end = 0.0002#50us
cfl = 0.1

dt = cfl * dx/c_l #5.3E-9
Nt = int(time_end/dt)#10000

print(json.dumps({
                    # Logistics ================================================
                    'case_dir'                     : '\'.\'',                  \
                    'run_time_info'                : 'F',                      \
                    # ==========================================================
                                                                               \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 :  0,                   \
                    'x_domain%end'                 :  0.25,                   \
                    'y_domain%beg'                 :  0,                  \
                    'y_domain%end'                 :  0.037,          
                    'm'                            : int(Nx),                      \
                    'n'                            : int(Ny),                        \
                    'p'                            : 0,                        \
                    'dt'                           : dt,                   \
                    't_step_start'                 : 0,                        \
                    't_step_stop'                  : 50000,                     \
                    't_step_save'                  : 1000,#(Nt/1000),                     \
		    # ==========================================================
                                                                               \
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 3,                        \
                    'model_eqns'                   : 2,                        \
                    'alt_soundspeed'               : 'F',                      \
                    'num_fluids'                   : 2,                        \
		    'adv_alphan'                   : 'T',                      \
		    'mpp_lim'                      : 'T',                      \
		    'mixture_err'                  : 'T',                      \
		    'time_stepper'                 : 3,                        \
                    'weno_vars'                    : 2,                        \
                    'weno_order'                   : 5,                        \
                    'weno_eps'                     : 1.E-16,                   \
                    'mapped_weno'                  : 'T',                      \
                    'null_weights'                 : 'F',                      \
                    'mp_weno'                      : 'F',                      \
                    'weno_Re_flux'                 : 'F',                      \
		    'riemann_solver'               : 2,                        \
                    'wave_speeds'                  : 1,                        \
                    'avg_state'                    : 2,                        \
                    'bc_x%beg'                     : -6,#11,                       \
                    'bc_x%end'                     : -6,#12                       \
                    'bc_y%beg'                     : -2,                      \
                    'bc_y%end'                     : -6,                      \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
                    'schlieren_wrt'                 : 'T',                    \
                    'schlieren_alpha(1)'                 : 40,                 \
                    'schlieren_alpha(2)'                 : 400,                    \
                    'fd_order'                 : 4,                    \
   
		    'parallel_io'                  :'T',                       \
		    # ==========================================================
                                                                                
		    # Patch 1: Background  ============================
                    'patch_icpp(1)%geometry'       : 3,                        \
                    'patch_icpp(1)%x_centroid'     : 0.25/2,                  \
                    'patch_icpp(1)%y_centroid'     : 0.037/2,          \
                    'patch_icpp(1)%length_x'       : 0.25,                   \
                    'patch_icpp(1)%length_y'       : 0.037,         \
                    'patch_icpp(1)%vel(1)'         : 0.,                   \
                    'patch_icpp(1)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(1)%pres'           : 101325.,                   \
                    'patch_icpp(1)%alpha_rho(1)'   : 0.,                \
                    'patch_icpp(1)%alpha_rho(2)'   : 1.17,                   \
                    'patch_icpp(1)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(1)%alpha(2)'       : 1.E+00,                   \
                    # ==========================================================

		    # Patch 2: Shocked state ============================
                    'patch_icpp(2)%geometry'       : 3,                        \
                    'patch_icpp(2)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(2)%x_centroid'     : 0.,                  \
                    'patch_icpp(2)%y_centroid'     : 0.037/2,          \
                    'patch_icpp(2)%length_x'       : 0.25-D,                   \
                    'patch_icpp(2)%length_y'       : 0.037,         \
                    'patch_icpp(2)%vel(1)'         : vel,                   \
                    'patch_icpp(2)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(2)%pres'           : ps,                   \
                    'patch_icpp(2)%alpha_rho(1)'   : 0.E+00,                \
                    'patch_icpp(2)%alpha_rho(2)'   : rho_post_a,                   \
                    'patch_icpp(2)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(2)%alpha(2)'       : 1.E+00,                   \
                    # ==========================================================
			
                    # Patch 3: Bubble  ======================================
                    'patch_icpp(3)%geometry'       : 2,                        \
                    'patch_icpp(3)%x_centroid'     : 0.25/2,                 \
                    'patch_icpp(3)%y_centroid'     : 0,                  \
                    'patch_icpp(3)%radius'         : D/2,                  \
                    'patch_icpp(3)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(3)%vel(1)'         : 0.,                   \
                    'patch_icpp(3)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(3)%pres'           : 101325.,                   \
                    'patch_icpp(3)%alpha_rho(1)'   : rho_w,                   \
                    'patch_icpp(3)%alpha_rho(2)'   : 0,#0.061,                  \
                    'patch_icpp(3)%alpha(1)'       : 1.,# 0.95                  \
                    'patch_icpp(3)%alpha(2)'       : 0.,#0.05,                   \
    
                    # ==========================================================
 
		    # Fluids Physical Parameters ===============================
                    'fluid_pp(1)%gamma'            : 1.E+00/(gam_w-1.E+00),  \
                    'fluid_pp(1)%pi_inf'           : pi_w*gam_w/(gam_w-1.E+00), \
                    'fluid_pp(2)%gamma'            : 1.E+00/(gam_a-1.E+00),  \
                    'fluid_pp(2)%pi_inf'           : 0.E+00,                   \
	            # ==========================================================

    }))

