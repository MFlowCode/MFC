#!/usr/bin/python
# Dependencies and Logistics ===================================================
# Command to navigate between directories
from os import chdir

# Command to acquire directory path
from os.path import dirname

# Command to acquire script name and module search path
from sys import argv, path
import math


ps = 248758.567
gam = 1.4
rho = 1.

c_l = math.sqrt( 1.4*ps/rho )

vel = 230.

leng = 1.


Ny = 100.
Nx = Ny*4
dx = leng/Nx

time_end = 4*leng/vel
cfl = 0.3

dt = cfl * dx/c_l 
Nt = int(time_end/dt)

print 'c_l, Ma', c_l, vel/c_l
print 'Nt', Nt


# Navigating to script directory
if len(dirname(argv[0])) != 0: chdir(dirname(argv[0]))

# Adding master_scripts directory to module search path
mfc_dir = '../../src'; path[:0] = [mfc_dir + '/master_scripts']

# Command to execute the MFC components
from m_python_proxy import f_execute_mfc_component

# Serial or parallel computational engine
engine = 'parallel'
# ==============================================================================


# Case Analysis Configuration ==================================================

# Selecting MFC component
comp_name = argv[1].strip()

# Configuring case dictionary
case_dict =                                                                    \
    {                                                                          \
                    # Logistics ================================================
                    'case_dir'                     : '\'.\'',                  \
                    'run_time_info'                : 'T',                      \
                    'nodes'                        : 1,                        \
                    'ppn'                          : 24,                        \
                    'queue'                        : 'normal',                 \
                    'walltime'                     : '24:00:00',               \
                    'mail_list'                    : '',                       \
                    # ==========================================================
                                                                               \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 :  -leng/2.,                   \
                    'x_domain%end'                 :  leng/2+3*leng,                   \
                    'y_domain%beg'                 :  0.,                  \
                    'y_domain%end'                 :  leng/2.,           \
                    'm'                            : int(Nx),                      \
                    'n'                            : int(Ny),                        \
                    'p'                            : 0,                        \
                    'cyl_coord'                    : 'T',                      \
                    'dt'                           : dt,                   \
                    't_step_start'                 : 0,                        \
                    't_step_stop'                  : Nt,                     \
                    't_step_save'                  : int(Nt/100.),                     \
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
                    'char_decomp'                  : 'F',                      \
                    'mapped_weno'                  : 'T',                      \
                    'null_weights'                 : 'F',                      \
                    'mp_weno'                      : 'F',                      \
                    'weno_avg'                     : 'F',                      \
                    'weno_Re_flux'                 : 'F',                      \
		    'riemann_solver'               : 2,                        \
                    'wave_speeds'                  : 1,                        \
                    'avg_state'                    : 2,                        \
                    'commute_err'                  : 'F',                      \
                    'split_err'                    : 'F',                      \
                    'reg_eps'                      : 1.E+00,                   \
                    'bc_x%beg'                     : -6,                       \
                    'bc_x%end'                     : -6,                       \
                    'bc_y%beg'                     : -2,                      \
                    'bc_y%end'                     : -6,                      \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
		    'parallel_io'                  :'T',                       \
		    # ==========================================================
                                                                                
		    # Patch 1: Background  ============================
                    'patch_icpp(1)%geometry'       : 3,                        \
                    'patch_icpp(1)%x_centroid'     : 0.,                  \
                    'patch_icpp(1)%y_centroid'     : leng*0.25,          \
                    'patch_icpp(1)%length_x'       : 10*leng,                   \
                    'patch_icpp(1)%length_y'       : leng*0.5,         \
                    'patch_icpp(1)%vel(1)'         : vel,                   \
                    'patch_icpp(1)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(1)%pres'           : 101325.,                   \
                    'patch_icpp(1)%alpha_rho(1)'   : 1.29,                \
                    'patch_icpp(1)%alpha_rho(2)'   : 0.E+00,                   \
                    'patch_icpp(1)%alpha(1)'       : 1.E+00,                   \
                    'patch_icpp(1)%alpha(2)'       : 0.E+00,                   \
                    # ==========================================================

		    # Patch 2: Shocked state ============================
                    'patch_icpp(2)%geometry'       : 3,                        \
                    'patch_icpp(2)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(2)%x_centroid'     : -3*leng/8.,                  \
                    'patch_icpp(2)%y_centroid'     : 0.25*leng,          \
                    'patch_icpp(2)%length_x'       : leng/4.,                   \
                    'patch_icpp(2)%length_y'       : leng*0.5,         \
                    'patch_icpp(2)%vel(1)'         : vel,                   \
                    'patch_icpp(2)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(2)%pres'           : ps,                   \
                    'patch_icpp(2)%alpha_rho(1)'   : 2.4,                \
                    'patch_icpp(2)%alpha_rho(2)'   : 0.E+00,                   \
                    'patch_icpp(2)%alpha(1)'       : 1.E+00,                   \
                    'patch_icpp(2)%alpha(2)'       : 0.E+00,                   \
                    # ==========================================================
			
                    # Patch 3: Bubble  ======================================
                    'patch_icpp(3)%geometry'       : 2,                        \
                    'patch_icpp(3)%x_centroid'     : 0.E+00,                 \
                    'patch_icpp(3)%y_centroid'     : 0.E+00,                  \
                    'patch_icpp(3)%radius'         : leng/5.,                  \
                    'patch_icpp(3)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(3)%vel(1)'         : 0.,                   \
                    'patch_icpp(3)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(3)%pres'           : 101325.,                   \
                    'patch_icpp(3)%alpha_rho(1)'   : 0.E+00,                   \
                    'patch_icpp(3)%alpha_rho(2)'   : 0.167,                  \
                    'patch_icpp(3)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(3)%alpha(2)'       : 1.E+00,                   \
                    # ==========================================================
 
		    # Fluids Physical Parameters ===============================
                    'fluid_pp(1)%gamma'            : 1.E+00/(1.4E+00-1.E+00),  \
                    'fluid_pp(1)%pi_inf'           : 0., \
                    'fluid_pp(2)%gamma'            : 1.E+00/(1.6666E+00-1.E+00),  \
                    'fluid_pp(2)%pi_inf'           : 0.E+00,                   \
	            # ==========================================================

    }

# Executing MFC component
f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine)

# ==============================================================================
