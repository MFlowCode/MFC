#!/usr/bin/python
# Dependencies and Logistics ===================================================
# Command to navigate between directories
from os import chdir

# Command to acquire directory path
from os.path import dirname

# Command to acquire script name and module search path
from sys import argv, path
import math

myeps=1.4/150.

p_l = 1E+06
p_g = 1E+06

rho_l = 1000.
rho_g = 1.

v_l = 500.
v_g = -500.

c_l = math.sqrt( ( p_l+4.E8)/rho_l )
c_g = math.sqrt( 1.4*p_g/rho_g )

print 'c_l, Mach_l', c_l, v_l/c_l
print 'c_g, Mach_g', c_g, v_g/c_g


# Navigating to script directory
if len(dirname(argv[0])) != 0: chdir(dirname(argv[0]))

# Adding master_scripts directory to module search path
mfc_dir = '../../src'; path[:0] = [mfc_dir + '/master_scripts']

# Command to execute the MFC components
from m_python_proxy import f_execute_mfc_component

# Serial or parallel computational engine
engine = 'serial'
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
                    'ppn'                          : 1,                        \
                    'queue'                        : 'normal',                 \
                    'walltime'                     : '24:00:00',               \
                    'mail_list'                    : '',                       \
                    # ==========================================================
                                                                               \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 : -0.5,                   \
                    'x_domain%end'                 : 0.5,                   \
                    'y_domain%beg'                 : -0.5,                  \
                    'y_domain%end'                 : 0.5,           \
                    'cyl_coord'                    : 'F',                      \
                    'm'                            : 50,                      \
                    'n'                            : 50,                        \
                    'p'                            : 0,                        \
                    'dt'                           : 5.E-10,                   \
                    't_step_start'                 : 0,                        \
                    't_step_stop'                  : 10,                     \
                    't_step_save'                  : 1,                     \
		            # ==========================================================
                                                                               \
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 2,                        \
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
                    'weno_Re_flux'                 : 'T',                      \
		    'riemann_solver'               : 2,                        \
                    'wave_speeds'                  : 1,                        \
                    'avg_state'                    : 2,                        \
                    'commute_err'                  : 'F',                      \
                    'split_err'                    : 'F',                      \
                    'regularization'               : 'F',                      \
                    'reg_eps'                      : 1.E+00,                   \
                    'We_riemann_flux'              : 'F',                      \
                    'We_rhs_flux'                  : 'F',                      \
                    'We_src'                       : 'F',                      \
                    'We_wave_speeds'               : 'F',                      \
                    'lsq_deriv'                    : 'F',                      \
                    'alt_crv'                      : 'F',                      \
                    'bc_x%beg'                     : -1,                       \
                    'bc_x%end'                     : -1,                       \
                    'bc_y%beg'                     : -6,                      \
                    'bc_y%end'                     : -6,                      \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
		    'parallel_io'                  :'T',                       \
		    # ==========================================================
                                                                                
		    # Patch 1: Top fluid, water ============================
                    'patch_icpp(1)%geometry'       : 3,                        \
                    'patch_icpp(1)%x_centroid'     : 0.,                  \
                    'patch_icpp(1)%y_centroid'     : 0,          \
                    'patch_icpp(1)%length_x'       : 1.E+00,                   \
                    'patch_icpp(1)%length_y'       : 1.,         \
                    'patch_icpp(1)%vel(1)'         : v_l,                   \
                    'patch_icpp(1)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(1)%pres'           : p_l,                   \
                    'patch_icpp(1)%alpha_rho(1)'   : rho_l,                \
                    'patch_icpp(1)%alpha_rho(2)'   : rho_l,                   \
                    'patch_icpp(1)%alpha(1)'       : 0.5E+00,                   \
                    'patch_icpp(1)%alpha(2)'       : 0.5E+00,                   \
                    # ==========================================================

                    # Patch 2: Main bottom fluid, water ======================================
                    'patch_icpp(2)%geometry'       : 3,                        \
                    'patch_icpp(2)%x_centroid'     : 0.,                 \
                    'patch_icpp(2)%y_centroid'     : 0.25,                  \
                    'patch_icpp(2)%length_x'       : 1.E+00,                  \
                    'patch_icpp(2)%length_y'       : 0.5,                  \
                    'patch_icpp(2)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(2)%vel(1)'         : v_g,                   \
                    'patch_icpp(2)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(2)%pres'           : p_g,                   \
                    'patch_icpp(2)%alpha_rho(1)'   : 0.,                   \
                    'patch_icpp(2)%alpha_rho(2)'   : rho_g,                  \
                    'patch_icpp(2)%alpha(1)'       : 0.,                   \
                    'patch_icpp(2)%alpha(2)'       : 1.E+00,                   \
                    # ==========================================================

                    # Fluids Physical Parameters ===============================
                    'fluid_pp(1)%gamma'            : 1.E+00/(4.4E+00-1.E+00),  \
                    'fluid_pp(2)%gamma'            : 1.E+00/(4.4E+00-1.E+00),  \
                    'fluid_pp(1)%pi_inf'           : 4.4E+00*6.E+08/(4.4E+00-1.E+00), \
                    'fluid_pp(2)%pi_inf'           : 4.4E+00*6.E+08/(4.4E+00-1.E+00), \
	            'fluid_pp(1)%Re(1)'            : 0.0001,      \
	            'fluid_pp(1)%Re(2)'            : 0.0001,      \
 	            'fluid_pp(2)%Re(1)'            : 0.0001,      \
	            'fluid_pp(2)%Re(2)'            : 0.0001,      \
                    # ==========================================================

    }

# Executing MFC component
f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine)

# ==============================================================================
