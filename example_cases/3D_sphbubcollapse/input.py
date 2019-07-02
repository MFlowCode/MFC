#!/usr/bin/python

# Dependencies and Logistics ===================================================

# Command to navigate between directories
from os import chdir

# Command to acquire directory path
from os.path import dirname

# Command to acquire script name and module search path
from sys import argv, path

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
                    'ppn'                          : 24,                       \
                    'queue'                        : 'normal',                   \
                    'walltime'                     : '120:00:00',              \
                    'mail_list'                    : '',                       \
                    # ==========================================================
                                                                               \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 : 0.E+00,                   \
                    'x_domain%end'                 : 4.E-03/1.E-03,            \
                    'y_domain%beg'                 : 0.E+00,                   \
                    'y_domain%end'                 : 4.E-03/1.E-03,            \
                    'z_domain%beg'                 : 0.E+00,                   \
                    'z_domain%end'                 : 4.E-03/1.E-03,            \
                    'stretch_x'                    : 'T',                      \
		    'a_x'                          : 4.E+00,                   \
		    'x_a'                          : -1.5E-03/1.E-03,          \
		    'x_b'                          : 1.5E-03/1.E-03,           \
		    'stretch_y'                    : 'T',                      \
		    'a_y'                          : 4.E+00,                   \
		    'y_a'                          : -1.5E-03/1.E-03,          \
		    'y_b'                          : 1.5E-03/1.E-03,           \
		    'stretch_z'                    : 'T',                      \
		    'a_z'                          : 4.E+00,                   \
		    'z_a'                          : -1.5E-03/1.E-03,          \
		    'z_b'                          : 1.5E-03/1.E-03,           \
                    'cyl_coord'                    : 'F',                      \
                    'm'                            : 99,                       \
                    'n'                            : 99,                       \
                    'p'                            : 99,                       \
                    'dt'                           : 1.5E-09/1.E-03,           \
                    't_step_start'                 : 0,                        \
                    't_step_stop'                  : 133300,                   \
                    't_step_save'                  : 100,                      \
		    # ==========================================================
                                                                               \
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 2,                        \
                    'model_eqns'                   : 3,                        \
                    'alt_soundspeed'               : 'F',                      \
                    'num_fluids'                   : 2,                        \
		    'adv_alphan'                   : 'T',                      \
		    'mpp_lim'                      : 'T',                      \
		    'mixture_err'                  : 'T',                      \
		    'time_stepper'                 : 3,                        \
                    'weno_vars'                    : 2,                        \
                    'weno_order'                   : 5,                        \
                    'weno_eps'                     : 1.E-16,                   \
                    'char_decomp'                  : 'T',                      \
                    'avg_state'                    : 2,                        \
                    'mapped_weno'                  : 'T',                      \
                    'null_weights'                 : 'F',                      \
                    'mp_weno'                      : 'F',                      \
                    'weno_avg'                     : 'F',                      \
                    'weno_Re_flux'                 : 'F',                      \
		    'riemann_solver'               : 2,                        \
                    'wave_speeds'                  : 1,                        \
                    'commute_err'                  : 'F',                      \
                    'split_err'                    : 'F',                      \
                    'regularization'               : 'F',                      \
                    'reg_eps'                      : 1.E+00,                   \
                    'anti_diffusion'               : 'F',                      \
                    'We_riemann_flux'              : 'F',                      \
                    'We_rhs_flux'                  : 'F',                      \
                    'We_src'                       : 'F',                      \
                    'We_wave_speeds'               : 'F',                      \
                    'lsq_deriv'                    : 'F',                      \
                    'alt_crv'                      : 'F',                      \
                    'bc_x%beg'                     : -2,                       \
                    'bc_x%end'                     : -6,                       \
                    'bc_y%beg'                     : -2,                       \
                    'bc_y%end'                     : -6,                       \
                    'bc_z%beg'                     : -2,                       \
                    'bc_z%end'                     : -6,                       \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
		    'parallel_io'                  :'T',                       \
		    # ==========================================================
                                                                                
		    # Patch 1: High pressured water ============================
                    'patch_icpp(1)%geometry'       : 9,                        \
                    'patch_icpp(1)%x_centroid'     : 80.E-03/1.E-03,           \
                    'patch_icpp(1)%y_centroid'     : 80.E-03/1.E-03,           \
                    'patch_icpp(1)%z_centroid'     : 80.E-03/1.E-03,           \
                    'patch_icpp(1)%length_x'       : 160.E-03/1.E-03,          \
                    'patch_icpp(1)%length_y'       : 160.E-03/1.E-03,          \
                    'patch_icpp(1)%length_z'       : 160.E-03/1.E-03,          \
                    'patch_icpp(1)%vel(1)'         : 0.E+00,                   \
                    'patch_icpp(1)%vel(2)'         : 0.E+00,                   \
                    'patch_icpp(1)%vel(3)'         : 0.E+00,                   \
                    'patch_icpp(1)%pres'           : 1.E+05,                   \
                    'patch_icpp(1)%alpha_rho(1)'   : 1000.E+00,                \
                    'patch_icpp(1)%alpha_rho(2)'   : 0.E+00,                   \
                    'patch_icpp(1)%alpha(1)'       : 1.E+00,                   \
                    'patch_icpp(1)%alpha(2)'       : 0.E+00,                   \
                    # ==========================================================

                    # Patch 3: Air bubble ======================================
                    'patch_icpp(2)%geometry'       : 8,                        \
                    'patch_icpp(2)%smoothen'       : 'T',                      \
                    'patch_icpp(2)%smooth_patch_id' : 1,                       \
                    'patch_icpp(2)%smooth_coeff'   : 0.5E+00,                  \
                    'patch_icpp(2)%x_centroid'     : 0.E+00,                   \
                    'patch_icpp(2)%y_centroid'     : 0.E+00,                   \
                    'patch_icpp(2)%z_centroid'     : 0.E+00,                   \
                    'patch_icpp(2)%radius'         : 1.E-03/1.E-03,            \
                    'patch_icpp(2)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(2)%vel(1)'         : 0.E+00,                   \
                    'patch_icpp(2)%vel(2)'         : 0.E+00,                   \
                    'patch_icpp(2)%vel(3)'         : 0.E+00,                   \
                    'patch_icpp(2)%pres'           : 1.E+03,                   \
                    'patch_icpp(2)%alpha_rho(1)'   : 0.E+00,                   \
                    'patch_icpp(2)%alpha_rho(2)'   : 0.19E+00,                 \
                    'patch_icpp(2)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(2)%alpha(2)'       : 1.E+00,                   \
                    # ==========================================================
 
		    # Fluids Physical Parameters ===============================
                    'fluid_pp(1)%gamma'            : 1.E+00/(4.4E+00-1.E+00),  \
                    'fluid_pp(1)%pi_inf'           : 4.4E+00*6.E+08/(4.4E+00-1.E+00), \
                    'fluid_pp(2)%gamma'            : 1.E+00/(1.4E+00-1.E+00),  \
                    'fluid_pp(2)%pi_inf'           : 0.E+00,                   \
	            # ==========================================================

    }

# Executing MFC component
f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine)

# ==============================================================================
