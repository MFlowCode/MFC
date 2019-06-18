#!/usr/bin/python


# MFC v3.0 - Cases - Kapila model - Shock tube - 1D
# Description: Shock tube.
# Author: Kevin Schmidmayer
# Date: 02/26/18


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
                    'run_time_info'                : 'F',                      \
                    'nodes'                        : 1,                        \
                    'ppn'                          : 1,                        \
                    'queue'                        : 'normal',                 \
                    'walltime'                     : '24:00:00',               \
                    'mail_list'                    : '',                       \
                    # ==========================================================
                                                                               \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 : 0.E+00,                   \
                    'x_domain%end'                 : 1.E+00,                   \
                    'm'                            : 299,                      \
                    'n'                            : 0,                        \
                    'p'                            : 0,                        \
                    'dt'                           : 4.E-08,                   \
                    't_step_start'                 : 0,                        \
                    't_step_stop'                  : 500,                     \
                    't_step_save'                  : 500,                     \
		    # ==========================================================
                                                                               \
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 2,                        \
                    'model_eqns'                   : 2,                        \
                    'alt_soundspeed'               : 'T',                      \
                    'num_fluids'                   : 2,                        \
		    'adv_alphan'                   : 'T',                      \
		    'mpp_lim'                      : 'F',                      \
		    'mixture_err'                  : 'F',                      \
		    'time_stepper'                 : 3,                        \
                    'weno_vars'                    : 2,                        \
                    'weno_order'                   : 5,                        \
                    'weno_eps'                     : 1.E-16,                   \
                    'char_decomp'                  : 'F',                      \
                    'mapped_weno'                  : 'T',                      \
                    'null_weights'                 : 'F',                      \
                    'mp_weno'                      : 'F',                      \
		    'riemann_solver'               : 2,                        \
                    'wave_speeds'                  : 1,                        \
                    'avg_state'                    : 2,                        \
                    'commute_err'                  : 'F',                      \
                    'split_err'                    : 'F',                      \
                    'bc_x%beg'                     : -3,                       \
                    'bc_x%end'                     : -3,                       \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        \
                    'precision'                    : 1,                        \
                    'prim_vars_wrt'                :'T',                       \
		    'parallel_io'                  :'F',                       \
                    # ==========================================================
                                                                                
		    # Patch 1: High pressured water ============================
                    'patch_icpp(1)%geometry'       : 1,                        \
                    'patch_icpp(1)%x_centroid'     : 0.5E+00,                  \
                    'patch_icpp(1)%length_x'       : 1.E+00,                   \
                    'patch_icpp(1)%vel(1)'         : 0.E+00,                   \
                    'patch_icpp(1)%pres'           : 1.E+09,                   \
                    'patch_icpp(1)%alpha_rho(1)'   : 1000.E+00,                \
                    'patch_icpp(1)%alpha_rho(2)'   : 0.,                   \
                    'patch_icpp(1)%alpha(1)'       : 1.E+00,                   \
                    'patch_icpp(1)%alpha(2)'       : 0.E+00,                   \
                    # ==========================================================

                    # Patch 2: Air bubble ======================================
                    'patch_icpp(2)%geometry'       : 1,                        \
                    'patch_icpp(2)%x_centroid'     : 0.85E+00,                 \
                    'patch_icpp(2)%length_x'       : 0.3E+00,                  \
                    'patch_icpp(2)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(2)%vel(1)'         : 0.E+00,                   \
                    'patch_icpp(2)%pres'           : 1.E+05,                   \
                    'patch_icpp(2)%alpha_rho(1)'   : 0.,                   \
                    'patch_icpp(2)%alpha_rho(2)'   : 50.E+0,                  \
                    'patch_icpp(2)%alpha(1)'       : 0,                   \
                    'patch_icpp(2)%alpha(2)'       : 1.,                   \
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
