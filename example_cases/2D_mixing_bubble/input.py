#!/usr/bin/env python2
import math

x0      = 10.E-06
p0      = 101325.
rho0    = 1.E+03
c0      = math.sqrt( p0/rho0 )
patm    = 1.

#water props
n_tait  = 7.1
B_tait  = 306.E+06 / p0
mul0    = 1.002E-03     #viscosity
ss      = 0.07275       #surface tension
pv      = 2.3388E+03    #vapor pressure

gamma_gas = 1.4
R0ref   = 10.E-06
pa      = 0.1 * 1.E+06 / 101325.

#Characteristic velocity
uu = math.sqrt( p0/rho0 )

#Cavitation number
Ca = (p0 - pv)/(rho0*(uu**2.))
#Ca = 1.

#Weber number
# We = rho0*(uu**2.)*R0ref/ss
We = p0*R0ref/ss

#Inv. bubble Reynolds number
Re_inv = mul0/(rho0*uu*R0ref)


Nx      = 49
Ny      = Nx

dx      = 1./(Nx+1.)
dy      = 1./(Ny+1.)

Nt      = 2550 #500*2

cfl     = 0.03E+00
mydt    = cfl*dx/uu
Ttot    = mydt*Nt

myr0    = 1.E+00
vv      = 0.1*uu

vf0 = 1.E-5

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
#engine = 'parallel'

# ==============================================================================

# Case Analysis Configuration ==================================================

# Selecting MFC component
comp_name = argv[1].strip()

# Configuring case dictionary
case_dict =                                                                     \
    {                                                                           \
                    # Logistics ================================================
                    'case_dir'                     : '\'.\'',                   \
                    'run_time_info'                : 'T',                       \
                    'nodes'                        : 1,                         \
                    'ppn'                          : 1,                         \
                    'queue'                        : 'normal',                  \
                    'walltime'                     : '24:00:00',                \
                    'mail_list'                    : '',                        \
                    # ==========================================================
                                                                                \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 : 0.E+00,                    \
                    'x_domain%end'                 : 1.E+00,                    \
                    'y_domain%beg'                 : 0.E+00,                   \
                    'y_domain%end'                 : 1.E+00,            \
                    'cyl_coord'                    : 'F',                       \
                    'm'                            : Nx,                        \
                    'n'                            : Ny,                         \
                    'p'                            : 0,                         \
                    'dt'                           : mydt,                      \
                    't_step_start'                 : 0,                         \
                    't_step_stop'                  : Nt,                        \
                    't_step_save'                  : int(math.ceil(Nt/10.)),         \
		    # ==========================================================
                                                                                \
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 3,                        \
                    'model_eqns'                   : 2,                        \
                    'alt_soundspeed'               : 'F',                      \
                    'num_fluids'                   : 1,                        \
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
                    'mp_weno'                      : 'T',                      \
                    'weno_avg'                     : 'F',                      \
                    'weno_Re_flux'                 : 'F',                      \
		    'riemann_solver'               : 2,                        \
                    'wave_speeds'                  : 1,                        \
                    'avg_state'                    : 2,                        \
                    'commute_err'                  : 'F',                      \
                    'split_err'                    : 'F',                      \
                    'reg_eps'                      : 1.E+00,                   \
                    'alt_crv'                      : 'F',                      \
                    'bc_x%beg'                     : -1,                       \
                    'bc_x%end'                     : -1,                       \
                    'bc_y%beg'                     : -5,                      \
                    'bc_y%end'                     : -5,                      \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
		    'parallel_io'                  :'F',                       \
		    # ==========================================================
                                                                                
		    # Patch 1 ==================================================
                    'patch_icpp(1)%geometry'       : 3,                     \
                    'patch_icpp(1)%x_centroid'     : 0.50,                   \
                    'patch_icpp(1)%y_centroid'     : 0.25,                   \
                    'patch_icpp(1)%length_x'       : 1.00,                    \
                    'patch_icpp(1)%length_y'       : 0.50,                    \
                    'patch_icpp(1)%alpha_rho(1)'   : (1.-1.E-12)*1.E+00,                    \
                    'patch_icpp(1)%vel(1)'         : vv,                    \
                    'patch_icpp(1)%vel(2)'         : 0.00,                    \
                    'patch_icpp(1)%pres'           : 1.00,                    \
                    'patch_icpp(1)%alpha(1)'       : 1.E-12,                \
                    'patch_icpp(1)%r0'             : 1.E+00,                       \
                    'patch_icpp(1)%v0'             : 0.0E+00,                       \
                    # ==========================================================

                    # Patch 2 ==================================================
                    'patch_icpp(2)%geometry'       : 3,                     \
                    'patch_icpp(2)%x_centroid'     : 0.50,                   \
                    'patch_icpp(2)%y_centroid'     : 0.75,                   \
                    'patch_icpp(2)%length_x'       : 1.00,                    \
                    'patch_icpp(2)%length_y'       : 0.50,                    \
                    'patch_icpp(2)%alpha_rho(1)'   : (1-1.E-12)*1.E+00,                    \
                    'patch_icpp(2)%vel(1)'         : -1.*vv,                    \
                    'patch_icpp(2)%vel(2)'         : 0.0,                    \
                    'patch_icpp(2)%pres'           : 1.0,                    \
                    'patch_icpp(2)%alpha(1)'       : 1.E-12,                \
                    'patch_icpp(2)%r0'             : 1.E+00,                       \
                    'patch_icpp(2)%v0'             : 0.0E+00,                       \
                    # ==========================================================

                    # ==========================================================
                    'patch_icpp(3)%geometry'       : 2,                     \
                    'patch_icpp(3)%x_centroid'     : 5.E-01,                   \
                    'patch_icpp(3)%y_centroid'     : 0.5,                   \
                    'patch_icpp(3)%radius'          : 0.1,                    \
                    'patch_icpp(3)%alter_patch(1)' : 'T',                       \
                    'patch_icpp(3)%alter_patch(2)' : 'T',                       \
                    'patch_icpp(3)%alpha_rho(1)'   : (1-vf0)*1.E+00,                    \
                    'patch_icpp(3)%vel(1)'         : 0.00,                    \
                    'patch_icpp(3)%vel(2)'         : 0.00,                    \
                    'patch_icpp(3)%pres'           : 1.00,                    \
                    'patch_icpp(3)%alpha(1)'       : vf0,                \
                    'patch_icpp(3)%r0'             : 1.E+00,                       \
                    'patch_icpp(3)%v0'             : 0.0E+00,                       \
                    # ==========================================================

                    # SHB: Bubbles =============================================
                    # 'perturb_flow'              : 'T',                  \
                    # 'perturb_flow_fluid'             : 1,                 \

                    # Fluids Physical Parameters ===============================
                    # Surrounding liquid
                    'fluid_pp(1)%gamma'             : 1.E+00/(n_tait-1.E+00),  \
                    'fluid_pp(1)%pi_inf'            : n_tait*B_tait/(n_tait-1.),   \

                    # Last fluid_pp is always reserved for bubble gas state ===
                    # if applicable  ==========================================
                    'fluid_pp(2)%gamma'             : 1./(gamma_gas-1.),      \
                    'fluid_pp(2)%pi_inf'            : 0.0E+00,      \
                    # ==========================================================
                    # SHB: Tait EOS ============================================
                    'pref'                  : p0,                  \
                    'rhoref'                : rho0,                \
	            # ==========================================================

                    # Bubbles ==================================================
                    'bubbles'               : 'T',                  \
                    'bubble_model'          : 2,                    \
                    'polytropic'            : 'T',                  \
                    # 'polydisperse'          : 'T',                  \
                    'R0_type'               : 1,                    \
                    #'polydisperse'          : 'F',                  \
                    # 'poly_sigma'            : 0.1,                  \
                    'thermal'               : 3,                    \
                    'R0ref'                 : myr0,                 \
                    # 'nb'                    : 3,                    \
                    'nb'                    : 1,                    \
                    'Ca'                    : Ca,                   \
                    # 'Web'                   : We,                   \
                    # 'Re_inv'                : Re_inv,               \
                    #'qbmm'               : 'T',                     \
                    #'nnode'              : 4,                       \
                    #'dist_type'          : 2,                       \
                    #'sigR'               : 0.1,                     \
                    #'sigV'               : 0.1,                     \
                    #'rhoRV'              : 0.0,                     \
                    # ==========================================================

    }

# Executing MFC component
f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine)

# ==============================================================================
