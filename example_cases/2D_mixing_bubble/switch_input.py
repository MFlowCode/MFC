#!/usr/bin/env python2
import math

# Command to navigate between directories
from os import chdir

# Command to acquire directory path
from os.path import dirname

# Command to acquire script name and module search path
from sys import argv, path, exit

# Command to execute the MFC components
# from m_python_proxy import f_execute_mfc_component

# Navigating to script directory
if len(dirname(argv[0])) != 0: chdir(dirname(argv[0]))

# Adding master_scripts directory to module search path
mfc_dir = '../../src'; path[:0] = [mfc_dir + '/master_scripts']

from m_python_proxy import f_execute_mfc_component_SHB

# Selecting MFC component
comp_name = argv[1].strip()
# Select type of simulation
restart_name = argv[2].strip()

# x0      = 10.E-06
x0      = 1.
p0      = 101325.
rho0    = 1000.
u0      = math.sqrt( p0/rho0 )
c0      = 1475.

n_tait  = 7.1
B_tait  = 306.E+06 / p0
gamma_gas = 1.4

# Velocity
uu = 4. / u0

#Cavitation number
Ca = 1.

Ly =  0.5/x0
Lx = Ly*2/x0

Ny      = 79
Nx      = Ny*2+1
dx      = Lx/float(Nx)
dy      = Ly/float(Ny)

# Time stepping parameters
cfl = 0.3
dt  = cfl*dx/(c0/u0)
T   = 15.
Ntfinal   = int(T/dt)
Ntrestart = int(Ntfinal/5.)

# Init
t_start = 0
Nfiles  = 5E1
t_save  = int(math.ceil(Ntrestart/float(Nfiles)))
Nt      = t_save*Nfiles
Ntrestart = Nt
bc_y    = 8

if restart_name == 'run':
    # Simulate
    t_start = Ntrestart
    Nfiles  = 5E2
    t_save  = int(math.ceil((Ntfinal-t_start)/float(Nfiles)))
    Nt      = t_save*Nfiles
    bc_y    = 5
elif restart_name != 'init':
    sys.exit("incorrect restart parameter")

ang = 1.

myr0 = 1.E+00
vf0 = 1.E-12

# ==============================================================================


# Case Analysis Configuration ==================================================
sub_name = restart_name + '+' + 'dimless_mixing_layer'

nnode = 1
np = 24

engine = 'parallel'
if (comp_name=='pre_process'): engine = 'serial'

# Configuring case dictionary
case_dict =                                                                     \
    {                                                                           \
                    # Logistics ================================================
                    'case_dir'                     : '\'.\'',                   \
                    'run_time_info'                : 'T',                       \
                    'nodes'                        : nnode,                     \
                    'ppn'                          : np,                        \
                    'queue'                        : 'normal',                  \
                    'walltime'                     : '24:00:00',                \
                    'mail_list'                    : '',                        \
                    # ==========================================================
                                                                                \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 : -Lx/2.,                    \
                    'x_domain%end'                 :  Lx/2.,                    \
                    'y_domain%beg'                 : -Ly/2.,                    \
                    'y_domain%end'                 :  Ly/2.,                    \
                    'cyl_coord'                    : 'F',                       \
                    'm'                            : Nx,                        \
                    'n'                            : Ny,                        \
                    'p'                            : 0,                         \
                    'dt'                           : dt,                      \
                    't_step_start'                 : t_start,                   \
                    't_step_stop'                  : Nt,                        \
                    't_step_save'                  : t_save ,   \
                    # ==========================================================
                                                                                \
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 2,                        \
                    'model_eqns'                   : 2,                        \
                    'alt_soundspeed'               : 'F',                      \
                    'num_fluids'                   : 3,                        \
                    'adv_alphan'                   : 'T',                      \
                    'mpp_lim'                      : 'T',                      \
                    'mixture_err'                  : 'T',                      \
                    'time_stepper'                 : 3,                        \
                    'weno_vars'                    : 2,                        \
                    'weno_order'                   : 3,                        \
                    'weno_eps'                     : 1.E-16,                   \
                    'char_decomp'                  : 'F',                      \
                    'mapped_weno'                  : 'T',                      \
                    'null_weights'                 : 'F',                      \
                    # 'mp_weno'                      : 'T',                      \
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
                    'bc_y%beg'                     : -bc_y,                       \
                    'bc_y%end'                     : -bc_y,                       \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
                    'parallel_io'                  :'T',                       \
                    'probe_wrt'                    :'T',                        \
                    'num_probes'                   : 1,                         \
                    'probe(1)%x'                   : 0.,                        \
                    'probe(1)%y'                   : 0.,                        \
                    'fd_order'                     : 1,                         \
                    # 'schlieren_wrt'                :'F',                        \
                    # ==========================================================
                                                                                
                    # Patch 1 ==================================================
                    'patch_icpp(1)%geometry'       : 3,                         \
                    'patch_icpp(1)%x_centroid'     : 0.,                        \
                    'patch_icpp(1)%y_centroid'     : 0.,                        \
                    'patch_icpp(1)%length_x'       : Lx,                        \
                    'patch_icpp(1)%length_y'       : Ly,                        \
                    'patch_icpp(1)%alpha_rho(1)'   : (1.-vf0)*1.,               \
                    'patch_icpp(1)%alpha_rho(2)'   : (1.-vf0)*1.E-12,           \
                    'patch_icpp(1)%alpha_rho(3)'   : vf0*1.E-3,                 \
                    'patch_icpp(1)%alpha(1)'       : (1-vf0)*1.,                \
                    'patch_icpp(1)%alpha(2)'       : (1-vf0)*1.E-12,            \
                    'patch_icpp(1)%alpha(3)'       : vf0,                       \
                    'patch_icpp(1)%vel(1)'         : uu,                        \
                    'patch_icpp(1)%vel(2)'         : 0.00,                      \
                    'patch_icpp(1)%pres'           : 1.00,                      \
                    'patch_icpp(1)%r0'             : 1.E+00,                    \
                    'patch_icpp(1)%v0'             : 0.0E+00,                   \
                    # ==========================================================

                    # Patch 2 ==================================================
                    'patch_icpp(2)%geometry'       : 4,                     \
                    'patch_icpp(2)%alter_patch(1)' : 'T',                     \
                    'patch_icpp(2)%x_centroid'     : 0.,                   \
                    'patch_icpp(2)%y_centroid'     : 0.,                   \
                    'patch_icpp(2)%normal(1)'      :  math.sin(math.pi*ang/180.),     \
                    'patch_icpp(2)%normal(2)'      : -math.cos(math.pi*ang/180.),     \
                    'patch_icpp(2)%alpha_rho(1)'   : (1.-vf0)*1.E-12,           \
                    'patch_icpp(2)%alpha_rho(2)'   : (1.-vf0)*1.,               \
                    'patch_icpp(2)%alpha_rho(3)'   : vf0*1.E-3,               \
                    'patch_icpp(2)%alpha(1)'       : (1-vf0)*1.E-12,            \
                    'patch_icpp(2)%alpha(2)'       : (1-vf0)*1.,                \
                    'patch_icpp(2)%alpha(3)'       : vf0,                \
                    'patch_icpp(2)%vel(1)'         : -1.*uu,                    \
                    'patch_icpp(2)%vel(2)'         : 0.0,                    \
                    'patch_icpp(2)%pres'           : 1.0,                    \
                    'patch_icpp(2)%r0'             : 1.E+00,                       \
                    'patch_icpp(2)%v0'             : 0.0E+00,                       \
                    # 'patch_icpp(2)%normal(1)'      : 0.00624987793326E+00,     \
                    # 'patch_icpp(2)%normal(2)'      :-0.99998046932219E+00,     \
                    # 'patch_icpp(2)%length_x'       : Lx,                   \
                    # 'patch_icpp(2)%length_y'       : Ly/2.,                    \
                    # ==========================================================

                    # # ==========================================================
                    # 'patch_icpp(3)%geometry'       : 2,                     \
                    # 'patch_icpp(3)%x_centroid'     : 5.E-01,                   \
                    # 'patch_icpp(3)%y_centroid'     : 0.5,                   \
                    # 'patch_icpp(3)%radius'          : 0.1,                    \
                    # 'patch_icpp(3)%alter_patch(1)' : 'T',                       \
                    # 'patch_icpp(3)%alter_patch(2)' : 'T',                       \
                    # 'patch_icpp(3)%alpha_rho(1)'   : (1-vf0)*1.E+00,                    \
                    # 'patch_icpp(3)%vel(1)'         : 0.00,                    \
                    # 'patch_icpp(3)%vel(2)'         : 0.00,                    \
                    # 'patch_icpp(3)%pres'           : 1.00,                    \
                    # 'patch_icpp(3)%alpha(1)'       : vf0,                \
                    # 'patch_icpp(3)%r0'             : 1.E+00,                       \
                    # 'patch_icpp(3)%v0'             : 0.0E+00,                       \
                    # # ==========================================================

                    # SHB: Bubbles =============================================
                    # 'perturb_flow'              : 'T',                      \
                    # 'perturb_flow_fluid'        : 1,                        \
                    # ==========================================================

                    # Fluids Physical Parameters ===============================
                    # Surrounding liquid
                    'fluid_pp(1)%gamma'             : 1.E+00/(n_tait-1.E+00),  \
                    'fluid_pp(1)%pi_inf'            : n_tait*B_tait/(n_tait-1.),   \
                    'fluid_pp(2)%gamma'             : 1.E+00/(n_tait-1.E+00),  \
                    'fluid_pp(2)%pi_inf'            : n_tait*B_tait/(n_tait-1.),   \
                    'fluid_pp(3)%gamma'             : 1./(gamma_gas-1.),      \
                    'fluid_pp(3)%pi_inf'            : 0.0E+00,      \
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
# f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine)
f_execute_mfc_component_SHB(comp_name, case_dict, mfc_dir, engine, sub_name)

# ==============================================================================
