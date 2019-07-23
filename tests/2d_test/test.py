#!/usr/bin/python

import math

x0      = 10.E-06
p0      = 101325.
rho0    = 998.
#c0      = 1475.
#p0      = rho0*c0*c0
c0      = math.sqrt( p0/rho0 )
#rho0    = p0/(c0*c0)
#patm    = 101325./p0
patm    = 1.

print 'x0', x0
print 'p0', p0
print 'c0', c0
print 'rho0', rho0

pl0     = p0 #101325.0
rhol0   = rho0 #998.2063

#water props
n_tait  = 7.15
B_tait  = 304.9E+06 / pl0
mul0    = 1.002E-03     #viscosity
ss      = 0.07275       #surface tension
pv      = 2.3388E+03    #vapor pressure

gamma_v = 1.33
M_v     = 18.02
mu_v    = 0.8816E-05
k_v     = 0.019426

#air props
gamma_n = 1.4
M_n     = 28.97
mu_n    = 1.8E-05
k_n     = 0.02556

#air props
gamma_gas = gamma_n

#reference bubble size
R0ref   = 10.E-06
#R0ref   = 1.E-03

dB      = 193
fac     = 6.666
pa      = fac*(1.E-6)*math.sqrt(10**(dB/10.)) / pl0

#Dimensionless groups
#Characteristic velocity
uu = math.sqrt( pl0/rhol0 )
#Cavitation number
Ca = (pl0 - pv)/(rhol0*(uu**2.))
#Weber number
We = rhol0*(uu**2.)*R0ref/ss
#Inv. bubble Reynolds number
Re_inv = mul0/(rhol0*uu*R0ref)

pl0 = 1.
rhol0 = 1.

#IC setup
vf0     = 1.E-5
vv      = 0.E+00 
cref    = math.sqrt( n_tait*(pl0+B_tait)/( (1-vf0)**2.E+00 ) )

nbubbles = 1 
myr0    = R0ref

#Size of the domain
x0   = 1.E-05
Rnet = 0.01*5./x0
Thickness = 0.01*1./x0
Lpulse = 0.01*0.375/x0
leng = 2.*(2.*Rnet)
monolocx = -0.75*0.5*leng

#Numerical setup
Nx      = 50
dx      = leng/float(Nx)
cfl     = 0.1

mydt    = cfl*dx/cref
Tend    = 0.4*leng/cref
Nt      = int(Tend/mydt)

print 'leng = ', leng
print 'Rnet = ', Rnet
print 'Thickness = ', Thickness

print 'maximum initial sound speed = ', cref
print 'dt, Nt = ', mydt, Nt
print 'Ma = ', vv/cref
print 'Re_inv, We, Ca = ', Re_inv, We, Ca

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
                    'x_domain%beg'                 : -leng/2.,                    \
                    'x_domain%end'                 :  leng/2.,                    \
                    'y_domain%beg'                 : -leng/2.,                   \
                    'y_domain%end'                 :  leng/2.,            \
                    'm'                            : Nx,                        \
                    'n'                            : Nx,                         \
                    'p'                            : 0,                         \
                    'dt'                           : mydt,                      \
                    't_step_start'                 : 0,                         \
                    #'t_step_stop'                  : 10,                        \
                    't_step_stop'                  : 10,                        \
                    #'t_step_save'                  : 1,    \
                    't_step_save'                  : 10,    \
		    # ==========================================================
                                                                                \
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 1,                        \
                    'model_eqns'                   : 2,                        \
                    'alt_soundspeed'               : 'F',                      \
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
                    'mp_weno'                      : 'T',                      \
                    'weno_avg'                     : 'F',                      \
                    'weno_Re_flux'                 : 'F',                      \
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
                    'bc_x%beg'                     : -2,                       \
                    'bc_x%end'                     : -2,                       \
                    'bc_y%beg'                     : -2,                      \
                    'bc_y%end'                     : -2,                      \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
                    'parallel_io'                  :'T',                       \
                    #'fd_order'                     : 1,                        \
                    #'schlieren_wrt'                :'T',                    \
                    #'schlieren_alpha(1)'           : 1.,                    \
                    #'probe_wrt'                    :'T',                    \
                    #'num_probes'                   : 3,                     \
                    #'probe(1)%x'                   : 0.,                    \
                    #'probe(1)%y'                   : 0.,                    \
                    #'probe(2)%x'                   : -Rnet,                 \
                    #'probe(2)%y'                   : 0.,                    \
                    #'probe(3)%x'                   : monolocx,              \
                    #'probe(3)%y'                   : 0.,                    \
                    #'integral_wrt'                 :'T',                    \
                    #'num_integrals'                : 3,                     \
                    #'integral(1)%xmin'             : Thickness,             \
                    #'integral(1)%xmax'             : Rnet,  \
                    #'integral(2)%xmin'             : 0,     \
                    #'integral(2)%xmax'             : 0,     \
                    #'integral(3)%xmin'             : 0,     \
                    #'integral(3)%xmax'             : 0,     \
                    # ==========================================================
                                                                               \
                    # Patch 1 _ Background =====================================
                    'patch_icpp(1)%geometry'       : 3,                     \
                    'patch_icpp(1)%x_centroid'     : 0,                   \
                    'patch_icpp(1)%y_centroid'     : 0,                   \
                    'patch_icpp(1)%length_x'       : leng,                    \
                    'patch_icpp(1)%length_y'       : leng,                    \
                    'patch_icpp(1)%vel(1)'         : 0.0,                    \
                    'patch_icpp(1)%vel(2)'         : 0.0,                    \
                    'patch_icpp(1)%pres'           : 1.0,                    \
                    'patch_icpp(1)%alpha_rho(1)'   : 1.E+00,                    \
                    'patch_icpp(1)%alpha_rho(2)'   : 1.E+00,                    \
                    'patch_icpp(1)%alpha(1)'       : 1.E-12,                \
                    'patch_icpp(1)%alpha(2)'       : 1.,                \
                    'patch_icpp(1)%r0'             : 1.,                       \
                    'patch_icpp(1)%v0'             : 0.0E+00,                       \
                    # ==========================================================

                    'fluid_pp(1)%gamma'             : 1.E+00/(n_tait-1.E+00),  \
                    'fluid_pp(1)%pi_inf'            : n_tait*B_tait/(n_tait-1.),   \

                    # Last fluid_pp is always reserved for the gaseous state ===
                    'fluid_pp(2)%gamma'             : 1./(gamma_gas-1.),      \
                    'fluid_pp(2)%pi_inf'            : 0.0E+00,      \
                    # ==========================================================

                    # SHB: Tait EOS ============================================
                    'pref'                  : 101325.,                  \
                    'rhoref'                : 998.2063,                \
                    # ==========================================================

    }

# Executing MFC component
f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine)

# ==============================================================================
