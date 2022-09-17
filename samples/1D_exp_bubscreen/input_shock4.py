#!/usr/bin/env python3

import math
import json

x0      = 613.E-06
p0      = 112900.
rho0    = 960
c0      = math.sqrt( p0/rho0 )
patm    = 1.

#water props
## AKA little \gamma (see coralic 2014 eq'n (13))
n_tait  = 10.0
## AKA little \pi(see coralic 2014 eq'n (13))
B_tait  = 924.E+05 / p0

mul0    = 0.048     #viscosity
# mul0    = 1.E-12
ss      = 20.8E-3      #surface tension
# ss      = 1.E-12 ## this would turn-off surface tension
pv      = 0    #vapor pressure

# water 
# These _v and _n parameters ONLY correspond to the bubble model of Preston (2010 maybe 2008)
# (this model would replace the usual Rayleigh-plesset or Keller-miksis model (it's more compilcated))
gamma_v = 1.33
M_v     = 18.02
mu_v    = 0.8816E-05
k_v     = 0.019426

##air props
gamma_n = 1.4
M_n     = 28.97
mu_n    = 1.8E-05
k_n     = 0.02556

#air props
gamma_gas = 1.09

#reference bubble size
R0ref   = 613.E-06

pa      = 0.1 * 1.E+06 / 112900.

#Characteristic velocity
uu = math.sqrt( p0/rho0 )
#Cavitation number
Ca = (p0 - pv)/(rho0*(uu**2.))
# Ca = 1
#Weber number
We = rho0*(uu**2.)*R0ref/ss
#Inv. bubble Reynolds number
Re_inv = mul0/(rho0*uu*R0ref)

#IC setup
vf0     = 0.0024
# vf0     = 1.E-6
n0      = vf0/(math.pi*4.E+00/3.E+00)

cphysical = 981.6
t0      = x0/c0

nbubbles = 1 
myr0    = R0ref


# CFL number should be < 1  for numerical stability
# CFL = speed of sound * dt/dx
cfl     = 0.1
Nx      = 500
Ldomain = 8347.5
L       = Ldomain
dx      = L/float(Nx)
dt      = 0.0131

Lpulse  = 0.3*Ldomain
Tpulse  = Lpulse/cphysical
# Tfinal  = 100
Tfinal  = 300
#Tfinal  = 0.3*0.25*120*Tpulse*c0/x0
Nt      = int(Tfinal/dt)

Nfiles = 50.
Nout = int(math.ceil(Nt/Nfiles))
Nt = int(Nout*Nfiles)

# Configuring case dictionary
print(json.dumps({                                                                           \
    # Logistics ================================================
    'case_dir'                     : '\'.\'',                   \
    'run_time_info'                : 'T',                       \
    # ==========================================================
                                                                \
    # Computational Domain Parameters ==========================
    'x_domain%beg'                 : 0/x0,               \
    'x_domain%end'                 : 8347.5,               \
    'stretch_x'                    : 'F',                       \
    'cyl_coord'                    : 'F',                       \
    'm'                            : Nx,                        \
    'n'                            : 0,                         \
    'p'                            : 0,                         \
    'dt'                           : dt,                      \
    't_step_start'                 : 0,                         \
    't_step_stop'                  : Nt,                        \
    't_step_save'                  : Nout,   \
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
    'mapped_weno'                  : 'T',                      \
    'null_weights'                 : 'F',                      \
    'mp_weno'                      : 'T',                      \
    'riemann_solver'               : 2,                        \
    'wave_speeds'                  : 1,                        \
    'avg_state'                    : 2,                        \
    'bc_x%beg'                     : -8,                       \
    'bc_x%end'                     : -8,                       \
    # ==========================================================
                                                                \
    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,                        \
    'precision'                    : 2,                        \
    'prim_vars_wrt'                :'T',                       \
    'parallel_io'                  :'T',                       \
    'fd_order'                     : 1,                       \
    'schlieren_wrt'                :'T',                      \
    'probe_wrt'                    :'T',                   \
    'num_probes'                   : 1,                    \
    'probe(1)%x'                   : 3577.5,             \
    # ==========================================================
                                                                
    # Patch 1 _ Background =====================================
    'patch_icpp(1)%geometry'       : 1,                         \
    'patch_icpp(1)%x_centroid'     : 4173.75,                      \
    'patch_icpp(1)%length_x'       : 8347.5,                      \
    'patch_icpp(1)%vel(1)'         : 0.0,                       \
    'patch_icpp(1)%pres'           : patm,                      \
    'patch_icpp(1)%alpha_rho(1)'   : (1.-1.E-12)*(1.E+03/rho0), \
    'patch_icpp(1)%alpha(1)'       : 1.E-12,                    \
    'patch_icpp(1)%r0'             : 1.,                        \
    'patch_icpp(1)%v0'             : 0.0E+00,                   \
    # ==========================================================

    # Patch 2 Screen ===========================================
    'patch_icpp(2)%geometry'       : 1,                         \
    'patch_icpp(2)%alter_patch(1)' : 'T',                       \
    'patch_icpp(2)%x_centroid'     : 4770,                        \
    'patch_icpp(2)%length_x'       : 7155,                 \
    'patch_icpp(2)%vel(1)'         : -.14/c0,                       \
    'patch_icpp(2)%pres'           : patm,                      \
    'patch_icpp(2)%alpha_rho(1)'   : (1.-vf0)*1.E+03/rho0,   \
    'patch_icpp(2)%alpha(1)'       : vf0,                       \
    'patch_icpp(2)%r0'             : 1.,                        \
    'patch_icpp(2)%v0'             : 0.0E+00,                   \
    # ==========================================================

    # Patch 3 Shock Wave
    'patch_icpp(3)%geometry'       : 1,                         \
    'patch_icpp(3)%alter_patch(1)' : 'T',                       \
    'patch_icpp(3)%x_centroid'     : 596.25,                     \
    'patch_icpp(3)%length_x'       : 1192.5,                     \
    'patch_icpp(3)%vel(1)'         : .0283,              \
    'patch_icpp(3)%pres'           : 2.157*patm,                \
    'patch_icpp(3)%alpha_rho(1)'   : (1.-.0011)*(1.E+03/rho0), \
    'patch_icpp(3)%alpha(1)'       : .0011,                    \
    'patch_icpp(3)%r0'             : 0.77401878091102305,       \
    'patch_icpp(3)%v0'             : 0.0E+00,                   \
    'patch_icpp(3)%p0'             : 2.1577765821708557,        \
    'patch_icpp(3)%m0'             : 0.0E+00,                   \
    # ==========================================================


    # Fluids Physical Parameters ===============================
    # Surrounding liquid
    'fluid_pp(1)%gamma'             : 1.E+00/(n_tait-1.E+00),  \
    'fluid_pp(1)%pi_inf'            : n_tait*B_tait/(n_tait-1.),   \
    'fluid_pp(1)%mul0'              : mul0,     \
    'fluid_pp(1)%ss'                : ss,       \
    'fluid_pp(1)%pv'                : pv,       \
    'fluid_pp(1)%gamma_v'           : gamma_v,  \
    'fluid_pp(1)%M_v'               : M_v,      \
    'fluid_pp(1)%mu_v'              : mu_v,     \
    'fluid_pp(1)%k_v'               : k_v,      \

    # Last fluid_pp is always reserved for bubble gas state ===
    # if applicable  ==========================================
    'fluid_pp(2)%gamma'             : 1./(gamma_gas-1.),      \
    'fluid_pp(2)%pi_inf'            : 0.0E+00,      \
    'fluid_pp(2)%gamma_v'           : gamma_n,      \
    'fluid_pp(2)%M_v'               : M_n,          \
    'fluid_pp(2)%mu_v'              : mu_n,         \
    'fluid_pp(2)%k_v'               : k_n,          \
    # ==========================================================

    # Non-polytropic gas compression model AND/OR Tait EOS =====
    'pref'                  : p0,                  \
    'rhoref'                : rho0,                \
    # ==========================================================

    # Bubbles ==================================================
    'bubbles'               : 'T',                  \
    'bubble_model'          : 2,                  \
    'polytropic'            : 'F',                  \
    'polydisperse'          : 'F',                  \
    'thermal'               : 1,           \
    'R0ref'                 : myr0,                 \
    'nb'                    : 1,             \
    'Ca'                    : Ca,                   \
    'Web'                   : We,                   \
    'Re_inv'                : Re_inv,               \
    # ==========================================================
}))

# ==============================================================================
