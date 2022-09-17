#!/usr/bin/env python2

import math
import json

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
# gamma_gas = gamma_n
gamma_gas = 1.4

#reference bubble size
R0ref   = 10.E-06

pa      = 0.1 * 1.E+06 / 101325.

#Characteristic velocity
uu = math.sqrt( p0/rho0 )
#Cavitation number
# Ca = (p0 - pv)/(rho0*(uu**2.))
Ca = 1.
#Weber number
# We = rho0*(uu**2.)*R0ref/ss
We = p0*R0ref/ss
#Inv. bubble Reynolds number
Re_inv = mul0/(rho0*uu*R0ref)

#IC setup
vf0     = 4.E-5
n0      = vf0/(math.pi*4.E+00/3.E+00)

cact    = 1475.
t0      = x0/c0

nbubbles = 1 
myr0    = R0ref

cfl     = 0.1
Nx      = 400
Ldomain = 20.E-03
L       = Ldomain/x0
dx      = L/float(Nx)
dt      = cfl*dx*c0/cact
Lpulse  = 0.3*Ldomain
Tpulse  = Lpulse/cact
Tfinal  = 0.25*10.*Tpulse*c0/x0
Nt      = int(Tfinal/dt)

dt = dt * 0.1
# print('dt: ',dt)

Nfiles = 20.
Nout = int(math.ceil(Nt/Nfiles))
Nt = int(Nout*Nfiles)

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================
    'case_dir'                     : '\'.\'',
    'run_time_info'                : 'F',
    # ==========================================================

    # Computational Domain Parameters ==========================
    'x_domain%beg'                 : -10.E-03/x0,
    'x_domain%end'                 :  10.E-03/x0,
    'stretch_x'                    : 'F',
    'cyl_coord'                    : 'F',
    'm'                            : Nx,
    'n'                            : 0,
    'p'                            : 0,
    'dt'                           : 0.002,
    't_step_start'                 : 0,
    't_step_stop'                  : 8000,
    # 't_step_stop'                  : 4,
    't_step_save'                  : 8000,
    # 't_step_save'                  : 1,
    # ==========================================================

    # Simulation Algorithm Parameters ==========================
    'num_patches'                  : 2,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 3,
    'weno_vars'                    : 2,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'T',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -8,
    'bc_x%end'                     : -8,
    # ==========================================================

    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'F',
    'fd_order'                     : 1,
    #'schlieren_wrt'                :'T',
    'probe_wrt'                    :'T',
    'num_probes'                   : 1,
    'probe(1)%x'                   : 0.,
    # ==========================================================
                                                                
    # Patch 1 _ Background =====================================
    'patch_icpp(1)%geometry'       : 1,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%length_x'       : 20.E-03/x0,
    'patch_icpp(1)%vel(1)'         : 0.0,
    'patch_icpp(1)%pres'           : 1,
    'patch_icpp(1)%alpha_rho(1)'   : (1.-1.E-12)*1.E+03/rho0,
    'patch_icpp(1)%alpha(1)'       : 1.E-12,
    'patch_icpp(1)%r0'             : 1.,
    'patch_icpp(1)%v0'             : 0.0E+00,
    # ==========================================================

    # Patch 2 Screen ===========================================
    'patch_icpp(2)%geometry'       : 1,
    'patch_icpp(2)%x_centroid'     : 0.,
    'patch_icpp(2)%length_x'       : 5.E-03/x0,
    # 'patch_icpp(1)%length_x'       : 20.E-03/x0,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%vel(1)'         : 0.0,
    'patch_icpp(2)%pres'           : 1,
    'patch_icpp(2)%alpha_rho(1)'   : (1.-vf0)*1.E+03/rho0,
    'patch_icpp(2)%alpha(1)'       : vf0,
    'patch_icpp(2)%r0'             : 1.,
    'patch_icpp(2)%v0'             : 0.,
    # ==========================================================

    # Fluids Physical Parameters ===============================
    # Surrounding liquid
    'fluid_pp(1)%gamma'             : 1.E+00/(n_tait-1.E+00),
    'fluid_pp(1)%pi_inf'            : n_tait*B_tait/(n_tait-1.),
    'fluid_pp(1)%mul0'              : mul0,
    'fluid_pp(1)%ss'                : ss,
    'fluid_pp(1)%pv'                : pv,
    'fluid_pp(1)%gamma_v'           : gamma_v,
    'fluid_pp(1)%M_v'               : M_v,
    'fluid_pp(1)%mu_v'              : mu_v,
    'fluid_pp(1)%k_v'               : k_v,

    # Last fluid_pp is always reserved for bubble gas state ===
    # if applicable  ==========================================
    'fluid_pp(2)%gamma'             : 1./(gamma_gas-1.),
    'fluid_pp(2)%pi_inf'            : 0.0E+00,
    'fluid_pp(2)%gamma_v'           : gamma_n,
    'fluid_pp(2)%M_v'               : M_n,
    'fluid_pp(2)%mu_v'              : mu_n,
    'fluid_pp(2)%k_v'               : k_n,
    # ==========================================================

    # Non-polytropic gas compression model AND/OR Tait EOS =====
    'pref'                  : p0,
    'rhoref'                : rho0,
    # ==========================================================

    # Bubbles ==================================================
    'bubbles'               : 'T',
    'bubble_model'          : 2,
    'polytropic'            : 'T',
    'polydisperse'          : 'T',
    'R0_type'               : 1,
    'poly_sigma'            : 0.3,
    'thermal'               : 3,
    'R0ref'                 : myr0,
    # 'nb'                    : 3,
    'nb'                    : 3,
    'Ca'                    : Ca,
    # 'Web'                   : We,
    # 'Re_inv'                : Re_inv,
    'qbmm'               : 'T',
    'dist_type'          : 2,
    'sigR'               : 0.1,
    'sigV'               : 0.1,
    'rhoRV'              : 0.0,
    # ==========================================================

    # Acoustic source ==========================================
    'Monopole'                  : 'T',
    'num_mono'                  : 1,
    'Mono(1)%loc(1)'            : -5.E-03/x0,
    'Mono(1)%npulse'            : 1,
    'Mono(1)%dir'               : 1.,
    'Mono(1)%pulse'             : 1,
    'Mono(1)%mag'               : 1*pa,
    'Mono(1)%length'            : (1./(300000.))*cact/x0,
    # ==========================================================
}))

# ==============================================================================
