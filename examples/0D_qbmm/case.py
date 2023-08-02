#!/usr/bin/env python3

import math, json

x0      = 10.E-06
p0      = 101325.
rho0    = 1000.
c0      = math.sqrt( p0/rho0 )

# water props ==================================================================
n_tait  = 7.1
B_tait  = 306.E+06 
mul0    = -10.002E-03

#air props
gamma_gas = 1.4

#reference bubble size
R0ref   = 10.E-06

#Characteristic velocity
uu = math.sqrt( p0/rho0 )
Ca = 1.

#Inv. bubble Reynolds number
Re_inv = mul0/(rho0*uu*R0ref)

# IC setup =====================================================================
vf0     = 1.E-4
n0      = vf0/(math.pi*4.E+00/3.E+00)

cact    = 1475.
t0      = x0/c0

nbubbles = 1 
myr0    = R0ref

cfl     = 0.1
Nx      = 30
Ldomain = 20.E-03
L       = Ldomain
dx      = L/float(Nx)
dt      = cfl*dx/cact
Lpulse  = 0.3*Ldomain
Tpulse  = Lpulse/cact
Tfinal  = 0.25*10.*Tpulse*c0/x0
Nt      = int(Tfinal/dt)

dt = dt * 0.1

Nfiles = 20.
Nout   = int(math.ceil(Nt/Nfiles))
Nt     = int(Nout*Nfiles)

# ==============================================================================

# Configuring case dictionary
print(json.dumps({
    # Logistics =================================================================
    'run_time_info'                : 'T',
    # ===========================================================================
    
    # Computational Domain Parameters ===========================================
    'x_domain%beg'                 : -10.E-03,
    'x_domain%end'                 :  10.E-03,
    'stretch_x'                    : 'F',
    'cyl_coord'                    : 'F',
    'm'                            : Nx,
    'n'                            : 0,
    'p'                            : 0,
    'dt'                           : dt,
    't_step_start'                 : 0,
    't_step_stop'                  : 30000,
    't_step_save'                  : 1000,
    # ===========================================================================
    
    # Simulation Algorithm Parameters ===========================================
    'num_patches'                  : 1,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 1,
    'weno_order'                   : 3,
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'F',  
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -1,
    'bc_x%end'                     : -1,
    # ===========================================================================
    
    # Formatted Database Files Structure Parameters =============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'F',
    'fd_order'                     : 1,
    # 'schlieren_wrt'                :'T',
    'probe_wrt'                    :'T',
    'num_probes'                   : 1,
    'probe(1)%x'                   : 0.,
    # ===========================================================================
    
    'patch_icpp(1)%geometry'       : 1,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%length_x'       : 20.E-03,
    'patch_icpp(1)%vel(1)'         : 0.0,
    'patch_icpp(1)%pres'           : p0,
    'patch_icpp(1)%alpha_rho(1)'   : (1.-vf0)*rho0,
    'patch_icpp(1)%alpha(1)'       : vf0,
    'patch_icpp(1)%r0'             : 1.,
    'patch_icpp(1)%v0'             : 0.,
    # ===========================================================================
    
    # Fluids Physical Parameters ================================================
    # Surrounding liquid
    'fluid_pp(1)%gamma'            : 1.E+00/(n_tait-1.E+00),
    'fluid_pp(1)%pi_inf'           : n_tait*B_tait/(n_tait-1.),
    # Last fluid_pp is always reserved for bubble gas state ====================
    # if applicable  ===========================================================
    'fluid_pp(2)%gamma'            : 1./(gamma_gas-1.),
    'fluid_pp(2)%pi_inf'           : 0.0E+00,
    # ==========================================================================
    
    # Non-polytropic gas compression model AND/OR Tait EOS =====================
    'pref'                         : p0,
    'rhoref'                       : rho0,
    # ==========================================================================
    
    # Bubbles ==================================================================
    'bubbles'                      : 'T',
    'bubble_model'                 : 3,
    'polytropic'                   : 'T',
    'thermal'                      : 3,
    'R0ref'                        : myr0,
    'nb'                           : 1,
    'Ca'                           : Ca,
    'Re_inv'                       : Re_inv,
    'qbmm'                         : 'T',
    'dist_type'                    : 1,
    'sigR'                         : 0.1,
    'sigV'                         : 0.1,
    'rhoRV'                        : 0.0,
    
}))

