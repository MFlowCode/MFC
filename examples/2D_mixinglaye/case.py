#!/usr/bin/env python3

import math
import json


# ORIGINAL PROPERTIES =========================================================
# Water
# nu0     = 1.E-06            # [m2/s] kinematic viscosity of water
rho0    = 1.E+03            # [kg/m3] density of water
gamma0  = 7.1               # [1] specific heat ratio
pi_inf0 = 3.06E+08          # [N/m2] water stiffness
mul0    = 1.E-03            # [kg/ms] dynamic viscosity of water 
ss      = 72.8E-03          # [N/m] surface tension 
pv      = 2.3388E+03        # [N/m2] vapor pressure of water at 20degC

# Vapor
# gamma_v = 1.33              # [1] specific heat ratio
# M_v     = 18.02             # [?] molecular weight
# mu_v    = 0.8816E-05        # [kg/ms] dynamic viscosity
# k_v     = 0.019426          # [?] ?

# Air
gamma_n = 1.4               # [1] specific heat ratio
# M_n     = 28.97             # [?] molecular weight
# mu_n    = 1.8E-05           # [kg/ms] dynamic viscosity
# k_n     = 0.02556           # [?] ?
# rho_n   = 1.2041            # [kg/m3] density

# Reynolds number
Re0     = 50.               # [1] Reynolds number based on the upper stream
                            # velocity and the initial vorticity thickness
x0      = 0.002475          # [m] initial vorticity thickness
u0      = 3.4343            # [m/s] upper stream velocity

# Cavitation number
Ca      = 1.0               # [1] Cavitation number = (p0-pv)/(0.5*rho0*u0^2)
pres0   = pv + Ca*(0.5*rho0*u0**2)   # [N/m2] pressure of water

# Sound speed
c0      = math.sqrt(gamma0*(pres0+pi_inf0)/rho0)

# Bubble dynamics
R0ref   = x0                  # [m] Initial bubble radius
vf0     = 10.E-04                   # [1] Initial void fraction
uref    = math.sqrt(pres0/rho0)     # [m/s] Reference velocity based on initial 
                                    #       bubble radius and approx bubble 
                                    #       natural frequency
Ca_bub  = (pres0-pv)/pres0          # [1] Cavitation number for bubble
We_bub  = rho0*uref**2*R0ref/ss     # [1] Weber number for bubble
Re_bub  = rho0*uref*R0ref/mul0      # [1] Reynolds number for bubble
uratio  = uref/u0                   # [1] 
rratio  = R0ref/x0                  # [1] 

# MODIFIED PROPERTIES FOR ARTIFICIAL MACH NUMBER ==============================
Ma_t    = 0.002               # [1] target Mach number
pi_fac  = (rho0*u0**2/(gamma0*Ma_t**2)-pres0)/pi_inf0
pi_inf0 = pi_inf0*pi_fac
c0      = math.sqrt(gamma0*(pres0+pi_inf0)/rho0)
# =============================================================================

# SIMULATION SETUP ============================================================
# Domain size
Lx = 59. / R0ref
Ly = 59. / R0ref
# Lz = 59.0

# Number of grid cells
Nx = 249
Ny = 249
# Nz = 191

# Grid spacing
dx      = Lx/float(Nx+1)
dy      = Ly/float(Ny+1)
# dz      = Lz/float(Nz+1)

# Time advancement
cfl     = 0.01
T       = 50.
#dt      = cfl*dx/(1.+c0/u0)
dt      = cfl * dx * uref / (c0 + u0)
Ntfinal = int(T/dt)
Ntstart = 0
Nfiles  = 50
t_save  = int(math.ceil((Ntfinal-0)/float(Nfiles)))
Nt      = t_save*Nfiles
t_step_start    = Ntstart
t_step_stop     = int(Nt)

# ==============================================================================


# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                 : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                  :  0.,
    'x_domain%end'                  :  Lx,
    'y_domain%beg'                  : -Ly/2.,
    'y_domain%end'                  :  Ly/2.,
    'm'                             : Nx,
    'n'                             : Ny,
    'p'                             : 0, #Nz,
    'dt'                            : dt,
    't_step_start'                  : Ntstart,
    't_step_stop'                   : t_step_stop,
    't_step_save'                   : t_save,
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                   : 1,
    'model_eqns'                    : 2,
    'num_fluids'                    : 1,
    'time_stepper'                  : 3,
    'weno_order'                    : 5,
    'weno_eps'                      : 1.E-16,
    'weno_Re_flux'                  : 'F',
    'mapped_weno'                   : 'T',
    'riemann_solver'                : 2,
    'wave_speeds'                   : 1,
    'avg_state'                     : 2,
    'bc_x%beg'                      : -1,
    'bc_x%end'                      : -1,
    'bc_y%beg'                      : -5,
    'bc_y%end'                      : -5,
    # 'bc_z%beg'                      : -1,
    # 'bc_z%end'                      : -1,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                        : 1,
    'precision'                     : 2,
    'cons_vars_wrt'                 :'T',
    'prim_vars_wrt'                 :'T',
    'parallel_io'                   :'T',
    'fd_order'                      : 1,
    # 'omega_wrt(1)'                  :'T',
    # 'omega_wrt(2)'                  :'T',
    'omega_wrt(3)'                  :'T',
    # 'qm_wrt'                        :'T',
    # ==========================================================================

    # Patch 1 ==================================================================
    'patch_icpp(1)%geometry'        : 3,
    'patch_icpp(1)%x_centroid'      : Lx/2.,
    'patch_icpp(1)%y_centroid'      : 0.,
    # 'patch_icpp(1)%z_centroid'      : Lz/2.,
    'patch_icpp(1)%length_x'        : Lx ,
    'patch_icpp(1)%length_y'        : Ly ,
    # 'patch_icpp(1)%length_z'        : Lz,
    'patch_icpp(1)%alpha_rho(1)'    : (1.0-vf0) * rho0 / rho0,
    'patch_icpp(1)%alpha(1)'        : vf0,
    'patch_icpp(1)%vel(1)'          : u0 / uref,
    'patch_icpp(1)%vel(2)'          : 0.,
    # 'patch_icpp(1)%vel(3)'          : 0.,
    'patch_icpp(1)%pres'            : pres0 / pres0,
    'patch_icpp(1)%r0'              : 1.,
    'patch_icpp(1)%v0'              : 0.,
    # ==========================================================================

    # Bubble ===================================================================
    'bubbles'                       : 'T',
    'bubble_model'                  : 3,
    'polytropic'                    : 'T',
    'polydisperse'                  : 'F',
    'thermal'                       : 1,
    'nb'                            : 1,
    'Ca'                            : Ca_bub,
    'Web'                           : We_bub,
    'Re_inv'                        : 1./Re_bub,
    'R0ref'                         : R0ref,
    'R0_type'					    : 1,
    #'uratio'                        : uratio,
    #'rratio'                        : rratio,
    # 'R0ref'                         : R0ref,  # For polydisperse = T
    # 'R0_type'                       : 1,      # For polydisperse = T
    # 'poly_sigma'                    : 0.3,    # For polydisperse = T
    # = ========================================================================
    
    # Mixing layer === =========================================================
    'vel_profile'                   : 'T',
    'instability_wave'              : 'T',
    # ==========================================================================
    
    #'pi_fac'                        : pi_fac,

    'pref'                  : pres0,
    'rhoref'                : rho0,

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(gamma0-1.E+00),
    'fluid_pp(1)%pi_inf'           : gamma0*(pi_inf0/(pres0))/(gamma0-1.E+00),
    'fluid_pp(1)%Re(1)'            : Re0,
   # 'fluid_pp(1)%mul0'             : mul0,
   # 'fluid_pp(1)%ss'               : ss,
   # 'fluid_pp(1)%pv'               : pv,
   # 'fluid_pp(1)%gamma_v'          : gamma_v,
   # 'fluid_pp(1)%M_v'              : M_v,
   # 'fluid_pp(1)%mu_v'             : mu_v,
   # 'fluid_pp(1)%k_v'              : k_v,

    'fluid_pp(2)%gamma'            : 1.E+00/(gamma_n-1.E+00),
    'fluid_pp(2)%pi_inf'           : 0.E+00,
   # 'fluid_pp(2)%gamma_v'          : gamma_n,
   # 'fluid_pp(2)%M_v'              : M_n,
   # 'fluid_pp(2)%mu_v'             : mu_n,
   # 'fluid_pp(2)%k_v'              : k_n,
    # =========================================================================
}))

# ==============================================================================
