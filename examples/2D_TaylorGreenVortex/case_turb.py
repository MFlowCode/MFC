#!/usr/bin/env python3

import math
import json

gamma = 1.4
Re_lit = 1600
Ma = 0.1

# Dimensions of the Vortex
l = 1.0
L = math.pi*l

pres_i = 1.E+05
alpha_rho1_i = 1.22

# Turbulent Mach Number
c = (gamma * pres_i/alpha_rho1_i) ** 0.5
vel1_i = Ma * c

Mu = alpha_rho1_i * vel1_i * l/Re_lit # Define the fluid's dynamic viscosity

# Numerical setup
Nx = 255
dx = 2 * L/(1.*(Nx+1))

C = 0.3
mydt = C * dx / c
Nt = 1000

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================
    'run_time_info'                : 'T',
    # ==========================================================

    # Computational Domain Parameters ==========================
    # Periodic Domain from -pi*l to pi*l
    'x_domain%beg'                 : -L, 
    'x_domain%end'                 : L,
    'y_domain%beg'                 : -L,
    'y_domain%end'                 : L,
    'm'                            : Nx,
    'n'                            : Nx,
    'p'                            : 0,
    'dt'                           : mydt,
    't_step_start'                 : 0,
    't_step_stop'                  : math.ceil(Nt),
    't_step_save'                  : math.ceil(Nt/10),
    # ==========================================================

    # Simulation Algorithm Parameters ==========================
    'num_patches'                  : 1,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'weno_Re_flux'                 : 'T',
    'weno_avg'                     : 'T',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'perturb_flow'                 : 'T',
    'perturb_flow_fluid'           : 1,
    'bc_x%beg'                     : -1,
    'bc_x%end'                     : -1,
    'bc_y%beg'                     : -1,
    'bc_y%end'                     : -1,
    # ==========================================================

    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'F',
    'omega_wrt(3)'                 :'T',
    'fd_order'                     :'1',
    'E_wrt'                        :'T',
    # ==========================================================

    # Patch 1: Base ============================================
    # Uncomment the configuration of the Taylor Green Vortex in 
    # 2D analytical patch under m_patch.f90
    'patch_icpp(1)%geometry'       : 20,
    'patch_icpp(1)%x_centroid'     : 0.0,
    'patch_icpp(1)%y_centroid'     : 0.0,
    'patch_icpp(1)%length_x'       : 2.0*L, 
    'patch_icpp(1)%length_y'       : 2.0*L,  
    'patch_icpp(1)%vel(1)'         : vel1_i, # Define the characteristic velocity of the vortex
    'patch_icpp(1)%vel(2)'         : l, # Define the characteristic length of the vortex
    'patch_icpp(1)%pres'           : 1.E+05,
    'patch_icpp(1)%alpha_rho(1)'   : 1.22,
    'patch_icpp(1)%alpha(1)'       : 1.,
    # ==========================================================
   

    # Fluids Physical Parameters ===============================
    'fluid_pp(1)%gamma'            : 1.E+00/(gamma-1.E+00),
    'fluid_pp(1)%pi_inf'           : 0.0,
    # Shear viscosity of STD air
    'fluid_pp(1)%Re(1)'            : 1/Mu,
    # ==========================================================
}))
# ==============================================================================