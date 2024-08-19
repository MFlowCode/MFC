#!/usr/bin/env python3

import math
import json

N = 256 

Re = 1600
L = 1
P0 = 101325
rho0 = 1
C0 = math.sqrt(1.4*P0)
V0 = 0.1*C0
mu = V0*L/Re

cfl = 0.5
dx = 2*math.pi*L/(N+1)

dt = cfl * dx / (C0)

tC = L/V0
tEnd = 20*tC

Nt = int(tEnd/dt)


# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================
    'run_time_info'                : 'T',
    # ==========================================================

    # Computational Domain Parameters ==========================
    'x_domain%beg'                 : -math.pi*L,
    'x_domain%end'                 : math.pi*L,
    'y_domain%beg'                 : -math.pi*L,
    'y_domain%end'                 : math.pi*L,
    'z_domain%beg'                 : -math.pi*L,
    'z_domain%end'                 : math.pi*L,
    'm'                            : N,
    'n'                            : N,
    'p'                            : N,
    'cyl_coord'                    : 'F',
    'dt'                           : dt,
    't_step_start'                 : 13529,
    't_step_stop'                  : Nt,
    't_step_save'                  : int(Nt/100),
    # ==========================================================

    # Simulation Algorithm Parameters ==========================
    'num_patches'                  : 1,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.0E-16,
    'weno_Re_flux'                 : 'F',
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'T',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'			           : 2,
    'bc_x%beg'                     : -1,
    'bc_x%end'                     : -1,
    'bc_y%beg'                     : -1,
    'bc_y%end'                     : -1,
    'bc_z%beg'                     : -1,
    'bc_z%end'                     : -1,
    # ==========================================================

    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 2,
    # 'prim_vars_wrt'                :'T',
    'omega_wrt(1)'                 :'T',
    'omega_wrt(2)'                 :'T',
    'omega_wrt(3)'                 :'T',
    'qm_wrt'                       :'T',
    'fd_order'                     : 4,
    'parallel_io'                  :'T',
    # ==========================================================

    # I will use 1 for WATER properties, and 2 for AIR properties
    # Patch 1: Background (AIR - 2) =============================
    'patch_icpp(1)%geometry'       : 9,
    'patch_icpp(1)%x_centroid'     : 0,
    'patch_icpp(1)%y_centroid'     : 0,
    'patch_icpp(1)%z_centroid'     : 0,
    'patch_icpp(1)%length_x'       : 2*math.pi*L,
    'patch_icpp(1)%length_y'       : 2*math.pi*L,
    'patch_icpp(1)%length_z'       : 2*math.pi*L,
    'patch_icpp(1)%vel(1)'         : f"{V0}*sin(x/{L})*cos(y/{L})*sin(z/{L})",
    'patch_icpp(1)%vel(2)'         : f"-{V0}*cos(x/{L})*sin(y/{L})*sin(z/{L})",
    'patch_icpp(1)%vel(3)'         : 0,
    'patch_icpp(1)%pres'           : f"{P0} + ({rho0}*{V0}**2/16)*(cos(2*x/{L}) + cos(2*y/{L}))*(cos(2*z/{L}) + 2)",
    'patch_icpp(1)%alpha_rho(1)'   : 1,
    'patch_icpp(1)%alpha(1)'       : 1,
    # ==========================================================

    # Fluids Physical Parameters ===============================
    'fluid_pp(1)%gamma'            : 1.0E+00/(1.4-1),
    'fluid_pp(1)%pi_inf'           : 0,
    'fluid_pp(1)%Re(1)'            : 1/mu,
    # ==========================================================
}))

# ==============================================================================
