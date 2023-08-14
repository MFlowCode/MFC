import math
import json

# Parameters
epsilon = '5d0'
alpha = '1d0'
gamma = '1.4'

# Initial conditions
vel1_i = '0d0'
vel2_i = '0d0'
T_i = '1d0'
pres_i = '1d0'
alpha_rho1_i = '1d0'

# Perturbations
vel1 = f'{vel1_i} + (y - yc)*({epsilon}/(2d0*pi))*' + \
    f'exp({alpha}*(1d0 - (x - xc)**2d0 - (y - yc)**2))'
vel2 = f'{vel2_i} - (x - xc)*({epsilon}/(2d0*pi))*' + \
    f'exp({alpha}*(1d0 - (x - xc)**2d0 - (y - yc)**2))'
T = f'{T_i} - (({gamma} - 1d0))/(16d0*{alpha}*{gamma}*pi**2)*' + \
    f'exp(2*{alpha}*(1d0 - (x - xc)**2 - (y - yc)**2))'
alpha_rho1 = f'{T}**(1d0/({gamma} - 1d0))'
pres = f'{alpha_rho1} ** {gamma}'

# Numerical setup
Nx = 399
dx = 1./(1.*(Nx+1))

c = 1.4**2
C = 0.3
mydt = C * dx / c
Nt = 1

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -5.,
    'x_domain%end'                 : 5.,
    'y_domain%beg'                 : -5.,
    'y_domain%end'                 : 5.,
    'm'                            : Nx,
    'n'                            : Nx,
    'p'                            : 0,
    'dt'                           : mydt,
    't_step_start'                 : 0,
    't_step_stop'                  : int(Nt),
    't_step_save'                  : int(Nt),
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 1,
    'model_eqns'                   : 3,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 3,
    'weno_order'                   : 3,
    'weno_eps'                     : 1.E-16,
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -1,
    'bc_x%end'                     : -1,
    'bc_y%beg'                     : -1,
    'bc_y%end'                     : -1,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    'omega_wrt(3)'                 :'T',
    'fd_order'                     : 2,
    # ==========================================================================

    # Patch 1 ==================================================================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : 0,
    'patch_icpp(1)%y_centroid'     : 0,
    'patch_icpp(1)%length_x'       : 10.,
    'patch_icpp(1)%length_y'       : 10.,
    'patch_icpp(1)%vel(1)'         : vel1,
    'patch_icpp(1)%vel(2)'         : vel2,
    'patch_icpp(1)%pres'           : pres,
    'patch_icpp(1)%rho'            : 1,
    'patch_icpp(1)%alpha_rho(1)'   : alpha_rho1,
    'patch_icpp(1)%alpha(1)'       : 1.,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(1.4-1.E+00),
    'fluid_pp(1)%pi_inf'           : 0.0,
    # ==========================================================================
}))

# ==============================================================================
