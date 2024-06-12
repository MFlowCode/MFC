#!/usr/bin/env python3

import math
import json

# material parameters
# material 1 :: gas
# Patterson 2018 et al.
gammag = 1.4  # unitless
Bg = 0.       # Pascals
rhog = 1.18   # kg/m^3
c_g = 347.2   # m/sec
Gg = 0.      # Pascals

# material 2: lung
# Patterson 2018 et al.
gammal = 5.5  # unitless
Bl = 492E+06  # Pascals
rhol = 996.0  # kg/m^3
c_l = 1648.7  # m/sec
Gl = 1E3     # Pascals, homework!

# primitive variables (if any)
patmos = 101325.0 # Pascals, at Standard temperature and pressure

# problem specific variable
lambda_wave = 1E-3 # meters

# non-dimensionalization
# define a characteristic density, length, time, and stress
rho_char = rhog
length_char = lambda_wave
vel_char = c_g
time_char = length_char/vel_char
stress_char = rho_char*vel_char*vel_char/gammag

# nondimensionalize the material properties
rhog_n = rhog/rho_char
c_g_n = c_g/vel_char
Bg_n = Bg/stress_char  
Gg_n = Gg/stress_char

rhol_n = rhol/rho_char  
c_l_n = c_l/vel_char
Bl_n = Bl/stress_char
Gl_n = Gl/stress_char

patmos_n = patmos/stress_char

# spatial geometry
dlengx = 1.0 
dlengy = 20.
Nx = 200
Ny = dlengy*Nx

dx = dlengx/Nx
dy = dlengy/Ny

alphal_back = 1.0
alphag_back = 0.0

alphal_lung = 0.0
alphag_lung = 1.0

interface_amp = 0.5

# time stepping requirements
time_end = 0.5
cfl = 0.1

dt = cfl * dx/c_l 
Nt = int(time_end/dt)
Nframes = 50000
tstart = 0
tstop = Nt
tsave = int(Nt/Nframes)

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 :  0,
    'x_domain%end'                 :  dlengx,
    'y_domain%beg'                 : -dlengy/2.,
    'y_domain%end'                 :  dlengy/2.,
    'm'                            : int(Nx),
    'n'                            : int(Ny),
    'p'                            : 0,
    'dt'                           : dt,
    't_step_start'                 : tstart,
    't_step_stop'                  : tstop,
    't_step_save'                  : tsave,
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 2,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 2,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'T',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'F',  
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'T',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -1,
    'bc_x%end'                     : -1,
    'bc_y%beg'                     : -6,
    'bc_y%end'                     : -6,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    # ==========================================================================
                                                                
    # Patch 1: Background ======================================================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : dlengx/2.0,
    'patch_icpp(1)%y_centroid'     : 0.,
    'patch_icpp(1)%length_x'       : dlengx,
    'patch_icpp(1)%length_y'       : dlengy,
    'patch_icpp(1)%vel(1)'         : 0.E+00,
    'patch_icpp(1)%vel(2)'         : 0.E+00,
    'patch_icpp(1)%pres'           : patmos_n,
    'patch_icpp(1)%alpha_rho(1)'   : rhol_n*alphal_back,
    'patch_icpp(1)%alpha_rho(2)'   : rhog_n*alphag_back,
    'patch_icpp(1)%alpha(1)'       : alphal_back,
    'patch_icpp(1)%alpha(2)'       : alphag_back,
    # ==========================================================================

    # Patch 2: Lung interface state ===================================================
    'patch_icpp(2)%geometry'       : 7,
    'patch_icpp(2)%hcid'           : 205,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%x_centroid'     : dlengx/2.0,
    'patch_icpp(2)%y_centroid'     : -dlengy/4.0,
    'patch_icpp(2)%length_x'       : dlengx,
    'patch_icpp(2)%length_y'       : dlengy/2.0 + 2.0,
    'patch_icpp(2)%a2'             : interface_amp, # this is the interface amplitude
    'patch_icpp(2)%vel(1)'         : 0.E+00,
    'patch_icpp(2)%vel(2)'         : 0.E+00,
    'patch_icpp(2)%pres'           : patmos_n,
    'patch_icpp(2)%alpha_rho(1)'   : rhol_n*alphal_lung,
    'patch_icpp(2)%alpha_rho(2)'   : rhog_n*alphag_lung,
    'patch_icpp(2)%alpha(1)'       : alphal_lung,
    'patch_icpp(2)%alpha(2)'       : alphag_lung,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(gammal-1.E+00),
    'fluid_pp(1)%pi_inf'           : gammal*Bl_n/(gammal-1.E+00),
    'fluid_pp(2)%gamma'            : 1.E+00/(gammag-1.E+00),
    'fluid_pp(2)%pi_inf'           : gammag*Bg_n/(gammag-1.E+00),
# ==============================================================================
}))

# ==============================================================================
