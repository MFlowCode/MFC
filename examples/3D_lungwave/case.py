#!/usr/bin/env python3

import math
import json

pi = 3.141592653589
# material parameters

#material1 :: gas
#patterson 2018

gammag = 1.4    #unitless
Bg =0           #pascals
rhog = 1.18     #kg/m^3
c_g = 347.2     #m/s
G_g = 0         #pa

#material2 :: water
gammal = 5.5
Bl = 492.E+06
rhol = 996.0
c_l = 1648.7
G_l = 1E+06


#primitive vartiables
patmos = 101325. #pa

#problem specific variable
lambda_wave = 200.E-6

#define pulse
P_amp = 10.E+6
P_len = 45                  #length of the impulse
theta = -math.pi/2          #direction of propagation

#non-dim

#define characteristic density, length, time, stress material                   #make it liquid
rho_char = rhol
length_char = lambda_wave
c_char = c_l                                                                    #should be liquid
time_char = length_char/c_char
stress_char = rho_char*c_char*c_char/gammal

#non-dim the properties
rhog_n  = rhog/rho_char
c_g_n = c_g/c_char
rhol_n = rhol/rho_char
c_l_n = c_l/c_char
Bg_n = Bg/stress_char
Bl_n = Bl/stress_char
G_g_n = G_g/stress_char
G_l_n = G_l/stress_char
patmos_n = patmos/stress_char
P_amp_n = P_amp/stress_char

#geometry
dlengx = 15.
dlengy = 1.
dlengz = 1.
Ny = 25
Nx = dlengx*Ny
Nz = dlengz*Ny
dx = dlengx/Nx
dy = dlengy/Ny
dz = dlengz/Nz
alphal_back = 1.0
alphag_back = 0.0
alphal_lung = 0.0
alphag_lung = 1.0

interface_amp = 0.03

# time stepping requirements
time_end = 50
cfl = 0.01

dt = cfl * dx/c_l_n
Nt = int(time_end/dt)
Nframes = 500
tstart = 0
tstop = Nt
tsave = int(Nt/Nframes)

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    #'sim_data'                     : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -dlengx/2.,
    'x_domain%end'                 : dlengx/2.,
    'y_domain%beg'                 : 0.,
    'y_domain%end'                 : dlengy,
    'z_domain%beg'                 : 0.,
    'z_domain%end'                 : dlengz,
    'm'                            : int(Nx),
    'n'                            : int(Ny),
    'p'                            : int(Nz),
    'stretch_x'                    : 'F',
    'a_x'                          : 4.0E+00,
    'x_a'                          : -5.,
    'x_b'                          : 5.,
    'dt'                           : dt,
    't_step_start'                 : tstart,
    't_step_stop'                  : tstop,
    't_step_save'                  : tsave,
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 2,
    'model_eqns'                   : 3,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 2,
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
    'bc_x%beg'                     : -6,
    'bc_x%end'                     : -6,
    'bc_y%beg'                     : -1,
    'bc_y%end'                     : -1,
    'bc_z%beg'                     : -1,
    'bc_z%end'                     : -1,
    # ==========================================================================

    # Turning on Hypoelasticity ================================================
    #'hypoelasticity'               : 'T',
    'hyperelasticity'              : 'F',
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    # ==========================================================================

    # acoustics setting =========================================================
    'acoustic_source'               : 'T',
    'num_source'                    : 1,
    'acoustic(1)%support'           : 3,
    'acoustic(1)%loc(1)'            : 4,
    'acoustic(1)%loc(2)'            : dlengy/2,
    #'acoustic(1)%loc(3)'            : dlengz/2,
    'acoustic(1)%pulse'             : 3,
    'acoustic(1)%npulse'            : 1,
    'acoustic(1)%wavelength'        : P_len,            #wavelength of the signal
    'acoustic(1)%mag'               : P_amp_n,
    'acoustic(1)%length'            : dlengy,            #length of the place ???
    'acoustic(1)%height'            : dlengz,
    'acoustic(1)%dir'               : -math.pi,
    #===========================================================================

    # Patch 1: Background ======================================================
    'patch_icpp(1)%geometry'       : 9,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%y_centroid'     : dlengy/2.,
    'patch_icpp(1)%z_centroid'     : dlengz/2.,
    'patch_icpp(1)%length_x'       : 2000, #dlengx,    #
    'patch_icpp(1)%length_y'       : dlengy,
    'patch_icpp(1)%length_z'       : dlengz,
    'patch_icpp(1)%vel(1)'         : 0.E+00,
    'patch_icpp(1)%vel(2)'         : 0.E+00,
    'patch_icpp(1)%vel(3)'         : 0.E+00,
    'patch_icpp(1)%pres'           : patmos_n,
    'patch_icpp(1)%alpha_rho(1)'   : rhol_n*alphal_back,
    'patch_icpp(1)%alpha_rho(2)'   : rhog_n*alphag_back,
    'patch_icpp(1)%alpha(1)'       : alphal_back,
    'patch_icpp(1)%alpha(2)'       : alphag_back,
    #'patch_icpp(1)%tau_e(1)'       : 0.0,

    # ==========================================================================

    # Patch 2: Lung ============================================================
    'patch_icpp(2)%geometry'       : 13,
    'patch_icpp(2)%hcid'           : 301,
    #'patch_icpp(2)%geometry'       : 9,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%x_centroid'     : -998,#-dlengx/4.,    #
    'patch_icpp(2)%y_centroid'     : dlengy/2.,
    'patch_icpp(2)%z_centroid'     : dlengz/2.,
    'patch_icpp(2)%length_x'       : 2000,#dlengx,#dlengx/2.+2,   #
    'patch_icpp(2)%length_y'       : dlengy,
    'patch_icpp(2)%length_z'       : dlengz,
    'patch_icpp(2)%a(2)'           : interface_amp,
    'patch_icpp(2)%vel(1)'         : 0.E+00,
    'patch_icpp(2)%vel(2)'         : 0.0,
    'patch_icpp(2)%vel(3)'         : 0.0,
    'patch_icpp(2)%pres'           : patmos_n,
    'patch_icpp(2)%alpha_rho(1)'   : rhol_n*alphal_lung,
    'patch_icpp(2)%alpha_rho(2)'   : rhog_n*alphag_lung,
    'patch_icpp(2)%alpha(1)'       : alphal_lung,
    'patch_icpp(2)%alpha(2)'       : alphag_lung,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(gammal-1.E+00),
    'fluid_pp(1)%pi_inf'           : gammal*Bl_n/(gammal-1.E+00),
    'fluid_pp(1)%G'                : G_l_n,
    'fluid_pp(2)%gamma'            : 1.E+00/(gammag-1.E+00),
    'fluid_pp(2)%pi_inf'           : gammag*Bg_n/(gammag-1.E+00),
    'fluid_pp(2)%G'                : G_g_n,
 #==============================================================================
}))

# ==============================================================================
