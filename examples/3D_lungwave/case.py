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

#gas Cv calculation
Ru = 8.3144598
Ra = Ru/(28.955E-3)
Cp_g = Ra*gammag/(gammag-1)
Cv_g = Cp_g/gammag
#

#material2 :: water
gammal = 5.5
Bl = 492.E+06
rhol = 996.0
c_l = 1540#1648.7
G_l = 1.0E3
Cv_l = 1816

#primitive vartiables
patmos = 101325. #pa

#problem specific variable
lambda_wave = 200.E-6

#define pulse
P_amp = 1.E+6
P_len = 45*lambda_wave                  #length of the impulse
theta = -math.pi/2          #direction of propagation 

#non-dim

#define characteristic density, length, time, stress material                   #make it liquid
rho_char = 1#rhol
length_char = 1#lambda_wave
c_char = 1#c_l                                                                    #should be liquid
time_char = 1#length_char/c_char
stress_char = 1#rho_char*c_char*c_char/gammal

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
dlengx = 10.*lambda_wave
dlengy = 1.*lambda_wave/2.
dlengz = 1.*lambda_wave/2.
Ny = 48
Nx = dlengx*Ny/dlengy
Nz = dlengz*Ny/dlengy
dx = dlengx/Nx
dy = dlengy/Ny
dz = dlengz/Nz

alphal_back = 0.99
alphag_back = 1.0 - alphal_back

alphag_lung = 0.99
alphal_lung = 1.0 - alphag_lung

interface_amp = 0.03*lambda_wave

# time stepping requirements
time_end = 5.0e-5#1.5E-04
cfl = 0.3

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
    'sim_data'                     : 'T',
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
    'stretch_x'                    : 'T',
    'a_x'                          : 4.0E+00,
    'x_a'                          : -3.*lambda_wave,
    'x_b'                          : 3.*lambda_wave,
    'loops_x'                      : 2,
    'dt'                           : dt,
    't_step_start'                 : tstart,
    't_step_stop'                  : tstop,
    't_step_save'                  : 10,#tsave,
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 2,
    'model_eqns'                   : 3,
    ####Change
    'relax'                        : 'T',
    'relax_model'                  : 5,
    'palpha_eps'                   : 1.0E-5,            #check smaller  -6/-8
    'ptgalpha_eps'                 : 0.999,
    ########Change
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
    'bc_y%beg'                     : -2,                # was -1 for all of them
    'bc_y%end'                     : -2,
    'bc_z%beg'                     : -2,
    'bc_z%end'                     : -2,
    # ==========================================================================

    # Turning on Hypoelasticity ================================================
    #'hypoelasticity'               : 'T',     
    'hyperelasticity'              : 'T',
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
    'acoustic(1)%loc(1)'            : 0.1*lambda_wave,
    'acoustic(1)%loc(2)'            : 0,#lambda_wave/2,
    'acoustic(1)%loc(3)'            : 0,
    'acoustic(1)%pulse'             : 3,
    'acoustic(1)%npulse'            : 1,
    'acoustic(1)%wavelength'        : P_len,            #wavelength of the signal
    'acoustic(1)%mag'               : P_amp_n,
    'acoustic(1)%length'            : 2*dlengy,            #length of the place ???
    'acoustic(1)%height'            : 2*dlengz,           #maybe 2dlengz
    'acoustic(1)%dir'               : -math.pi,
    #===========================================================================
    
    # Patch 1: Background ======================================================
    'patch_icpp(1)%geometry'       : 9,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%y_centroid'     : 0,#dlengy/2.,
    'patch_icpp(1)%z_centroid'     : 0,#dlengz/2.,
    'patch_icpp(1)%length_x'       : 5,#100*lambda_wave,#5*dlengx,#30*lambda_wave,
    'patch_icpp(1)%length_y'       : 2*dlengy,
    'patch_icpp(1)%length_z'       : 2*dlengz,
    'patch_icpp(1)%vel(1)'         : 0.E+00,
    'patch_icpp(1)%vel(2)'         : 0.E+00,            
    'patch_icpp(1)%vel(3)'         : 0.E+00,
    'patch_icpp(1)%pres'           : patmos_n,
    'patch_icpp(1)%alpha_rho(1)'   : rhol_n*alphal_back,
    'patch_icpp(1)%alpha_rho(2)'   : rhog_n*alphag_back,                                                        #   make non 0
    'patch_icpp(1)%alpha(1)'       : alphal_back,
    'patch_icpp(1)%alpha(2)'       : alphag_back,

    # ==========================================================================

    # Patch 2: Lung ============================================================
    'patch_icpp(2)%geometry'       : 13,
    'patch_icpp(2)%hcid'           : 302,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%x_centroid'     : 0,#-5*lambda_wave,#-dlengx/2., #-lambda_wave*5,    #
    'patch_icpp(2)%y_centroid'     : 0,#dlengy/2.,
    'patch_icpp(2)%z_centroid'     : 0,#dlengz/2.,
    'patch_icpp(2)%length_x'       : 5,#5*dlengx,#lambda_wave*30,   #
    'patch_icpp(2)%length_y'       : 2*dlengy,  
    'patch_icpp(2)%length_z'       : 2*dlengz,                 
    'patch_icpp(2)%a2'             : interface_amp,
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
    'fluid_pp(1)%cv'               : Cv_l,
    'fluid_pp(2)%gamma'            : 1.E+00/(gammag-1.E+00),
    'fluid_pp(2)%pi_inf'           : gammag*Bg_n/(gammag-1.E+00),   
    'fluid_pp(2)%G'                : G_g_n,
    'fluid_pp(2)%cv'               : Cv_g,
 #==============================================================================
}))

# ==============================================================================
