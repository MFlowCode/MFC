#!/usr/bin/env python3

import math
import json

# Overview: A planar acoustic wave interacts with a bubble cloud in water.

# Reference length - m
x0 = 1.e-03

# Reference density - kg/m3
rho0 = 1.e+03

# Reference speed of sound - m/s
c0 = 1475.

# Reference pressure - Pa
p0 = rho0*c0*c0

# Reference temperature - K
T0 = 298

# Gamma
Gam = 7.1

# pi infinity - Pa
pi_inf = 306.e+06

# Atmospheric pressure - Pa
patm = 101325.


# Amplitude of the sinusoidal acoustic wave - Pa
pamp = 2*(1.e5)

# Frequency of the sinusoidal acoustic wave - Hz
freq = 300e+03

# Wavelength of the sinusoidal acoustic wave - m
wlen = c0/freq


## Properties that govern the dynamics of the lagrangian bubbles

# Universal gas constant - kJ/mol/K
R = 8314

# gamma, gas and vapor
gamma_g = 1.4
gamma_v = 1.333

# vapor pressure of water - Pa
pv = 2350

# cp, gas and vapor - kJ/g/K
cp_g = 1.e3
cp_v = 2.1e3

# thermal conductivity, gas and vapor - W/m/K
k_g = 0.025
k_v = 0.02

# Molar weigth, gas and vapor - g/mol
MW_g = 28.0
MW_v = 18.0

# Diffusivity coefficient of the vapor
diffVapor = 2.5e-5

# Surface tension of the bubble - N/m
sigBubble = 0.069

# Water viscosity - Pa.s
mu_water = 1e-3


## Domain boundaries - m

# x direction
xb = -2.5e-3
xe = 2.5e-3

# y direction
yb = -2.5e-3
ye = 2.5e-3

# z direction
zb = -12.e-3
ze = 12.e-3

# number of elements into x direction
Nx = 50

# number of elements into y direction
Ny = 50

# number of elements into z direction
Nz = 500

# time-step - sec
dt = 7.5e-9

# ==============================================================================

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================
    'run_time_info'                : 'T',
    # ==========================================================

    # Computational Domain Parameters ==========================
    'x_domain%beg'                 : xb/x0,
    'x_domain%end'                 : xe/x0,
    'y_domain%beg'                 : yb/x0,
    'y_domain%end'                 : ye/x0,
    'z_domain%beg'                 : zb/x0,
    'z_domain%end'                 : ze/x0,
    'stretch_z'                    : 'T',
    'a_z'                          : 2.0E-00,
    'z_a'                          : -2.5,
    'z_b'                          : 2.5,
    'stretch_y'                    : 'F',
    'stretch_x'                    : 'F',
    'm'                            : Nx,
    'n'                            : Ny,
    'p'                            : Nz,
    'dt'                           : dt*c0/x0,
    't_step_start'                 : 0,
    't_step_stop'                  : 3000,
    't_step_save'                  : 500,
    # ==========================================================

    # Simulation Algorithm Parameters ==========================
    'model_eqns'                   : 2,
    'num_fluids'                   : 1,
    'num_patches'                  : 1,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'F',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.0E-16,
    'mapped_weno'                  :'T',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     :-6,
    'bc_x%end'                     :-6,
    'bc_y%beg'                     :-6,
    'bc_y%end'                     :-6,
    'bc_z%beg'                     :-6,
    'bc_z%end'                     :-6,
    # ==========================================================

    # Acoustic source ==========================================
    'Monopole'                     : 'T',
    'num_mono'                     : 1,
    'Mono(1)%support'              : 4,
    'Mono(1)%pulse'                : 1,
    'Mono(1)%npulse'               : 1,
    'Mono(1)%mag'                  : pamp/p0,
    'Mono(1)%length'               : wlen/x0,
    'Mono(1)%loc(1)'               : 0.,
    'Mono(1)%loc(2)'               : 0.,
    'Mono(1)%loc(3)'               : -7.E-03/x0,
    'Mono(1)%dir'                  : 0.,
    'Mono(1)%delay'                : 0.,
    # ==========================================================

    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'c_wrt'                        :'T',
    'parallel_io'                  :'T',
    'fd_order'                     : 1,
    'omega_wrt(3)'                 :'T',                 
    'rho_wrt'                      :'T',
    'mom_wrt'                      :'T',
    'E_wrt'                        :'T',
    'gamma_wrt'                    :'T',
    'heat_ratio_wrt'               :'T',
    'pi_inf_wrt'                   :'T',
    'probe_wrt'                    :'T',
    'num_probes'                   : 2,
    'probe(1)%x'                   : 0.,
    'probe(1)%y'                   : 0.,
    'probe(1)%z'                   : 0.,
    'probe(2)%x'                   : 0.,
    'probe(2)%y'                   : 0.,
    'probe(2)%z'                   : -6.E-03/x0,
    # ==========================================================

    # Patch 1: Water (left) ====================================
    'patch_icpp(1)%geometry'       : 9,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%y_centroid'     : 0.,
    'patch_icpp(1)%z_centroid'     : 0.,
    'patch_icpp(1)%length_x'       : 2*(xe-xb)/x0,
    'patch_icpp(1)%length_y'       : 2*(ye-yb)/x0,
    'patch_icpp(1)%length_z'       : 2*(ze-zb)/x0,
    'patch_icpp(1)%vel(1)'         : 0.,
    'patch_icpp(1)%vel(2)'         : 0.,
    'patch_icpp(1)%vel(3)'         : 0.,
    'patch_icpp(1)%pres'           : patm/p0,
    'patch_icpp(1)%alpha_rho(1)'   : 1.,
    'patch_icpp(1)%alpha(1)'       : 1.,
    # ==========================================================

    # Lagrangian Particles (bubbles) ===========================
     'particleflag'                : 'T',
     'avgdensFlag'                 : 'T',
     'particleoutFlag'             : 'T',
     'particlestatFlag'            : 'F',
     'RPflag'                      : 'T',
     'clusterflag'                 : '3',
     'stillparticlesflag'          : 'T',
     'heatflag'                    : 1,
     'massflag'                    : 1,
     'csonref'                     : c0,
     'rholiqref'                   : rho0,
     'Lref'                        : x0,
     'Tini'                        : T0,
     'Runiv'                       : R,
     'gammagas'                    : gamma_g,
     'gammavapor'                  : gamma_v,
     'pvap'                        : pv,
     'cpgas'                       : cp_g,
     'cpvapor'                     : cp_v,
     'kgas'                        : k_g,
     'kvapor'                      : k_v,
     'MWgas'                       : MW_g,
     'MWvap'                       : MW_v,
     'diffcoefvap'                 : diffVapor,
     'sigmabubble'                 : sigBubble,
     'viscref'                     : mu_water,
     'RKeps'                       : 1.E-05,
     'ratiodt'                     : 1,
     'projectiontype'              : 0,
     'smoothtype'                  : 1,
     'epsilonb'                    : 1.33,
     'coupledFlag'                 : 'T',
     'solverapproach'              : 2,
     'correctpresFlag'             : 'T',
     'charwidth'                   : 5.0,
     'valmaxvoid'                  : 0.9,
     'dtmaxpart'                   : 0.0,
    # ==========================================================

    # Fluids Physical Parameters ===============================
    'fluid_pp(1)%gamma'            : 1.0/(Gam-1.0),
    'fluid_pp(1)%pi_inf'           : Gam*(pi_inf/p0)/(Gam-1.0),
    # ==========================================================
 }))

# ==============================================================================
