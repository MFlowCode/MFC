#!/usr/bin/env python3
import argparse
import math
import json
import numpy as np

# Domain extents
xb = -2
xe = 2
yb = -2
ye = 2
zb = -2
ze = 2

# Reference values for nondimensionalization
L0 = 1e-3 # length - m (min bubble radius)
rho0 = 950  # density - kg/m3
c0 = 1449.0  # speed of sound - m/s
p0 = rho0 * c0 * c0  # pressure - Pa
T0 = 298  # temperature - K

# Host properties (water)
gamma_host = 6.12  # Specific heat ratio
pi_inf_host = 3.43e8  # Stiffness - Pa
mu_host = 0.001
c_host = 1449.0  # speed of sound - m/s
rho_host = 950 # density kg/m3
T_host = 298  # temperature K

# Lagrangian bubbles' properties
R_uni = 8314  # Universal gas constant - J/kmol/K
MW_g = 28.0  # Molar weight of the gas - kg/kmol
MW_v = 18.0  # Molar weight of the vapor - kg/kmol
gamma_g = 1.4  # Specific heat ratio of the gas
gamma_v = 1.333  # Specific heat ratio of the vapor
pv = 2350  # Vapor pressure of the host - Pa
cp_g = 1.0e3  # Specific heat of the gas - J/kg/K
cp_v = 2.1e3  # Specific heat of the vapor - J/kg/K
k_g = 0.025  # Thermal conductivity of the gas - W/m/K
k_v = 0.02  # Thermal conductivity of the vapor - W/m/K
diffVapor = 2.5e-5  # Diffusivity coefficient of the vapor - m2/s
sigBubble = 0.020  # Surface tension of the bubble - N/m
mu_g = 1.48e-5
rho_g = 1 # density kg/m3

# Nondimmensionalization of domain size
xb = xb / L0
xe = xe / L0
yb = yb / L0
ye = ye / L0
zb = zb / L0
ze = ze / L0

# patm = 1.0e05  # Atmospheric pressure - Pa
patm = 1e5
g0 = 9.81 / (c0*c0/L0)

# Patch prim vars
rho_host = rho_host / rho0
rho_g = rho_g / rho0
pres = patm / p0

# Timing
tend = 0.01 * c0 / L0
tsave = tend / 50 # save time - sec
dt = (tend / 2000) * c0 / L0

#Configuration case dictionary
data = {
    # Logistics
    'run_time_info'     : 'T',
    # Computational Domain
    'x_domain%beg'      : xb,
    'x_domain%end'      : xe,
    'y_domain%beg'      : yb,
    'y_domain%end'      : ye,
    'z_domain%beg'      : zb,
    'z_domain%end'      : ze,
    'm'                 : 99,
    'n'                 : 99,
    'p'                 : 99,
    "dt"                : 7.5,
    "t_step_start"      : 0,
    "t_step_stop"       : 2000,
    "t_step_save"       : 25,
    # Simulation Algorithm
    'model_eqns'        : 2,
    'alt_soundspeed'    : 'F',
    'mixture_err'       : 'F',
    'mpp_lim'           : 'T',
    'time_stepper'      : 3,
    'weno_order'        : 5,
    'mapped_weno'       : 'T',
    'mp_weno'           : 'F',
    'avg_state'         : 2,
    'weno_eps'          : 1e-16,
    'riemann_solver'    : 2,
    'wave_speeds'       : 1,
    'bc_x%beg'          : -1,
    'bc_x%end'          : -1,
    'bc_y%beg'          : -1,
    'bc_y%end'          : -1,
    'bc_z%beg'          : -1,
    'bc_z%end'          : -1,
    'num_patches'       : 2,
    'num_fluids'        : 2,
    # Database Structure Parameters
    'format'            : 1,
    'precision'         : 2,
    'prim_vars_wrt'     : 'T',
    'parallel_io'       : 'T',
    'lag_db_wrt' : 'T',
    'lag_rad_wrt' : 'T',
    # Fluid Parameters Host
    "fluid_pp(1)%gamma": 1.0 / (gamma_host - 1.0),
    "fluid_pp(1)%pi_inf": gamma_host * (pi_inf_host / p0) / (gamma_host - 1.0),
    "fluid_pp(1)%mul0": mu_host,
    "fluid_pp(1)%ss": sigBubble,
    "fluid_pp(1)%pv": pv,
    "fluid_pp(1)%gamma_v": gamma_v,
    "fluid_pp(1)%M_v": MW_v,
    "fluid_pp(1)%k_v": k_v,
    "fluid_pp(1)%cp_v": cp_v,
    # Fluid Parameters Gas
    "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
    "fluid_pp(2)%pi_inf": 0.0e00,
    "fluid_pp(2)%gamma_v": gamma_g,
    "fluid_pp(2)%M_v": MW_g,
    "fluid_pp(2)%k_v": k_g,
    "fluid_pp(2)%cp_v": cp_g,
    # Viscosity
    'viscous'         : 'T',
    "fluid_pp(1)%Re(1)": 1.0 / (mu_host / (rho0 * c0 * L0)),
    "fluid_pp(2)%Re(1)": 1.0 / (mu_g / (rho0 * c0 * L0)),
    # Patch for background flow
    'patch_icpp(1)%geometry'    : 9,
    'patch_icpp(1)%x_centroid'  : (xb + xe) / 2,
    'patch_icpp(1)%y_centroid'  : (yb + ye) / 2,
    'patch_icpp(1)%z_centroid'  : (zb + ze) / 2,
    'patch_icpp(1)%length_x'    : (xe - xb),
    'patch_icpp(1)%length_y'    : (ye - yb),
    'patch_icpp(1)%length_z'    : (ze - zb),
    'patch_icpp(1)%vel(1)'      : 0,
    'patch_icpp(1)%vel(2)'      : 0,
    'patch_icpp(1)%vel(3)'      : 0,
    'patch_icpp(1)%pres'        : pres,
    'patch_icpp(1)%alpha_rho(1)': rho_host,
    'patch_icpp(1)%alpha_rho(2)': rho_g,
    'patch_icpp(1)%alpha(1)'    : 1,
    'patch_icpp(1)%alpha(2)'    : 0,
    # High pressure
    'patch_icpp(2)%geometry'    : 8,
    'patch_icpp(2)%alter_patch(1)': "T",
    'patch_icpp(2)%x_centroid'  : 0, #1.15/L0,
    'patch_icpp(2)%y_centroid'  : 0, #0.33/L0,
    'patch_icpp(2)%z_centroid'  : 0, #0.33/L0,
    'patch_icpp(2)%radius'      : 0.2/L0,
    'patch_icpp(2)%vel(1)'      : 0,
    'patch_icpp(2)%vel(2)'      : 0,
    'patch_icpp(2)%vel(3)'      : 0,
    'patch_icpp(2)%pres'        : 200*pres,
    'patch_icpp(2)%alpha_rho(1)': rho_host,
    'patch_icpp(2)%alpha_rho(2)': rho_g,
    'patch_icpp(2)%alpha(1)'    : 1,
    'patch_icpp(2)%alpha(2)'    : 0,
    # Lagrangian Bubbles
    "bubbles_lagrange": "T",
    "fd_order": 4,
    "bubble_model": 0,
    "lag_params%nBubs_glb": 10000,
    "lag_params%vel_model": 2,
    "lag_params%drag_model": 1,
    "lag_params%solver_approach": 2,
    "lag_params%cluster_type": 1,
    "lag_params%pressure_corrector": "F",
    "lag_params%smooth_type": 1,
    "lag_params%epsilonb": 1.0,
    "lag_params%valmaxvoid": 0.9,
    "lag_params%write_bubbles": "T",
    "lag_params%write_bubbles_stats": "F",
    "lag_params%c0": c0,
    "lag_params%rho0": rho0,
    "lag_params%T0": T0,
    "lag_params%x0": L0,
    "lag_params%diffcoefvap": diffVapor,
    "lag_params%Thost": T_host,
}

mods = {}

print(json.dumps({**data,**mods},indent=4))
