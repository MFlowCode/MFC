#!/usr/bin/env python3
import argparse
import math
import json
import numpy as np


#Configuration case dictionary
data = {
    # Logistics
    'run_time_info'     : 'T',
    # Computational Domain
    'x_domain%beg'      : -0.5,
    'x_domain%end'      : 0.5,
    'y_domain%beg'      : -0.125,
    'y_domain%end'      : 0.125,
    'm'                 : 256,
    'n'                 : 64,
    'cyl_coord'        : 'F',
    # "dt"                : dt,
    "t_step_start"      : 0,
    "t_step_stop"       : 30000, #1250, #t_step_stop,
    "t_step_save"       : 200, #t_step_save,
    "adap_dt" : "F", #adaptive time stepping for bubble solver

    "cfl_adap_dt"       : "T",
    "cfl_target"        : 0.4,
    "n_start"           : 0,
    "t_stop"            : 0.002,
    "cfl_const_dt"      : "F",
    "t_save"            : 1e-5, #t_step_save,
    # Simulation Algorithm
    'model_eqns'        : 2,
    'alt_soundspeed'    : 'F',
    'mixture_err'       : 'F',
    'mpp_lim'           : 'F',
    'time_stepper'      : 3,
    'weno_order'        : 5,
    'mapped_weno'       : 'T',
    'mp_weno'           : 'T',
    'avg_state'         : 2,
    'weno_eps'          : 1e-16,
    'riemann_solver'    : 2,
    'wave_speeds'       : 1,
    'bc_x%beg'          : -16, #No slip
    'bc_x%end'          : -16, #No slip
    'bc_y%beg'          : -16, #No slip
    'bc_y%end'          : -16, #No slip
    'num_patches'       : 2,
    'num_fluids'        : 1,
    # Database Structure Parameters
    'format'            : 1,
    'precision'         : 2,
    'prim_vars_wrt'     : 'T',
    'parallel_io'       : 'T',
    'lag_db_wrt': "T",
    # Fluid Parameters Ambient gas
    "fluid_pp(1)%gamma"  : 1.4, #Sp. heat ratio
    "fluid_pp(1)%pi_inf" : 0.0, #Liquid stiffness
    "fluid_pp(1)%cv"     : 717.5,      # J/(kg·K) ,heat capacity 
    "fluid_pp(1)%qv"     : 0.0, #reference energy per unit mass for SGEOS, q (see Le Metayer (2004))
    "fluid_pp(1)%qv"     : 0.0, #reference entropy per unit mass for SGEOS, q' (see Le Metayer (2004))
    "fluid_pp(1)%G"      : 1.4-1,
    # Viscosity
    'viscous'         : 'F',
    # "fluid_pp(1)%Re(1)": 1.48e-5, #1.0 / (mu_host / (rho0 * c0 * L0)),
    # "fluid_pp(1)%Re(2)": 1.48e-5,
    # Patch for background flow
    'patch_icpp(1)%geometry'    : 3,
    'patch_icpp(1)%x_centroid'  : 0.,
    'patch_icpp(1)%y_centroid'  : 0.,
    'patch_icpp(1)%length_x'    : 1.,
    'patch_icpp(1)%length_y'    : 0.25,
    'patch_icpp(1)%vel(1)'      : 0.,
    'patch_icpp(1)%vel(2)'      : 0.,
    'patch_icpp(1)%pres'        : 1.e4, #pres,
    'patch_icpp(1)%alpha_rho(1)': 0.125, #rho_host,
    'patch_icpp(1)%alpha(1)'    : 1,
    #Patch for High Pressure Region
    'patch_icpp(2)%geometry'    : 3,
    'patch_icpp(2)%alter_patch(1)': "T",
    'patch_icpp(2)%x_centroid'  : -0.25,
    'patch_icpp(2)%y_centroid'  : 0,
    'patch_icpp(2)%length_x'    : 0.5,
    'patch_icpp(2)%length_y'    : 0.25,
    'patch_icpp(2)%vel(1)'      : 0.,
    'patch_icpp(2)%vel(2)'      : 0.,
    'patch_icpp(2)%pres'        : 1.e5,
    'patch_icpp(2)%alpha_rho(1)': 1.,
    'patch_icpp(2)%alpha(1)'    : 1,


    # Lagrangian Bubbles/Particles
    "bubbles_lagrange": "F",
    # "fd_order": 4,
    # "bubble_model": 0,  # (0) Particle (2) KM (3) RP
    # "lag_params%nBubs_glb": 0, # Number of bubbles
    # "lag_params%solver_approach": 2,
    # "lag_params%cluster_type": 2,
    # "lag_params%pressure_corrector": "T",
    # "lag_params%smooth_type": 1,
    # "lag_params%epsilonb": 1.0,
    # "lag_params%valmaxvoid": 0.9,
    # "lag_params%write_bubbles": "F",
    # "lag_params%write_bubbles_stats": "F",

    "particles_lagrange": "T",
    "fd_order": 4, #4th order can be unstable due to velocity interpolation for drag
    "lag_params%epsilonb" : 1.,
    "lag_params%charwidth": 0.00390625, #0.03515625, #virtual depth in 2D sims
    "lag_params%solver_approach": 2, #1 for one way coupling, 2 for two way coupling 
    "lag_params%nParticles_glb": 1000, # Number of Particles
    "lag_params%valmaxvoid": 0.9,
    "lag_params%vel_model": 1, #Set to 1 for moving particles, 0 for stationary
    "lag_params%smooth_type": 1, #1 gaussian, 2 delta
    "lag_params%input_path" : 'input/lag_particles.dat',
    "lag_params%qs_drag_model" : 2, # 1 = Parmar, 2 = Osnes, 3 = Modified Parmar, 4 = Gidaspow
    "lag_params%stokes_drag" : 0, # 0 = none, 1 = free slip stokes, 2 = no slip stokes
    "lag_params%pressure_force" : "F",
    "lag_params%gravity_force" : "F",
    "particle_pp%rho0ref_particle": 2700.,
    "particle_pp%cp_particle": 1000. #not currently used


}
mods = {}

print(json.dumps({**data,**mods},indent=4))
