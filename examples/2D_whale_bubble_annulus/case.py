#!/usr/bin/env python2

import math
import json

# Select type of simulation
# restart_name = argv[2].strip()

# x0      = 10.E-06
x0      = 1.
p0      = 101325.
rho0    = 1000.
u0      = math.sqrt( p0/rho0 )
c0      = 1475.



n_tait  = 7.1
B_tait  = 306.E+06/p0 
gamma_gas = 1.4

pv      = 2.3388E+03 
#Cavitation number
Ca = (p0-pv)/(0.5*rho0*u0**2)

Ly =  6./x0
Lx = 6./x0

Ny = 249
Nx = Ny
dx = Lx/float(Nx)
dy = Ly/float(Ny)

# Time stepping parameters
cfl = 0.3
dt  = cfl*dx*u0/(c0)
T   = 20.
Ntfinal   = int(T/dt)
Ntrestart = int(Ntfinal/5.)

# Init
# t_start = 0
# Nfiles  = 5E1
# t_save  = int(math.ceil(Ntrestart/float(Nfiles)))
# Nt      = t_save*Nfiles
# Ntrestart = Nt
# bc_y    = 8

# if restart_name == 'run':
    # Simulate
# t_start = Ntrestart
t_start = 0
Nfiles  = 1E2
t_save  = int(math.ceil((Ntfinal-t_start)/float(Nfiles)))
Nt      = t_save*Nfiles
bc_y    = 3
# elif restart_name != 'init':
#     sys.exit("incorrect restart parameter")

ang = 1.

myr0 = 1.E+00
vf0  = 1.E-12
alf = 4.E-3
# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'F',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -Lx/2.,
    'x_domain%end'                 :  Lx/2.,
    'y_domain%beg'                 : -Ly/2.,
    'y_domain%end'                 :  Ly/2.,
    'cyl_coord'                    : 'F',
    'm'                            : Nx,
    'n'                            : Ny,
    'p'                            : 0,
    'dt'                           : dt,
    't_step_start'                 : 0,
    't_step_stop'                  : 1000,
    't_step_save'                  : 10 ,
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 3,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'adv_alphan'                   : 'T',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_order'                   : 3,
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'F',  
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'T',
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
    # ==========================================================================
                                                                
    # Patch 1 ==================================================================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%y_centroid'     : 0.,
    'patch_icpp(1)%length_x'       : Lx,
    'patch_icpp(1)%length_y'       : Ly,
    'patch_icpp(1)%alpha_rho(1)'   : (1.-vf0)*1.,
    'patch_icpp(1)%alpha(1)'       : vf0,
    'patch_icpp(1)%vel(1)'         : 0.0,
    'patch_icpp(1)%vel(2)'         : 0.00,
    'patch_icpp(1)%pres'           : 1.,
    'patch_icpp(1)%r0'             : 1.E+00,
    'patch_icpp(1)%v0'             : 0.0E+00,
    # ==========================================================================

    # Patch 2 ==================================================================
    'patch_icpp(2)%geometry'       : 2,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%x_centroid'     : 0.,
    'patch_icpp(2)%y_centroid'     : 0.,
    'patch_icpp(2)%radius'          : 1.2,
    'patch_icpp(2)%alpha_rho(1)'   : (1.-alf)*1.,
    'patch_icpp(2)%alpha(1)'       : alf,
    'patch_icpp(2)%vel(1)'         : 0.0,
    'patch_icpp(2)%vel(2)'         : 0.0,
    'patch_icpp(2)%pres'           : 1.,
    'patch_icpp(2)%r0'             : 1.E+00,
    'patch_icpp(2)%v0'             : 0.0E+00,
    # ==========================================================================

    # # ========================================================================
     'patch_icpp(3)%geometry'       : 2,
     'patch_icpp(3)%x_centroid'     : 0.,
     'patch_icpp(3)%y_centroid'     : 0.,
     'patch_icpp(3)%radius'          : 0.8,
     'patch_icpp(3)%alter_patch(1)' : 'T',
     'patch_icpp(3)%alter_patch(2)' : 'T',
     'patch_icpp(3)%alpha_rho(1)'   : (1-vf0)*1.,
     'patch_icpp(3)%vel(1)'         : 0.00,
     'patch_icpp(3)%vel(2)'         : 0.00,
     'patch_icpp(3)%pres'           : 1.,
     'patch_icpp(3)%alpha(1)'       : vf0,
     'patch_icpp(3)%r0'             : 1.E+00,
     'patch_icpp(3)%v0'             : 0.0E+00,
    # # ========================================================================

    # SHB: Bubbles =============================================================
    # 'perturb_flow'              : 'T',
    # 'perturb_flow_fluid'        : 1,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    # Surrounding liquid
    'fluid_pp(1)%gamma'             : 1.E+00/(n_tait-1.E+00),
    'fluid_pp(1)%pi_inf'            : n_tait*B_tait/(n_tait-1.),
    'fluid_pp(2)%gamma'             : 1./(gamma_gas-1.),
    'fluid_pp(2)%pi_inf'            : 0.0E+00,
    # ==========================================================================
    # SHB: Tait EOS ============================================================
    # ==========================================================================

    # Bubbles ==================================================================
    'bubbles'               : 'T',
    'bubble_model'          : 3,
    'polytropic'            : 'T',
    'R0_type'               : 1,
    'thermal'               : 3,
    'R0ref'                 : myr0,
    'nb'                    : 1,
    'Ca'                    : Ca,
    # ==========================================================================

    'Monopole'                     : 'T',
    'num_mono'                     : 1,
    'Mono(1)%loc(1)'               : -1.5,
    'Mono(1)%loc(2)'               : 0.0,
    'Mono(1)%npulse'               : 4,
    'Mono(1)%dir'                  : 0.78539816339,
    'Mono(1)%pulse'                : 1,
    'Mono(1)%mag'                  : 1.,
    'Mono(1)%length'               : 0.2,
    'rdma_mpi'                     : 'F',
}))
