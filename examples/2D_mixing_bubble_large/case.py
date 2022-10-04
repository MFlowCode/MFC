#!/usr/bin/env python2

import math
import json

# x0 = 10.E-06
x0, p0, rho0, c0 = 1.E-3, 101325., 1000., 1475.
u0 = math.sqrt( p0/rho0 )

n_tait    = 7.1
B_tait    = 306.E+06 / p0
gamma_gas = 1.4

# Velocity
#udim = 4. # (used when length scale x0 -> 1)
udim = 4.*x0
uu = udim / u0

#Cavitation number
Ca = 1.

Ly =  0.5/x0
Lx = Ly*2/x0

Ny     = 79
Nx     = Ny*2+1
dx, dy = x0*Lx/float(Nx), x0*Ly/float(Ny)

# Select type of simulation
restart_name = 'init'

Ntend   = 1E5
if restart_name == 'init':
    # Init
    t_start = 0
    Nt      = Ntend
    Nfiles  = 1E1
    bc_y    = 8
elif restart_name == 'run':
    # Simulate
    t_start = Ntend
    Nt      = 5E5
    Nfiles  = 1E2
    bc_y    = 5
else:
    sys.exit("incorrect restart parameter")

# Time stepping parameters
cfl     = 0.7
mydt    = cfl*dx/(c0/u0)
Ttot    = mydt*Nt

ang = 91.

myr0 = 1.E+00
vf0  = 1.E-12

# ==============================================================================

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================
    'case_dir'                     : '\'.\'',
    'run_time_info'                : 'T',
    # ==========================================================

    # Computational Domain Parameters ==========================
    'x_domain%beg'                 : -Lx/2.,
    'x_domain%end'                 :  Lx/2.,
    'y_domain%beg'                 : -Ly/2.,
    'y_domain%end'                 :  Ly/2.,
    'cyl_coord'                    : 'F',
    'm'                            : Nx,
    'n'                            : Ny,
    'p'                            : 0,
    'dt'                           : mydt,
    't_step_start'                 : t_start,
    't_step_stop'                  : int(Nt),
    't_step_save'                  : int(math.ceil((Nt-t_start)/float(Nfiles))),
    # ==========================================================

    # Simulation Algorithm Parameters ==========================
    'num_patches'                  : 2,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 3,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 3,
    'weno_vars'                    : 2,
    'weno_order'                   : 3,
    'weno_eps'                     : 1.E-16,
    'mapped_weno'                  : 'F',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'weno_Re_flux'                 : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -1,
    'bc_x%end'                     : -1,
    'bc_y%beg'                     : -bc_y,
    'bc_y%end'                     : -bc_y,
    # ==========================================================

    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'F',
    'probe_wrt'                    :'T',
    'num_probes'                   : 1,
    'probe(1)%x'                   : 0.,
    'probe(1)%y'                   : 0.,
    'fd_order'                     : 1,
    # 'schlieren_wrt'                :'F',
    # ==========================================================
                                                                
    # Patch 1 ==================================================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%y_centroid'     : 0.,
    'patch_icpp(1)%length_x'       : Lx,
    'patch_icpp(1)%length_y'       : Ly,
    'patch_icpp(1)%alpha_rho(1)'   : (1.-vf0)*1.,
    'patch_icpp(1)%alpha_rho(2)'   : (1.-vf0)*1.E-12,
    'patch_icpp(1)%alpha_rho(3)'   : vf0*1.E-3,
    'patch_icpp(1)%alpha(1)'       : (1-vf0)*1.,
    'patch_icpp(1)%alpha(2)'       : (1-vf0)*1.E-12,
    'patch_icpp(1)%alpha(3)'       : vf0,
    'patch_icpp(1)%vel(1)'         : uu,
    'patch_icpp(1)%vel(2)'         : 0.00,
    'patch_icpp(1)%pres'           : 1.00,
    'patch_icpp(1)%r0'             : 1.E+00,
    'patch_icpp(1)%v0'             : 0.0E+00,
    # ==========================================================

    # Patch 2 ==================================================
    'patch_icpp(2)%geometry'       : 4,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%x_centroid'     : 0.,
    'patch_icpp(2)%y_centroid'     : 0.,
    'patch_icpp(2)%normal(1)'      : math.sin(math.pi*ang/180.),
    'patch_icpp(2)%normal(2)'      : math.cos(math.pi*ang/180.),
    'patch_icpp(2)%alpha_rho(1)'   : (1.-vf0)*1.E-12,
    'patch_icpp(2)%alpha_rho(2)'   : (1.-vf0)*1.,
    'patch_icpp(2)%alpha_rho(3)'   : vf0*1.E-3,
    'patch_icpp(2)%alpha(1)'       : (1-vf0)*1.E-12,
    'patch_icpp(2)%alpha(2)'       : (1-vf0)*1.,
    'patch_icpp(2)%alpha(3)'       : vf0,
    'patch_icpp(2)%vel(1)'         : -1.*uu,
    'patch_icpp(2)%vel(2)'         : 0.0,
    'patch_icpp(2)%pres'           : 1.0,
    'patch_icpp(2)%r0'             : 1.E+00,
    'patch_icpp(2)%v0'             : 0.0E+00,
    # 'patch_icpp(2)%normal(1)'      : 0.00624987793326E+00,
    # 'patch_icpp(2)%normal(2)'      :-0.99998046932219E+00,
    # 'patch_icpp(2)%length_x'       : Lx,
    # 'patch_icpp(2)%length_y'       : Ly/2.,
    # ==========================================================

    # # ==========================================================
    # 'patch_icpp(3)%geometry'       : 2,
    # 'patch_icpp(3)%x_centroid'     : 5.E-01,
    # 'patch_icpp(3)%y_centroid'     : 0.5,
    # 'patch_icpp(3)%radius'          : 0.1,
    # 'patch_icpp(3)%alter_patch(1)' : 'T',
    # 'patch_icpp(3)%alter_patch(2)' : 'T',
    # 'patch_icpp(3)%alpha_rho(1)'   : (1-vf0)*1.E+00,
    # 'patch_icpp(3)%vel(1)'         : 0.00,
    # 'patch_icpp(3)%vel(2)'         : 0.00,
    # 'patch_icpp(3)%pres'           : 1.00,
    # 'patch_icpp(3)%alpha(1)'       : vf0,
    # 'patch_icpp(3)%r0'             : 1.E+00,
    # 'patch_icpp(3)%v0'             : 0.0E+00,
    # # ==========================================================

    # SHB: Bubbles =============================================
    # 'perturb_flow'              : 'T',
    # 'perturb_flow_fluid'        : 1,
    # ==========================================================

    # Fluids Physical Parameters ===============================
    # Surrounding liquid
    'fluid_pp(1)%gamma'             : 1.E+00/(n_tait-1.E+00),
    'fluid_pp(1)%pi_inf'            : n_tait*B_tait/(n_tait-1.),
    'fluid_pp(2)%gamma'             : 1.E+00/(n_tait-1.E+00),
    'fluid_pp(2)%pi_inf'            : n_tait*B_tait/(n_tait-1.),
    'fluid_pp(3)%gamma'             : 1./(gamma_gas-1.),
    'fluid_pp(3)%pi_inf'            : 0.0E+00,
    # ==========================================================
    # SHB: Tait EOS ============================================
    'pref'                  : p0,
    'rhoref'                : rho0,
    # ==========================================================

    # Bubbles ==================================================
    'bubbles'               : 'T',
    'bubble_model'          : 2,
    'polytropic'            : 'T',
    # 'polydisperse'          : 'T',
    'R0_type'               : 1,
    #'polydisperse'          : 'F',
    # 'poly_sigma'            : 0.1,
    'thermal'               : 3,
    'R0ref'                 : myr0,
    # 'nb'                    : 3,
    'nb'                    : 1,
    'Ca'                    : Ca,
    # 'Web'                   : We,
    # 'Re_inv'                : Re_inv,
    #'qbmm'               : 'T',
    #'dist_type'          : 2,
    #'sigR'               : 0.1,
    #'sigV'               : 0.1,
    #'rhoRV'              : 0.0,
    # ==========================================================
}))

# ==============================================================================
