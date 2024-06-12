#!/usr/bin/env python3

import json
import math

# Adjustung the parameters for the 2D case of U-Sound-lung tissue interaction
"""
 TODO's
1. Runtime Parameters
    (done) run_time_info
    rdma_mpi
2. Computational Domain Parameters
    (done) domain - [x,y] = [{0,1},{-20,15}]
    streching?
    (done) gridcells
    (done) dt
    (done) t_step_start
    (done) t_step_stop
    t_step_save
    t_step_print
3. Patch Parameters
    (done) geometry, density and volumetric fraction
4. Immersed Boundary Patches
    geometry
5. Fluid Material's Parameters
    Re(1) - sheer viscosity of the fluid (5 eq model only)
    Re(2) - volume viscosity of the fluid (5 eq model only)
    sigma - surface tension
6. Simulation Algorithm Parameters
    BC
    (done) model eq - eq.model
    TBD
7. Formatted Database and Structure Parameters
    TBD
8. (Optional) Acoustic Source Parameters
    TBD - might be useful to implemetn US imaging
9. (Optional) Ensemble-Averaged Bubble Model Parameters
10. (Optional) Velocity Field Setup Parameters
11. (Optional) Phase Change Parameters
12. (Optional) Artificial Mach Number Parameters
"""

#Define problem specific variables

    #refference values (water, characteristic lenght)
rho_0 = 996.                #kg/m3
c_0 = 1648.7                #m/s spped of sound
l_0 = 200.E-6               #length scale um
p_0 = rho_0*c_0*c_0         #characteristic pressure

    #define non-dim
N = 100                     #points per l
dx = 1/(N-1)                #dx of the grid
l_t = 15                    #y-scaling
l_b = 20                    #y-scaling

    #dimensional parameters (air at 300K)
p_atm = 101325              #Pa
    #air
rho_a = 1.18/rho_0          #density air
c_a = 347.2/c_0             #speed of sound in water
n_a = 1.4                   #stiffened EoS constant
B_a = 0                     #stiffened EoS constant
    #water
rho_w = 996./rho_0          #density water
c_w = 1648.7/c_0            #spped of sound in air
n_w = 5.5                   #stiffened EoS constant
B_w = 492.E+6/p_0           #stiffened EoS constant

    #time settings (followed 2D_whale_bubble_annulus)
cfl = 0.25                  #cfl condition
t_char = l_0/c_0            #s characteristic timescale
dt = cfl*t_char             #s time step
L = 1000*l_0                 #m total distance travelled by the wave
Tfinal = L/c_0              #s final time of the simulation
N_steps = int(Tfinal/dt)    #number of steps the simulation will run (4000)

'''    
    #time settings (followed 2D_whale_bubble_annulus)
u0 = math.sqrt(p_atm/rho_w) #refference velocity
cfl = 0.25                  #cfl condition
dt = cfl*dx*u0/c_w          #time step calculated
Tfinal = 5                  #final time of the simulation
N_steps = int(Tfinal/dt)    #number of steps the simulation will run
'''

    #DUS settings
P_amp = 10.E+6/p_0
P_len = 45                  #length of the impulse
theta = -math.pi/2          #direction of propagation  
    
    #membrane setting
a_0 = 0.03                  #amplitude
y_l_top = f"{l_t}-{a_0}*sin(2*pi*x/{1}-pi/2)"
#y_l_bot = f"{l_b}+{a_0}*sin(2*pi*x/{1}-pi/2)"


# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================

    # The domain consists of 2 regions with a membrane inwetween. The origin is
    # set on the left end of the membrane. The domain of interest. y: -20l to 15l
    'x_domain%beg'                 : 0,                                       # update
    'x_domain%end'                 : 1,                                         # update    
    'y_domain%beg'                 : -l_b,                                    # update
    'y_domain%end'                 : l_t,                                     # update
        # Grid stretching is used in the all coordinate directions 
        # to minimize computational costs. The grid is coarsened
        # away from the bubble / origin
    #might not need
    'stretch_x'                    : 'F',
    'a_x'                          : 4.E+00,
    'x_a'                          : -1.5E-03/1.E-03,
    'x_b'                          : 1.5E-03/1.E-03,
    'stretch_y'                    : 'F',
    'a_y'                          : 4.E+00,
    'y_a'                          : -1.5E-03/1.E-03,
    'y_b'                          : 1.5E-03/1.E-03,
    # grid sells setup - m=x,n=y,p=z
    'm'                            : int(N-1),                                  # update 
    'n'                            : int((l_t+l_b)*N-1),                      # update
    # time setup
    'dt'                           : dt,                                        # update : see above
    't_step_start'                 : 0,                                         # update : start at 0
    't_step_stop'                  : N_steps,                                   # update : number of iterations
    't_step_save'                  : 10,
    # ==========================================================================
    
    # Simulation Algorithm Parameters ==========================================
    # Only two patches are necesssary, the liquid(tissue) and the
    # gas(lung)
    'num_patches'                  : 2,                                         # update: Tissue and Lung
    # Use the 5 equation model
    'model_eqns'                   : 2,                                         # update: number of equations is 5
        # 6 equations model does not need the K \div(u) term                    
    'alt_soundspeed'               : 'F',

    # num_fluids defines the total number of fluids defined in each of the 
    # patches.
    'num_fluids'                   : 2,                                         # update: each patch has 2 fluid associated with it
        # Advect both volume fractions
    'adv_alphan'                   : 'T',
        # Ensure the volume fractions sum to unity at the end of each
        # time step
    'mpp_lim'                      : 'T',
        # Correct errors when computing speed of sound
    'mixture_err'                  : 'T',
        # Use TVD RK3 for time marching
    'time_stepper'                 : 3,                                         # update: tvdrk3
        # Use WENO5
    'weno_order'                   : 5,                                         # update
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'F',  
    'weno_avg'                     : 'F',
    'avg_state'                    : 2,
        # Use the mapped WENO weights to maintain monotinicity
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
        # Use the HLL  Riemann solver
    'riemann_solver'               : 1,                                          # update : 1=HLL
    'wave_speeds'                  : 1,
    
    # We will use symmetric BC at the x-boundaries. THe bottom booundary needs
    # to be at a zero gradient. Top - non-reflective boundary conditions
    'bc_x%beg'                     : -1,                                        # update : used periodec
    'bc_x%end'                     : -1,                                        # update : used periodic
    'bc_y%beg'                     : -6,                                        # update : used non-reflecting subsonic buffer
    'bc_y%end'                     : -6,                                        # update : used non-reflecting subsonic buffer
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    # Export primitive variables in double precision with parallel
    # I/O to minimize I/O computational time during large simulations
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    # ==========================================================================

        # For now setting the the patches geometry to be rectangles with a 
        # varible y lenght:
        # Patch 1: yl = 15-a_0*sin(2*pi*x/l-pi/2)
        # Patch 2: yb = 20+a_0*sin(2*pi*x/l-pi/2)

    # Patch 1: Air (Lung) ======================================================
    # Specify the gas grid geometry
    'patch_icpp(1)%geometry'       : 3,                                         # update : assumes to be a rectangle
    'patch_icpp(1)%x_centroid'     : 1/2,                                       # update : x_centroid = 0.5l
    'patch_icpp(1)%y_centroid'     : (l_t-l_b)/2,                                  # update : y_centroid = -10l
    'patch_icpp(1)%length_x'       : 1,                                         # update : x_l = l
    'patch_icpp(1)%length_y'       : l_t+l_b,                                   # update : y_l = function of position
    # Specify the patch primitive variables 
    'patch_icpp(1)%vel(1)'         : 0.E+00,                                    # update : no initial velocity
    'patch_icpp(1)%vel(2)'         : 0.E+00,                                    # update : no initial velocity
    'patch_icpp(1)%pres'           : p_atm/p_0,                                     # update : assume atmospheric pressure
    'patch_icpp(1)%alpha_rho(1)'   : 0.E+00,                                     # update : Partial density of fluid 1 in patch 2 (no water)
    'patch_icpp(1)%alpha_rho(2)'   : rho_a,                                    # update : Partial density of fluid 2 in patch 2 (air only)
    'patch_icpp(1)%alpha(1)'       : 0.E+00,                                    # update : volume fraction of fluid 1 in patch 2 (no water)       
    'patch_icpp(1)%alpha(2)'       : 1.E+00,                                    # update : volume fraction of fluid 2 in patch 2 (air only)
    # ==========================================================================

    # Patch 2: Water(Tisue) ====================================================
    # Specify the water background grid geometry
    'patch_icpp(2)%geometry'       : 3,                                         # update : assumes to be a rectangle
    #'patch_icpp(2)%hcid'           : 205,                                       # update : hardcoded geometry of the patch
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%x_centroid'     : 1/2,                                       # update : x_centroid = 0.5l
    'patch_icpp(2)%y_centroid'     : l_t/2,                                     # update : y_centorid = 7.5l
    'patch_icpp(2)%length_x'       : 1,                                         # update : x_l = l
    'patch_icpp(2)%length_y'       : l_t,                                        # update : y_l = fucntion of position, y_l_top - didnot accept str
    # Specify the patch primitive variables
    'patch_icpp(2)%vel(1)'         : 0.E+00,                                    # update : no initial velocity
    'patch_icpp(2)%vel(2)'         : 0.E+00,                                    # update : no initial velocity
    'patch_icpp(2)%pres'           : p_atm/p_0,                                     # update : assume atmospheric pressure
    'patch_icpp(2)%alpha_rho(1)'   : rho_w,                                     # update : Partial density of fluid 1 in patch 1 (water only)
    'patch_icpp(2)%alpha_rho(2)'   : 0.E+00,                                    # update : Partial density of fluid 2 in patch 1 (no air)
    'patch_icpp(2)%alpha(1)'       : 1.E+00,                                    # update : volume fraction of fluid 1 in patch 1 (water only)
    'patch_icpp(2)%alpha(2)'       : 0.E+00,                                    # update : volume fraction of fluid 2 in patch 1 (no air)  
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    # Fluid 1 - water ; Fluid 2 - air
    'fluid_pp(1)%gamma'            : 1.E+00/(n_w-1.E+00),                       # update : refferenced 2D_whale_bubble_annulus example for the formula
    'fluid_pp(1)%pi_inf'           : n_w*B_w/(n_w-1.E+00),                      # update : refferenced 2D_whale_bubble_annulus example for the formula
    'fluid_pp(2)%gamma'            : 1.E+00/(n_a-1.E+00),                       # update : refferenced 2D_whale_bubble_annulus example for the formula
    'fluid_pp(2)%pi_inf'           : 0.E+00,                                    # update : air has liquid stiffness = 0
    # ==========================================================================

    # Acoustic Wave source =====================================================
    # The acoustic wave is placed at y = 15, ?at each node along the boundary?
    'Monopole'                      : 'T',                                      # update : creating an acoustic wave
    'num_mono'                      : 1,                                        # update : place in the middle and expand
    'Mono(1)%pulse'                 : 3,                                        # update : square wave
    'Mono(1)%npulse'                : 1,                                        # update : 1 impulse
    'Mono(1)%mag'                   : P_amp,                                    # update : magnitude
    'Mono(1)%length'                : P_len,                                    # update : impulse length
    'Mono(1)%support'               : 2,                                        # update : 2D semi infinite plane (x: -inf,inf; y:-len/2, len/2)
    'Mono(1)%loc(1)'                : 0.5,                                    # update : x_center of the domain
    'Mono(1)%loc(2)'                : 15,                                     # update : upper boundary of the domain
    'Mono(1)%dir'                   : theta,                                    # update : direction: -pi/2
    'Mono(1)%support_width'         : 49,                                       # update : 49 cells in each direction

}))

# ==============================================================================