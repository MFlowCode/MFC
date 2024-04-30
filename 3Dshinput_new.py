import math
import json


## Deine characteristic values for the sim
Ri = 500E-6
rhol = 1051
rhog = 0.027
Pb = 3550 
Pl = 5.0E+06
#mulc = 9E-04
#mubc = 1.0E-05
domain_length = 4*Ri


## fluid properties
Pi_inf_l = 702.8E+06
Pi_inf_b = 0.0E+00
nl = 1.19
nb = 1.47
bl = 6.61E-04
bb = 0.0E+00
alpha1 = 1.0E+00
alpha2 = alpha1
cl = (nl*(Pl+Pi_inf_l)/(rhol*(1-rhol*bl)))**(0.5)

## Defining Characteristic Values 
Rc = Ri
rhoc = rhol
uc = (Pl/rhol)**(0.5E+00)
tc = Ri/uc
Pc = rhoc*uc**(2.0E+00)

## Non-Dimensionalizing values using characteristic values
Rin = Ri/Rc
rholn = rhol/rhoc
rhogn = rhog/rhoc
Pln = Pl/Pc
Pbn = Pb/Pc
tcn = tc/tc
Pi_inf_ln = Pi_inf_l/Pc
leng = domain_length/Rc
clc = cl/uc


## Non-Dimensional Numbers
#Rel = rhol*uc*2*Rin/mulc
#Reb = rhog*uc*2*Rin/mubc
Ma = uc/cl

## Grid Specifications
CFL = 0.005	 
PpBr = 48 
Nx = PpBr*leng
Ny = Nx
Nz = Nx
x_beg = 0.0E+00*leng  
x_end = 0.5E+00*leng
y_beg = 0.0E+00
y_end = 0.5E+00*leng
z_beg = 0.0E+00
z_end = 0.5E+00*leng
delta_x = leng/Nx
delta_t = CFL*delta_x/clc
Nt = int(1.0*tcn/delta_t)
x_centroidl = (x_end+x_beg)/(2.0E+00)
y_centroidl = (y_end+y_beg)/(2.0E+00)
z_centroidl = (z_end+z_beg)/(2.0E+00)
x_centroidb = 0.0E+00
y_centroidb = 0.0E+00
z_centroidb = 0.0E+00


# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================
    'run_time_info'                : 'T',
    'sim_data'                     : 'F',
    # ==========================================================

    # Computational Domain Parameters ==========================
    'x_domain%beg'                 : x_beg,
    'x_domain%end'                 : x_end,
    'y_domain%beg'                 : y_beg,
    'y_domain%end'                 : y_end,
    'z_domain%beg'                 : z_beg,
    'z_domain%end'                 : z_end,
    'm'                            : int(Nx),
    'n'                            : int(Ny),
    'p'                            : int(Nz),
    'cyl_coord'                    : 'F',
    'dt'                           : delta_t,
    't_step_start'                 : 0, 
    't_step_stop'                  : Nt,
    't_step_save'                  : 1,
# ==========================================================

    # Simulation Algorithm Parameters ==========================
    'num_patches'                  : 2,
    'model_eqns'                   : 3,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 2,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'T',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-100,
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'weno_Re_flux'                 : 'F',
    'weno_avg'                     : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -2,
    'bc_x%end'                     : -6,
    'bc_y%beg'                     : -2,
    'bc_y%end'                     : -6,
    'bc_z%beg'                     : -2,
    'bc_z%end'                     : -6,
    'stretch_x'                    : 'T',
    'stretch_y'                    : 'T',
    'stretch_z'                    : 'T',
    'a_x'                          : 4.0E+00,
    'x_a'                          : -1.5E+00,                         
    'x_b'                          : 1.5E+00,
    'a_y'                          : 4.0E+00,
    'y_a'                          : -1.5E+00,
    'y_b'                          : 1.5E+00,
    'a_z'                          : 4.0E+00,
    'z_a'                          : -1.5E+00,
    'z_b'                          : 1.5E+00,
    # ==========================================================

    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    'fd_order'			   :'1',
    'probe_wrt'			   :'T',
    'num_probes'		   : 1,
    'probe(1)%x'                   : 0.,
    'probe(1)%y'		   : 0.,
    'probe(1)%z'                   : 0.,
    # ==========================================================

    # Patch 1: Background  ============================
    'patch_icpp(1)%geometry'       : 9,
    'patch_icpp(1)%x_centroid'     : 40*x_centroidl,
    'patch_icpp(1)%y_centroid'     : 40*y_centroidl,
    'patch_icpp(1)%z_centroid'     : 40*z_centroidl,
    'patch_icpp(1)%length_x'       : 40*leng,
    'patch_icpp(1)%length_y'       : 40*leng,
    'patch_icpp(1)%length_z'       : 40*leng,
    'patch_icpp(1)%vel(1)'         : 0.E+00,
    'patch_icpp(1)%vel(2)'         : 0.E+00,
    'patch_icpp(1)%vel(3)'         : 0.E+00,
    'patch_icpp(1)%pres'           : Pln,
    'patch_icpp(1)%alpha_rho(1)'   : alpha1*rholn,
    'patch_icpp(1)%alpha_rho(2)'   : 0.0E+00,
    'patch_icpp(1)%alpha(1)'       : alpha1,
    'patch_icpp(1)%alpha(2)'       : 0.0E+00,
    # ==========================================================
    # Patch 2: Bubble  ======================================
    'patch_icpp(2)%geometry'       : 8,#14,
    'patch_icpp(2)%x_centroid'     : x_centroidb,
    'patch_icpp(2)%y_centroid'     : y_centroidb,
    'patch_icpp(2)%z_centroid'     : z_centroidb, 
    'patch_icpp(2)%radius'         : Rin,
    'patch_icpp(2)%smoothen'       : 'T',
    'patch_icpp(2)%smooth_patch_id': 1,
    'patch_icpp(2)%smooth_coeff'   : 0.5E+00,
#    'patch_icpp(2)%non_axis_sym'   : 'F',
#    'patch_icpp(2)%a2'             : 0.0E+00,
#    'patch_icpp(2)%a3'             : 0.0E+00,
#    'patch_icpp(2)%a4'             : 0.0E+00,
#    'patch_icpp(2)%a5'             : 0.0E+00,
#    'patch_icpp(2)%a6'             : 0.0E+00,
#    'patch_icpp(2)%a7'             : 0.0E+00,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%vel(1)'         : 0.E+00,
    'patch_icpp(2)%vel(2)'         : 0.E+00,
    'patch_icpp(2)%vel(3)'         : 0.E+00,
    'patch_icpp(2)%pres'           : Pbn,
    'patch_icpp(2)%alpha_rho(1)'   : 0.0E+00,
    'patch_icpp(2)%alpha_rho(2)'   : alpha2*rhogn,
    'patch_icpp(2)%alpha(1)'       : 0.0E+00,
    'patch_icpp(2)%alpha(2)'       : alpha2,
    # ==========================================================

    # Fluids Physical Parameters ===============================
    'fluid_pp(1)%gamma'            : 1/(nl-1),
    'fluid_pp(1)%pi_inf'           : nl*Pi_inf_ln/(nl-1),
    'fluid_pp(2)%gamma'            : 1/(nb-1),
    'fluid_pp(2)%pi_inf'           : 0.0E+00,
#    'fluid_pp(1)%Re(1)'            : Rel,
#    'fluid_pp(2)%Re(1)'            : Reb,
    # ==========================================================
}))

# ==============================================================================

