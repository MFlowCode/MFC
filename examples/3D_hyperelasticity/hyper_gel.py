import math
import json
#need quadrant and collapse near hypoelastic wall

## Define characteristic values for the sim
Ri = 230.4E-6
# temperature
T = (10*101325+1.0E+09)/(1000*(2.35-1)*1816)
# print("T :: ",T)
rhog = 1
Pb = 3550 #Pb = 101325 
Pl = 101325 #Pl = 20*101325
Po = Pl 
Pi_inf_l = 1.0E+09
Pi_inf_b = 0.0E+00
Pi_inf_o = 1.1754E+09

## fluid properties
# liquid 
nl = 2.35E+00
cv_l = 1816
rhol = (Pl+Pi_inf_l)/((nl-1)*cv_l*T)
# object 
obj_rhol = rhol
cv_o = cv_l
rhoo = 1060
no = 2.35
#no = 1.19E+00
muo = 0.060E+00
# gas 
ng = 1.47E+00

# alpha seeding fractions
bub_wl = 1.0E-12
bub_wo = 1.0E-12
bub_wg = 1 - bub_wl - bub_wo
liq_wo = 1.0E-12
liq_wg = 1.0E-12
liq_wl = 1 - liq_wo - liq_wg
obj_wl = 1.0E-12
obj_wg = 1.0E-12
obj_wo = 1 - obj_wl - obj_wg

Gl = 0.
Gg = 0.
Go = 0.57E+03 #1.0933E+04

## mixture values in the liquid
#rhoml = (alpha1-alph_eps)*rhol+alph_eps*rhog
#pi_inf_m = (alpha1-alph_eps)*Pi_inf_l
#nml = (alpha1-alph_eps)*nl+alph_eps*ng
cl = (nl*(Pl+Pi_inf_l)/rhol)**(0.5)
rhoml = liq_wl*rhol + liq_wg*rhog + liq_wo*rhoo 

## Defining Characteristic Values 
Rc = Ri
rhoc = rhol
uc = (Pl/rhol)**(0.5E+00)
ucc = (Pl/rhoml)**(0.5E+00)
# characteristic collapse time, change later to prevent horrific confusion
tc = Ri/uc
Pc = rhoc*uc**(2.0E+00)

## Non-Dimensionalizing values using characteristic values
Rin = Ri/Rc
rholn = rhol/rhoc
rhogn = rhog/rhoc
rhoon = rhoo/rhoc
Pln = Pl/Pc
Pbn = Pb/Pc
Pon = Po/Pc
tcn = tc/tc
Pi_inf_ln = Pi_inf_l/Pc
clc = cl/uc
sod_nd = 2.17
#clmc = clm/ucc
Pi_inf_on = Pi_inf_o/Pc

# Un comment if dimensional
#Rin = Ri
#rholn = rhol
#rhogn = rhog
#rhoon = rhoo
#Pln = Pl
#Pbn = Pb/rhoo
#Pon = Po
#tcn = tc
#Pi_inf_ln = Pi_inf_l
#leng = domain_length
#clc = cl
#sod_nd = sod

## Non-Dimensional Numbers
#Rel = rhol*uc*2*Rin/mulc
#Reb = rhog*uc*2*Rin/mubc
Rel = 0
Reb = 0
Reo = rhoo*uc*2/muo
Ma = uc/cl
Co = Pl/Go
iCo = 1/Co 

## GEOMETRY:: Grid Specifications
lengx = 8.0E+00
lengy = 3.0E+00
lengz = 3.0E+00
CFL = 0.3
PPBR = 16 #92
x_beg = -5.0E+00
x_end = 3.0E+00
y_beg = 0.0E+00
y_end = lengy
z_beg = 0.0E+00
z_end = lengz
Nx = PPBR*lengx
Ny = PPBR*lengy
Nz = PPBR*lengz
delta_x = lengx/Nx
delta_t = CFL*delta_x/clc

#print("NX :: ",Nx,", NY :: ",Ny,", NZ :: ",Nz)

#delta_t_c = 0.185*delta_x/clmc
#print(delta_t_c)

Nt = int(1.7*tcn/delta_t)
# liquid centroid, patch 1
x_centroidl = x_beg/(2.0E+00)
y_centroidl = (y_end+y_beg)/(2.0E+00)
z_centroidl = (z_end+z_beg)/(2.0E+00)
# bubble centroid, patch 2
x_centroidb = -sod_nd
y_centroidb = 0.0
z_centroidb = 0.0
# objective centroid, patch 3
x_centroido = x_end/2
y_centroido = (y_end+y_beg)/(2.0E+00)
z_centroido = (z_end+z_beg)/(2.0E+00)

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================
    'run_time_info'                : 'T',
    'sim_data'                     : 'T',
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
    't_step_save'                  : int(5),#int(Nt/150),
# ==========================================================

    # Simulation Algorithm Parameters ==========================
    'num_patches'                  : 3,
    'model_eqns'                   : 2,
    'hypoelasticity'               : 'F', 
    'hyperelasticity'              : 'T',
    'pre_stress'                   : 'F',
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 3,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'T',
    'weno_Re_flux'                 : 'F',
    'weno_avg'                     : 'F',
    'riemann_solver'               : 1,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -6, #-16,#-2
    'bc_x%end'                     : -6,
    'bc_y%beg'                     : -2,
    'bc_y%end'                     : -6,
    'bc_z%beg'                     : -2,
    'bc_z%end'                     : -6,
    'stretch_x'                    : 'F',
    'stretch_y'                    : 'F',
    'stretch_z'                    : 'F',
    'a_x'                          : 4.0E+00,
    'x_a'                          : -1.5E+00-sod_nd,                   
    'x_b'                          : 2.5E+00,
    #'loops_x'                      : 0, 
    'a_y'                          : 4.0E+00,
    'y_a'                          : -1.5E+00,
    'y_b'                          : 1.5E+00,
    #'loops_y'                      : 0, 
    'a_z'                          : 4.0E+00,
    'z_a'                          : -1.5E+00,
    'z_b'                          : 1.5E+00,
    #'loops_z'                      : 0,
    # ==========================================================

    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    'probe_wrt'                    :'T',
    'fd_order'                     : 1,
    'num_probes'                   : 1,
    'probe(1)%x'  		   : 0., 		       
    'probe(1)%y'		   : 0., 	    	       
    'probe(1)%z'		   : 0., 	    	       
    # ==========================================================

    # Patch 1: Background  ============================
    'patch_icpp(1)%geometry'       : 9,# for 3D
    'patch_icpp(1)%x_centroid'     : x_centroidl, #100*x_centroidl,
    'patch_icpp(1)%y_centroid'     : y_centroidl, #100*y_centroidl,
    'patch_icpp(1)%z_centroid'     : z_centroidl, #100*z_centroidl,
    'patch_icpp(1)%length_x'       : lengx, #200*lengx,
    'patch_icpp(1)%length_y'       : lengy, #200*lengy,
    'patch_icpp(1)%length_z'       : lengz, #200*lengz,
    'patch_icpp(1)%vel(1)'         : 0.E+00,
    'patch_icpp(1)%vel(2)'         : 0.E+00,
    'patch_icpp(1)%vel(3)'         : 0.E+00,
    'patch_icpp(1)%pres'           : Pln,
    'patch_icpp(1)%alpha_rho(1)'   : liq_wl*rholn,
    'patch_icpp(1)%alpha_rho(2)'   : liq_wg*rhogn,
    'patch_icpp(1)%alpha_rho(3)'   : liq_wo*rhoon,
    'patch_icpp(1)%alpha(1)'       : liq_wl,
    'patch_icpp(1)%alpha(2)'       : liq_wg,
    'patch_icpp(1)%alpha(3)'       : liq_wo,
    # ==========================================================
    # Patch 2: Bubble  ======================================
     # Specify the spherical gas bubble grid geometry
    'patch_icpp(2)%geometry'       : 8,# for 3D
    'patch_icpp(2)%smoothen'       : 'T',
    'patch_icpp(2)%smooth_patch_id' : 1,
    'patch_icpp(2)%smooth_coeff'   : 4.0E+00,
    'patch_icpp(2)%x_centroid'     : x_centroidb,
    'patch_icpp(2)%y_centroid'     : y_centroidb,
    'patch_icpp(2)%z_centroid'     : z_centroidb,
    'patch_icpp(2)%radius'         : Rin,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    # Specify the patch primitive variables 
    'patch_icpp(2)%vel(1)'         : 0.E+00,
    'patch_icpp(2)%vel(2)'         : 0.E+00,
    'patch_icpp(2)%vel(3)'         : 0.E+00,
    'patch_icpp(2)%pres'           : Pbn,
    'patch_icpp(2)%alpha_rho(1)'   : bub_wl*rholn,
    'patch_icpp(2)%alpha_rho(2)'   : bub_wg*rhogn,
    'patch_icpp(2)%alpha_rho(3)'   : bub_wo*rhoon,
    'patch_icpp(2)%alpha(1)'       : bub_wl,
    'patch_icpp(2)%alpha(2)'       : bub_wg,
    'patch_icpp(2)%alpha(3)'       : bub_wo,
    # ==========================================================
     # Patch 3: Gel  ===========================================
     # Specify the gel grid geometry
    'patch_icpp(3)%geometry'       : 9,# for 3D
    'patch_icpp(3)%x_centroid'     : x_centroido, #100*x_centroido,
    'patch_icpp(3)%y_centroid'     : y_centroido, #100*y_centroido,
    'patch_icpp(3)%z_centroid'     : z_centroido, #100*z_centroido,
    'patch_icpp(3)%length_x'       : 3.0E+00, #100*lengx,
    'patch_icpp(3)%length_y'       : lengy, #200*lengy,
    'patch_icpp(3)%length_z'       : lengz, #200*lengz,
    'patch_icpp(3)%alter_patch(1)' : 'T',
    # Specify the patch primitive variables 
    'patch_icpp(3)%vel(1)'         : 0.E+00,
    'patch_icpp(3)%vel(2)'         : 0.E+00,
    'patch_icpp(3)%vel(3)'         : 0.E+00,
    'patch_icpp(3)%pres'           : Pon,
    'patch_icpp(3)%alpha_rho(1)'   : obj_wl*rholn,
    'patch_icpp(3)%alpha_rho(2)'   : obj_wg*rhogn,
    'patch_icpp(3)%alpha_rho(3)'   : obj_wo*rhoon,
    'patch_icpp(3)%alpha(1)'       : obj_wl,
    'patch_icpp(3)%alpha(2)'       : obj_wg,
    'patch_icpp(3)%alpha(3)'       : obj_wo,
    # ==========================================================


    # Fluids Physical Parameters ===============================
    'fluid_pp(1)%gamma'            : 1./(nl-1.),
    'fluid_pp(1)%pi_inf'           : nl*Pi_inf_ln/(nl-1.),
    'fluid_pp(1)%G'                : Gl,
    'fluid_pp(2)%gamma'            : 1./(ng-1.),
    'fluid_pp(2)%pi_inf'           : 0.0E+00,
    'fluid_pp(2)%G'                : Gg,
    'fluid_pp(3)%gamma'            : 1./(no-1.),
    'fluid_pp(3)%pi_inf'           : no*Pi_inf_on/(no-1.),
    'fluid_pp(3)%G'                : iCo,
    'fluid_pp(1)%qv'               : 0.0E+00,
    'fluid_pp(1)%qvp'              : 0.0E+00,
    'fluid_pp(2)%qv'               : 0.0E+00,
    'fluid_pp(2)%qvp'              : 0.0E+00,
    'fluid_pp(3)%qv'               : 0.0E+00,
    'fluid_pp(3)%qvp'              : 0.0E+00,
#    'fluid_pp(1)%Re(1)'            : Rel,
#    'fluid_pp(2)%Re(1)'            : Reb,
#    'fluid_pp(3)%Re(1)'            : Reo,
    # ===========++=============================================
}))

# ==============================================================================


