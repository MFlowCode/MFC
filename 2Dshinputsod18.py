import math
import json

Ri = 1.0E-01
Rc = Ri
Rin = Ri/Rc
rhol = 1000
rhog = 1
rhoc = rhol
rholn = rhol/rhoc
rhogn = rhog/rhoc
Pb = 3550
Pl = 5.0E+06
delta_P = Pl-Pb
tc = 0.915*Ri*(rhol/delta_P)**(0.5)
mulc = 9E-04
mubc = 1.0E-05
uc = Ri/tc
Pc = rhoc*uc**(2)
Pln = Pl/Pc
Pbn = Pb/Pc
tcn = tc/tc
CFL = 0.2	 
sod = 18/16
leng = 4.5
PpBr = 256 
Nx = PpBr*leng
Ny = Nx 
x_end = leng
y_end = leng
x_beg = 0.0 
y_beg = 0.0
Pi_inf_l = 702.8E+06
Pi_inf_b = 0.0E+00
nl = 1.47
nb = 1.19
bl = 6.61E-04
bb = 0.0E+00
cl = (nl*(Pl+Pi_inf_l)/(rhol*(1-rhol*bl)))**(0.5)
clc = cl/uc
delta_x = leng/Nx
delta_t = CFL*delta_x/clc
Nt = int(1.255555*tcn/delta_t)
alpha1 = 1.0
alpha2 = 1.0
x_centroidl = (x_end-x_beg)/2.0
y_centroidl = (y_end-y_beg)/2.0
x_centroidb = sod
y_centroidb = 0.0E+00
Rel = rhol*uc*2*Rin/mulc
Reb = rhog*uc*2*Rin/mubc


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
    'm'                            : int(Nx),
    'n'                            : int(Ny),
    'p'                            : 0,
    'cyl_coord'                    : 'T',
    'dt'                           : delta_t,
    't_step_start'                 : 0, 
    't_step_stop'                  : Nt,
    't_step_save'                  : int(Nt/360),
# ==========================================================

    # Simulation Algorithm Parameters ==========================
    'num_patches'                  : 2,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 2,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'T',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'weno_Re_flux'                 : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -5,
    'bc_x%end'                     : -6,
    'bc_y%beg'                     : -2,
    'bc_y%end'                     : -6,
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
    # ==========================================================

    # Patch 1: Background  ============================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : x_centroidl,
    'patch_icpp(1)%y_centroid'     : y_centroidl,
    'patch_icpp(1)%length_x'       : leng,
    'patch_icpp(1)%length_y'       : leng,
    'patch_icpp(1)%vel(1)'         : 0.E+00,
    'patch_icpp(1)%vel(2)'         : 0.E+00,
    'patch_icpp(1)%pres'           : Pln,
    'patch_icpp(1)%alpha_rho(1)'   : alpha1*rholn,
    'patch_icpp(1)%alpha_rho(2)'   : 0.0E+00,
    'patch_icpp(1)%alpha(1)'       : alpha1,
    'patch_icpp(1)%alpha(2)'       : 0.0E+00,
    # ==========================================================
    # Patch 2: Bubble  ======================================
    'patch_icpp(2)%geometry'       : 14,
    'patch_icpp(2)%x_centroid'     : x_centroidb,
    'patch_icpp(2)%y_centroid'     : y_centroidb,
    'patch_icpp(2)%radius'         : Rin,
    'patch_icpp(2)%non_axis_sym'   : 'T',
    'patch_icpp(2)%a2'             : 0.0,
    'patch_icpp(2)%a3'             : 0.0,
    'patch_icpp(2)%a4'             : 0.0,
    'patch_icpp(2)%a5'             : 0.0,
    'patch_icpp(2)%a6'             : 0.0,
    'patch_icpp(2)%a7'             : 0.0,
    'patch_icpp(2)%a8'             : 0.0,
    'patch_icpp(2)%a9'             : 0.0,
    'patch_icpp(2)%a10'            : 0.0,
    'patch_icpp(2)%a11'            : 0.0,
    'patch_icpp(2)%a12'            : 0.0,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%vel(1)'         : 0.E+00,
    'patch_icpp(2)%vel(2)'         : 0.E+00,
    'patch_icpp(2)%pres'           : Pbn,
    'patch_icpp(2)%alpha_rho(1)'   : 0.0E+00,
    'patch_icpp(2)%alpha_rho(2)'   : alpha2*rhogn,
    'patch_icpp(2)%alpha(1)'       : 0.0E+00,
    'patch_icpp(2)%alpha(2)'       : alpha2,
    # ==========================================================

    # Fluids Physical Parameters ===============================
    'fluid_pp(1)%gamma'            : nl,
    'fluid_pp(1)%pi_inf'           : Pi_inf_l/Pc,
    'fluid_pp(2)%gamma'            : nb,
    'fluid_pp(2)%pi_inf'           : Pi_inf_b/Pc,
    'fluid_pp(1)%Re(1)'            : Rel,
    'fluid_pp(2)%Re(1)'            : Reb,
    # ==========================================================
}))

# ==============================================================================

