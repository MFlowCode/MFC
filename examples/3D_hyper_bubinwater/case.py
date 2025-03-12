#!/usr/bin/env python3
import math, json

## 1 FOR BACKGROUND, 2 FOR BUBBLE, 3 FOR GEL
# Pressure [Pa]
p01 = 5e6
p02 = 3550
p03 = p01

# Temperature [K]
T01 = 298.15
T02 = 298.15
T03 = T01

#### FLUID PROPERTIES ####

### liquid water ###
# pi infty
piwl = 1.0e09
# qv
qvwl = -1167000
# qv'
qvpwl = 0.0e0
# cv
cvwl = 1816
# cp
cpwl = 4267
# gamma
gamwl = cpwl / cvwl

## FOR PATCHES 1 & 2 ##

# density
rho0wl1 = (p01 + piwl) / ((gamwl - 1) * cvwl * T01)
rho0wl2 = (p02 + piwl) / ((gamwl - 1) * cvwl * T02)
rho0wl3 = (p03 + piwl) / ((gamwl - 1) * cvwl * T03)

# speed of sound FOR
c_wl1 = math.sqrt(gamwl * (p01 + piwl) / rho0wl1)
c_wl2 = math.sqrt(gamwl * (p02 + piwl) / rho0wl2)
c_wl3 = math.sqrt(gamwl * (p03 + piwl) / rho0wl3)

# part for Gases - relations from IMR
Ru = 8.3144598  # Universal gas constant (J/mol-K)

### Vapor water ###
Rv = Ru / (18.01528e-3)  # Gas constant for vapor (Ru/molecular weight) (J/kg-K)
# gamma
gamwv = 1.4
# cp
cpwv = Rv * gamwv / (gamwv - 1)
# cv
cvwv = cpwv / gamwv
# pi infinity
piwv = 0.0e0
# qv
qvwv = 2030000
# qv'
qvpwv = -23400

## FOR PATCHES 1 & 2 ##

# density
rho0wv1 = (p01 + piwv) / ((gamwv - 1) * cvwv * T01)
rho0wv2 = (p02 + piwv) / ((gamwv - 1) * cvwv * T02)
rho0wv3 = (p03 + piwv) / ((gamwv - 1) * cvwv * T03)

# speed of sound
c_wv1 = math.sqrt(gamwv * (p01 + piwv) / rho0wv1)
c_wv2 = math.sqrt(gamwv * (p02 + piwv) / rho0wv2)
c_wv3 = math.sqrt(gamwv * (p03 + piwv) / rho0wv3)

### Air ###

Ra = Ru / (28.966e-3)  # Gas constant for air (Ru/molecular weight) (J/kg-K)
gamwa = 1.4
# cp
cpa = Ra * gamwa / (gamwa - 1)
# cv
cva = cpa / gamwa
# pi infinity
pia = 0.0e0
# qv
qvwa = 0.0e0
# qv'
qvpwa = 0.0e0

## FOR PATCHES 1 & 2 ##

# density
rho0wa1 = (p01 + pia) / ((gamwa - 1) * cva * T01)
rho0wa2 = (p02 + pia) / ((gamwa - 1) * cva * T02)
rho0wa3 = (p03 + pia) / ((gamwa - 1) * cva * T03)

# Speed of sound
c_a1 = math.sqrt(gamwa * (p01 + pia) / rho0wa1)
c_a2 = math.sqrt(gamwa * (p02 + pia) / rho0wa2)
c_a3 = math.sqrt(gamwa * (p03 + pia) / rho0wa3)

### 3% polyacrylamide gel ###

# gamma
gamwg = gamwl  # 2.35
# pi infty
pig = piwl  # 1.1754E+09
# qv
qvwg = qvwl  # 0.0E0
# qv'
qvpwg = qvpwl  # 0.0E0
# cv
cvg = cvwl
# cp
cpg = gamwg * cvg

## FOR PATCHES 1 & 2 & 3 ##

# density
rho0wg1 = (p01 + pig) / ((gamwg - 1) * cvg * T01)
rho0wg2 = (p02 + pig) / ((gamwg - 1) * cvg * T02)
rho0wg3 = (p03 + pig) / ((gamwg - 1) * cvg * T03)

# Speed of sound
c_g1 = math.sqrt(gamwg * (p01 + pig) / rho0wg1)
c_g2 = math.sqrt(gamwg * (p02 + pig) / rho0wg2)
c_g3 = math.sqrt(gamwg * (p03 + pig) / rho0wg3)

## SHOCK RELATIONS
p02Op01 = p02 / p01

# Mach number of the shocked region - this should agree with Min, if everything is correct
Ms = math.sqrt((gamwa + 1.0) / (2.0 * gamwa) * (p02Op01 - 1.0) * (p02 / (p02 + pia)) + 1.0)

# shock speed
ss = Ms * c_a1

### volume fractions for each of the patches ###
C0 = 0.25  # vapor concentration for IMR

# patch 1: liquid water
liq_wg = 0
liq_wa = 1.00e-15
liq_wv = 1.00e-15
liq_wl = 1.00e00 - liq_wv - liq_wa - liq_wg
# water vapor
vap_wl = 1.00e-15
vap_wv = 1 / ((1 - C0) / C0 * rho0wv2 / rho0wa2 + 1)
vap_wa = 1.00e-15
vap_wg = 0
vap_tot = vap_wl + vap_wv + vap_wa + vap_wg
# air
air_wl = 1.00e-15
air_wv = vap_tot
air_wg = 0
air_wa = 1.00e00 - air_wl - air_wv - air_wg
# bubble
bub_wl = 1e-15
bub_wv = vap_tot
bub_wg = 0
bub_wa = 1 - bub_wl - bub_wv - bub_wg
# gel
gel_wl = 0
gel_wv = 0
gel_wa = 0
gel_wg = 1.00e00 - gel_wl - gel_wv - gel_wa

## Elasticity
Gl = 0
Gv = 0
Ga = 0
Gg = 0.57e03

## SIMULATION PARAMETERS

# CFL
cfl = 0.50

# Bubble Initial Radius
R0 = 230.4e-06

# number of elements
Nx = 249  # 404 #249
Ny = 124  # 179 #124
Nz = 124  # 179 #124
Nx0 = Nx

# domain boundaries
lref = 4 * R0
xb = -5 * R0
xe = lref

yb = 0.00
ye = lref

zb = 0.00
ze = lref

lenx = xe - xb
leny = ye - yb
lenz = ze - zb

xcenl = 0.0
ycenl = leny / 2.0
zcenl = lenz / 2.0

# xdist = 6.51E-10 #2.17E-5
# sod = xdist/R0
sod = 2.17
xcenb = -sod * R0  # -sod #-sod*R0
ycenb = 0.00
zcenb = 0.00

xceng = xe / 2.0
yceng = ycenl
zceng = zcenl

# typical cell size
dx = (xe - xb) / Nx
dy = (ye - yb) / Ny
dz = (ze - zb) / Nz
# print(dx)
# time step

# save frequency = SF + 1 (because the initial state, 0.dat, is also saved)
SF = 100

# Critical time-step
tc = 0.915 * R0 * math.sqrt(rho0wl1 / p01)

# making Nt divisible by SF
# tendA = 1.5 * tc
tend = 1.2 * tc

# 1 - ensure NtA is sufficient to go a little beyond tendA
# NtA = int( tendA // dt + 1 )

# Array of saves. it is the same as Nt/Sf = t_step_save
# AS = int( NtA // SF + 1 )

# Nt = total number of steps. Ensure Nt > NtA (so the total tendA is covered)
# Nt = AS * SF
Nt = int(2e3 * tend // tc * Nx / Nx0 + 1)
# print(Nt)
dt = tend / Nt

AS = int(Nt // SF)
tstart = 0  # 2304
# Total physical time
# tend = Nt * dt

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": xb,
            "x_domain%end": xe,
            "y_domain%beg": yb,
            "y_domain%end": ye,
            "z_domain%beg": zb,
            "z_domain%end": ze,
            "stretch_x": "T",
            "loops_x": 1,
            "a_x": 4.0e0,
            "x_a": -1.75 * R0 * (abs(sod) + 1),
            "x_b": 5 * R0,
            "stretch_y": "T",
            "loops_y": 1,
            "a_y": 1.0e0,
            "y_a": -1.5 * R0 * abs(sod),
            "y_b": 1.5 * R0 * abs(sod),
            "stretch_z": "T",
            "loops_z": 1,
            "a_z": 1.0e0,
            "z_a": -1.5 * R0 * abs(sod),
            "z_b": 1.5 * R0 * abs(sod),
            "cyl_coord": "F",
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "dt": dt,
            "t_step_start": tstart,
            "t_step_stop": Nt,
            "t_step_save": AS,
            # Simulation Algorithm Parameters
            "num_patches": 3,
            "model_eqns": 3,
            "num_fluids": 4,
            "hypoelasticity": "F",
            "hyperelasticity": "T",
            "mpp_lim": "T",
            "mixture_err": "T",
            "relax": "T",
            "relax_model": 6,
            "palpha_eps": 1.0e-6,
            "ptgalpha_eps": 1.0e-2,
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -6,  # -2,
            "bc_x%end": -6,
            "bc_y%beg": -2,
            "bc_y%end": -6,
            "bc_z%beg": -2,
            "bc_z%end": -6,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "probe_wrt": "T",
            "fd_order": 1,
            "num_probes": 1,
            "probe(1)%x": 0.0,
            "probe(1)%y": 0.0,
            "probe(1)%z": 0.0,
            # Patch 1: High pressured water
            # Specify the cubic water background grid geometry
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 20 * xcenl,
            "patch_icpp(1)%y_centroid": 20 * ycenl,
            "patch_icpp(1)%z_centroid": 20 * zcenl,
            "patch_icpp(1)%length_x": 20 * lenx,
            "patch_icpp(1)%length_y": 20 * leny,
            "patch_icpp(1)%length_z": 20 * lenz,
            "patch_icpp(1)%vel(1)": 0.0e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%vel(3)": 0.0e00,
            "patch_icpp(1)%pres": p01,
            "patch_icpp(1)%alpha_rho(1)": liq_wl * rho0wl1,
            "patch_icpp(1)%alpha_rho(2)": liq_wv * rho0wv1,
            "patch_icpp(1)%alpha_rho(3)": liq_wa * rho0wa1,
            "patch_icpp(1)%alpha_rho(4)": liq_wg * rho0wg1,
            "patch_icpp(1)%alpha(1)": liq_wl,
            "patch_icpp(1)%alpha(2)": liq_wv,
            "patch_icpp(1)%alpha(3)": liq_wa,
            "patch_icpp(1)%alpha(4)": liq_wg,
            # Patch 2: (Vapor) Bubble
            "patch_icpp(2)%geometry": 8,
            "patch_icpp(2)%x_centroid": xcenb,
            "patch_icpp(2)%y_centroid": ycenb,
            "patch_icpp(2)%z_centroid": zcenb,
            "patch_icpp(2)%radius": R0,
            "patch_icpp(2)%vel(1)": 0.0e00,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%vel(3)": 0.0e00,
            "patch_icpp(2)%pres": p02,
            "patch_icpp(2)%alpha_rho(1)": bub_wl * rho0wl2,
            "patch_icpp(2)%alpha_rho(2)": bub_wv * rho0wv2,
            "patch_icpp(2)%alpha_rho(3)": bub_wa * rho0wa2,
            "patch_icpp(2)%alpha_rho(4)": bub_wg * rho0wg2,
            "patch_icpp(2)%alpha(1)": bub_wl,
            "patch_icpp(2)%alpha(2)": bub_wv,
            "patch_icpp(2)%alpha(3)": bub_wa,
            "patch_icpp(2)%alpha(4)": bub_wg,
            "patch_icpp(2)%alter_patch(1)": "T",
            # Patch 3: Gel Object
            "patch_icpp(3)%geometry": 9,
            "patch_icpp(3)%x_centroid": 20 * xceng,
            "patch_icpp(3)%y_centroid": 20 * yceng,
            "patch_icpp(3)%z_centroid": 20 * zceng,
            "patch_icpp(3)%length_x": 20 * xe,
            "patch_icpp(3)%length_y": 20 * leny,
            "patch_icpp(3)%length_z": 20 * lenz,
            "patch_icpp(3)%vel(1)": 0.0e00,
            "patch_icpp(3)%vel(2)": 0.0e00,
            "patch_icpp(3)%vel(3)": 0.0e00,
            "patch_icpp(3)%pres": p03,
            "patch_icpp(3)%alpha_rho(1)": gel_wl * rho0wl3,
            "patch_icpp(3)%alpha_rho(2)": gel_wv * rho0wv3,
            "patch_icpp(3)%alpha_rho(3)": gel_wa * rho0wa3,
            "patch_icpp(3)%alpha_rho(4)": gel_wg * rho0wg3,
            "patch_icpp(3)%alpha(1)": gel_wl,
            "patch_icpp(3)%alpha(2)": gel_wv,
            "patch_icpp(3)%alpha(3)": gel_wa,
            "patch_icpp(3)%alpha(4)": gel_wg,
            "patch_icpp(3)%alter_patch(1)": "T",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gamwl - 1),
            "fluid_pp(1)%pi_inf": gamwl * piwl / (gamwl - 1),
            "fluid_pp(1)%cv": cvwl,
            "fluid_pp(1)%qv": qvwl,
            "fluid_pp(1)%qvp": qvpwl,
            "fluid_pp(1)%G": Gl,
            "fluid_pp(2)%gamma": 1.0e00 / (gamwv - 1),
            "fluid_pp(2)%pi_inf": gamwv * piwv / (gamwv - 1),
            "fluid_pp(2)%cv": cvwv,
            "fluid_pp(2)%qv": qvwv,
            "fluid_pp(2)%qvp": qvpwv,
            "fluid_pp(2)%G": Gv,
            "fluid_pp(3)%gamma": 1.0e00 / (gamwa - 1),
            "fluid_pp(3)%pi_inf": gamwa * pia / (gamwa - 1),
            "fluid_pp(3)%cv": cva,
            "fluid_pp(3)%qv": qvwa,
            "fluid_pp(3)%qvp": qvpwa,
            "fluid_pp(3)%G": Ga,
            "fluid_pp(4)%gamma": 1.0e00 / (gamwg - 1),
            "fluid_pp(4)%pi_inf": gamwg * pig / (gamwg - 1),
            "fluid_pp(4)%cv": cvg,
            "fluid_pp(4)%qv": qvwg,
            "fluid_pp(4)%qvp": qvpwg,
            "fluid_pp(4)%G": Gg,
        }
    )
)
