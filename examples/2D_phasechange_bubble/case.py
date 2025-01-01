#!/usr/bin/env python3
import math, json

## 1 FOR BACKGROUND, 2 FOR BUBBLE
# Pressure [Pa]
p01 = 5e6
p02 = 3550

# Temperature [K]
T01 = 298.15
T02 = 298.15
# T02 = 7.914

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

# speed of sound FOR
c_wl1 = math.sqrt(gamwl * (p01 + piwl) / rho0wl1)
c_wl2 = math.sqrt(gamwl * (p02 + piwl) / rho0wl2)

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

# speed of sound
c_wv1 = math.sqrt(gamwv * (p01 + piwv) / rho0wv1)
c_wv2 = math.sqrt(gamwv * (p02 + piwv) / rho0wv2)

### Air ###

Ra = Ru / (28.966e-3)  # Gas constant for air (Ru/molecular weight) (J/kg-K)
# gamma
gama = 1.4
# cp
cpa = Ra * gama / (gama - 1)
# cv
cva = cpa / gama
# pi infinity
pia = 0.0e0
# qv
qva = 0.0e0
# qv'
qvpa = 0.0e0

## FOR PATCHES 1 & 2 ##

# density
rho0a1 = (p01 + pia) / ((gama - 1) * cva * T01)
rho0a2 = (p02 + pia) / ((gama - 1) * cva * T02)

# Speed of sound
c_a1 = math.sqrt(gama * (p01 + pia) / rho0a1)
c_a2 = math.sqrt(gama * (p02 + pia) / rho0a2)

## SHOCK RELATIONS
p02Op01 = p02 / p01

# Mach number of the shocked region - this should agree with Min, if everything is correct
Ms = math.sqrt((gama + 1.0) / (2.0 * gama) * (p02Op01 - 1.0) * (p02 / (p02 + pia)) + 1.0)

# shock speed
ss = Ms * c_a1

### volume fractions for each of the patches ###
C0 = 0.25  # vapor concentration for IMR

# water liquid
awl1 = 1.00e00 - 2.00e-12
awl2 = 1.00e-12
# water vapor
awv1 = 1.00e-12
awv2 = 1 / ((1 - C0) / C0 * rho0wv2 / rho0a2 + 1)
# air
aa1 = 1.0 - awl1 - awv1
aa2 = 1.0 - awl2 - awv2

## SIMULATION PARAMETERS

# CFL
cfl = 0.50

# Bubble Initial Radius
R0 = 30e-06

# number of elements
Nx0 = 400
Nx = 1600
Ny = 1600
Nz = 1600

# domain boundaries
xb = 0.00
xe = 120e-6

yb = 0.00
ye = 120e-6

zb = 0.00
ze = 120e-6

# typical cell size
dx = (xe - xb) / Nx
dy = (ye - yb) / Ny
dz = (ze - zb) / Nz

# time step

# save frequency = SF + 1 (because the initial state, 0.dat, is also saved)
SF = 200

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
Nt = int(18e3 * tend // tc * Nx / Nx0 + 1)

dt = tend / Nt

AS = int(Nt // SF)

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
            "stretch_x": "T",
            "loops_x": 3,
            "a_x": 4.0e0,
            "x_a": -2.0 * R0,
            "x_b": 2.0 * R0,
            "stretch_y": "T",
            "loops_y": 3,
            "a_y": 4.0e0,
            "y_a": -2.0 * R0,
            "y_b": 2.0 * R0,
            "cyl_coord": "T",
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": AS,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 3,
            "num_fluids": 3,
            "mpp_lim": "T",
            "mixture_err": "T",
            "relax": "T",
            "relax_model": 6,
            "palpha_eps": 1.0e-8,
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
            "bc_x%beg": -2,
            "bc_x%end": -6,
            "bc_y%beg": -2,
            "bc_y%end": -6,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: High pressured water
            # Specify the cubic water background grid geometry
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": (xe + xb) * 500000 / 100,
            "patch_icpp(1)%y_centroid": (ye + yb) * 500000 / 100,
            "patch_icpp(1)%length_x": (xe - xb) * 1000000 / 100,
            "patch_icpp(1)%length_y": (ye - yb) * 1000000 / 100,
            "patch_icpp(1)%vel(1)": 0.0e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": p01,
            "patch_icpp(1)%alpha_rho(1)": awl1 * rho0wl1,
            "patch_icpp(1)%alpha_rho(2)": awv1 * rho0wv1,
            "patch_icpp(1)%alpha_rho(3)": aa1 * rho0a1,
            "patch_icpp(1)%alpha(1)": awl1,
            "patch_icpp(1)%alpha(2)": awv1,
            "patch_icpp(1)%alpha(3)": aa1,
            # Patch 2: (Vapor) Bubble
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%x_centroid": xb,
            "patch_icpp(2)%y_centroid": yb,
            "patch_icpp(2)%radius": R0,
            "patch_icpp(2)%vel(1)": 0.0e00,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%pres": p02,
            "patch_icpp(2)%alpha_rho(1)": awl2 * rho0wl2,
            "patch_icpp(2)%alpha_rho(2)": awv2 * rho0wv2,
            "patch_icpp(2)%alpha_rho(3)": aa2 * rho0a2,
            "patch_icpp(2)%alpha(1)": awl2,
            "patch_icpp(2)%alpha(2)": awv2,
            "patch_icpp(2)%alpha(3)": aa2,
            "patch_icpp(2)%alter_patch(1)": "T",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gamwl - 1),
            "fluid_pp(1)%pi_inf": gamwl * piwl / (gamwl - 1),
            "fluid_pp(1)%cv": cvwl,
            "fluid_pp(1)%qv": qvwl,
            "fluid_pp(1)%qvp": qvpwl,
            "fluid_pp(2)%gamma": 1.0e00 / (gamwv - 1),
            "fluid_pp(2)%pi_inf": gamwv * piwv / (gamwv - 1),
            "fluid_pp(2)%cv": cvwv,
            "fluid_pp(2)%qv": qvwv,
            "fluid_pp(2)%qvp": qvpwv,
            "fluid_pp(3)%gamma": 1.0e00 / (gama - 1),
            "fluid_pp(3)%pi_inf": gama * pia / (gama - 1),
            "fluid_pp(3)%cv": cva,
            "fluid_pp(3)%qv": qva,
            "fluid_pp(3)%qvp": qvpa,
        }
    )
)
