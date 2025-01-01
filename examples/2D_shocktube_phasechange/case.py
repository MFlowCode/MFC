#!/usr/bin/env python3
import math, json

## 1 FOR BACKGROUND, 2 FOR SHOKED STATE, 3 FOR WATER REGION (WHEN NEEDED)
# Pressure
p01 = 1.0843e05
p02 = 2.1114e05
p03 = 1.0685e05

## FLUID PROPERTIES FOR EACH PATCH

## liquid water

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

# density
rho0wl1 = 1.08599e03
rho0wl2 = 1.02883e03
rho0wl3 = 1.08714e03

# volume fraction
awl1 = 2.2640e-05
awl2 = 5.7742e-05
awl3 = 6.5969e-03

# speed of sound
c_wl1 = math.sqrt(gamwl * (p01 + piwl) / rho0wl1)
c_wl2 = math.sqrt(gamwl * (p02 + piwl) / rho0wl2)

## Vapor water

# pi infinity
piwv = 0
# qv
qvwv = 2030000
# qv'
qvpwv = -23400
# cv
cvwv = 1040
# cp
cpwv = 1487
# gamma
gamwv = cpwv / cvwv

# density
rho0wv1 = 6.45561e-01
rho0wv2 = 1.19086e00
rho0wv3 = 6.36894e-01

# volume fraction
awv1 = 9.3851e-02
awv2 = 2.9270e-02
awv3 = 5.9082e-01

# speed of sound
c_wv1 = math.sqrt(gamwv * (p01 + piwv) / rho0wv1)
c_wv2 = math.sqrt(gamwv * (p02 + piwv) / rho0wv2)

# Air

# pi infinity
pia = 0
# qv
qva = 0e0
# qv'
qvpa = 0e0
# cv
cva = 717.5
# cp
cpa = 1006
# gamma
gama = cpa / cva

# density
rho0a1 = 1.000227e00
rho0a2 = 1.845120e00
rho0a3 = 9.867987e-01

# volume fractions
aa1 = 1 - awl1 - awv1
aa2 = 1 - awl2 - awv2
aa3 = 1 - awl3 - awv3

# Speed of sound
c_a1 = math.sqrt(gama * (p01 + pia) / rho0a1)
c_a2 = math.sqrt(gama * (p02 + pia) / rho0a2)

## SHOCK RELATIONS
p02Op01 = p02 / p01

# Mach number of the shocked region - this should agree with Min, if everything is correct
Ms = math.sqrt((gama + 1.0) / (2.0 * gama) * (p02Op01 - 1.0) * (p02 / (p02 + pia)) + 1.0)

# shock speed
ss = Ms * c_a1

## SIMULATION PARAMETERS

# CFL
cfl = 0.25

# discretization parameters for the x direction
# number of elements
Nx = 250

# domain boundaries
xb = 0.00
xe = 0.25

# typical cell size
dx = (xe - xb) / Nx

# discretization parameters for the y direction
# number of elements
Ny = 50

# domain boundaries
yb = 0.0
ye = 0.1

# typical cell size
dy = (ye - yb) / Ny

# time step
dt = cfl * dx / ss

# save frequency = SF + 1 (because the initial state, 0.dat, is also saved)
SF = 40

# making Nt divisible by SF
tendA = (xe - xb) / ss * 5 / 24

# 1 - ensure NtA is sufficient to go a little beyond tendA
NtA = int(tendA // dt + 1)

# Array of saves. it is the same as Nt/Sf = t_step_save
AS = int(NtA // SF + 1)

# Nt = total number of steps. Ensure Nt > NtA (so the total tendA is covered)
Nt = AS * SF

# Total physical time
tend = Nt * dt

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
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": AS,
            # Simulation Algorithm Parameters
            "num_patches": 3,
            "model_eqns": 3,
            "num_fluids": 3,
            "mpp_lim": "T",
            "mixture_err": "T",
            "relax": "T",
            "relax_model": 6,
            "palpha_eps": 1.0e-2,
            "ptgalpha_eps": 1.0e-2,
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -6,
            "bc_x%end": -6,
            "bc_y%beg": -2,
            "bc_y%end": -6,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1 - Background
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": (xe + xb) * 2 / 4,
            "patch_icpp(1)%y_centroid": (ye + yb) * 2 / 4,
            "patch_icpp(1)%length_x": (xe - xb) * 2 / 2,
            "patch_icpp(1)%length_y": (ye - yb) * 2 / 2,
            "patch_icpp(1)%vel(1)": 0.0e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": p01,
            "patch_icpp(1)%alpha_rho(1)": awl1 * rho0wl1,
            "patch_icpp(1)%alpha_rho(2)": awv1 * rho0wv1,
            "patch_icpp(1)%alpha_rho(3)": aa1 * rho0a1,
            "patch_icpp(1)%alpha(1)": awl1,
            "patch_icpp(1)%alpha(2)": awv1,
            "patch_icpp(1)%alpha(3)": aa1,
            # Patch 2 - Shocked State
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": (xe + xb) * 1 / 4,
            "patch_icpp(2)%y_centroid": (ye + yb) * 2 / 4,
            "patch_icpp(2)%length_x": (xe - xb) * 1 / 2,
            "patch_icpp(2)%length_y": (ye - yb) * 2 / 2,
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
            # patch 3: Water
            "patch_icpp(3)%geometry": 3,
            "patch_icpp(3)%x_centroid": (xe + xb) * 2.5 / 4,
            "patch_icpp(3)%y_centroid": (ye + yb) * 2 / 4,
            "patch_icpp(3)%length_x": (xe - xb) * 1 / 8,
            "patch_icpp(3)%length_y": (ye - yb) * 2 / 2,
            "patch_icpp(3)%vel(1)": 0.0e00,
            "patch_icpp(3)%vel(2)": 0.0e00,
            "patch_icpp(3)%pres": p03,
            "patch_icpp(3)%alpha_rho(1)": awl3 * rho0wl3,
            "patch_icpp(3)%alpha_rho(2)": awv3 * rho0wv3,
            "patch_icpp(3)%alpha_rho(3)": aa3 * rho0a3,
            "patch_icpp(3)%alpha(1)": awl3,
            "patch_icpp(3)%alpha(2)": awv3,
            "patch_icpp(3)%alpha(3)": aa3,
            "patch_icpp(3)%alter_patch(1)": "T",
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
