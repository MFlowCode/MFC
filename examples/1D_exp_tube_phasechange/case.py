#!/usr/bin/env python3
import math, json

# Pressure
p01 = 1.175435854855077e05

# AIR
# density
rho0a1 = 1.077440586592175

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
# Speed of sound
c_a = math.sqrt(gama * (p01 + pia) / rho0a1)

# liquid water
# density
rho0wl1 = 1.079065851351639e03

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
# speed of sound
c_wl = math.sqrt(gamwl * (p01 + piwl) / rho0wl1)

# Vapor water
# density
rho0wv1 = 0.695395098952667

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
# speed of sound
c_wv = math.sqrt(gamwv * (p01 + piwv) / rho0wv1)

# Mach number of the shocked region - this should agree with Min, if everything is correct
Ms = 2.0

ss = Ms * c_a

vel = 2.0e0

# domain boundaries
xb = 0.0
xe = 1.0

# CFL
cfl = 0.40

# Number of elements into the x direction
Nx = 200

dx = (xe - xb) / Nx
dt = cfl * dx / ss

# save frequency = SF + 1 (because the initial state, 0.dat, is also saved)
SF = 200

# making Nt divisible by SF
tendA = (xe - xb) / ss * 0.25

# 1 - ensure NtA is sufficient to go a little beyond tendA
NtA = int(tendA // dt + 1)

# Array of saves. it is the same as Nt/Sf = t_step_save
AS = int(NtA // SF + 1)

# Nt = total number of steps. Ensure Nt > NtA (so the total tendA is covered)
Nt = AS * SF

tend = Nt * dt

# patch properties: aFN in which F is the Fluid and N is the patch number
awl1 = 0.849964400537437

awv1 = 0.029353906292532

aa1 = 1 - awl1 - awv1

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": xb,
            "x_domain%end": xe,
            "m": Nx,
            "n": 0,
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
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1 - Shocked Stated
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": (xe + xb) * 1 / 4,
            "patch_icpp(1)%length_x": (xe - xb) * 1 / 2,
            "patch_icpp(1)%vel(1)": -vel,
            "patch_icpp(1)%pres": p01,
            "patch_icpp(1)%alpha_rho(1)": awl1 * rho0wl1,
            "patch_icpp(1)%alpha_rho(2)": awv1 * rho0wv1,
            "patch_icpp(1)%alpha_rho(3)": aa1 * rho0a1,
            "patch_icpp(1)%alpha(1)": awl1,
            "patch_icpp(1)%alpha(2)": awv1,
            "patch_icpp(1)%alpha(3)": aa1,
            # Patch 2 R
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": (xe + xb) * 3 / 4,
            "patch_icpp(2)%length_x": (xe - xb) * 1 / 2,
            "patch_icpp(2)%vel(1)": vel,
            "patch_icpp(2)%pres": p01,
            "patch_icpp(2)%alpha_rho(1)": awl1 * rho0wl1,
            "patch_icpp(2)%alpha_rho(2)": awv1 * rho0wv1,
            "patch_icpp(2)%alpha_rho(3)": aa1 * rho0a1,
            "patch_icpp(2)%alpha(1)": awl1,
            "patch_icpp(2)%alpha(2)": awv1,
            "patch_icpp(2)%alpha(3)": aa1,
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
