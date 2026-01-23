#!/usr/bin/env python3
import math
import json

# athmospheric pressure - Pa (used as reference value)
patm = 101325

# Initial Droplet Diameter / Reference length - m
D0 = 2.0e-3

# cavity to droplet ratio
CtD = 0.06

# cavity relative eccentricity (distance between radii)
ecc = 0.564

# initial shock distance from the y axis. Note that the droplet center is located at y = 0. Thus, the distance from the shock to
# the droplet is about D0/8
ISD = 5.0 / 8 * D0

## pre-shock properties - AIR

# pressure - Pa
p0a = patm

# density - kg/m3
rho0a = 1.204

# gamma
gama = 1.40

# pi infinity - Pa
pia = 0

# speed of sound - M/s
c_a = math.sqrt(gama * (p0a + pia) / rho0a)

## Droplet - WATER

# surface tension - N / m
st = 0.00e0

# Delta Pressure - Pa
DP = -st * 4 / D0

# initial pressure inside the droplet - Pa
p0w = p0a - DP

# density - kg/m3
rho0w = 1000

# gama
gamw = 6.12

# pi infty - Pa
piw = 3.43e08

# speed of sound - m/s
c_w = math.sqrt(gamw * (p0w + piw) / rho0w)

# Shock Mach number of interest. Note that the post-shock properties can be defined in terms of either
# Min or psOp0a. Just comment/uncomment appropriately
Min = 2.146

## Pos to pre shock ratios - AIR

# pressure
psOp0a = (Min**2 - 1) * 2 * gama / (gama + 1) + 1
# psOp0a = 4.5

# density
rhosOrho0a = (1 + (gama + 1) / (gama - 1) * psOp0a) / ((gama + 1) / (gama - 1) + psOp0a)

# Mach number of the shocked region - just a checker, as it must return "Min"
Ms = math.sqrt((gama + 1.0) / (2.0 * gama) * (psOp0a - 1.0) * (p0a / (p0a + pia)) + 1.0)

# shock speed of sound - m/s
ss = Ms * c_a

## post-shock - AIR

# pressure - Pa
ps = psOp0a * p0a

# density - kg / m3
rhos = rhosOrho0a * rho0a

# post shock speed of sound - m/s
c_s = math.sqrt(gama * (ps + pia) / rhos)

# velocity at the post shock - m/s
vel = c_a / gama * (psOp0a - 1.0) * p0a / (p0a + pia) / Ms

## Domain boundaries - m

# x direction
xb = -2.4707 * D0
xe = 3.6226 * D0

# y direction
yb = 0 * D0
ye = 1.6358 * D0

# Stretching factor, to make sure the domaing is sufficiently large after the mesh stretch
StF = 4.0

# number of elements into y direction
Ny = 1712

# number of elements into x direction
Nx = Ny * 2

# grid delta x if mesh were uniform in x direction - m. Note that I do not need a measure for dy
dx = (xe - xb) / Nx

# I calculating tend twice; first is an estimate, second is
# the actual value used. This is because I am getting errors in the
# post process part every time I approximate the actual Nt by an integer
# number (think of a smarter way).

# dimensionless time
ttilde = 1.92

# auxiliary simulation physical time - s. This is not YET the total simulation time, as it will be corrected so as to avoid
# mismatches in simulation and post_process parts. Note that I wrote it this way so I have better control over the # of autosaves
tendA = ttilde * D0 / vel

# "CFL" number that I use to control both temporal and spatial discretizations, such that the ratio dx/dt remains constant for a given
# simulation
cfl = 0.05

# time-step - s
dt = cfl * dx / ss

# Save Frequency. Note that the number of autosaves will be SF + 1, as th IC (0.dat) is also saved
SF = 400

## making Nt divisible by SF
# 1 - ensure NtA goes slightly beyond tendA
NtA = int(tendA // dt + 1)

# Array of saves. It is the same as Nt/Sf = t_step_save
AS = int(NtA // SF + 1)

# Nt = total number of steps. Note that Nt >= NtA (so at least tendA is completely simulated)
Nt = AS * SF

# total simulation time - s. Note that tend >= tendA
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
            "stretch_x": "T",
            "a_x": 20,
            "x_a": -1.2 * D0,
            "x_b": 1.2 * D0,
            "stretch_y": "T",
            "a_y": 20,
            "y_a": -0.0 * D0,
            "y_b": 1.2 * D0,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "cyl_coord": "T",
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": AS,
            # Simulation Algorithm Parameters
            "num_patches": 4,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
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
            "alpha_wrt": "T",
            "parallel_io": "T",
            # I will use 1 for WATER properties, and 2 for AIR properties
            # Patch 1: Background (AIR - 2)
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": (xb + xe) / 2 * StF,
            "patch_icpp(1)%y_centroid": (yb + ye) / 2 * StF,
            "patch_icpp(1)%length_x": (xe - xb) * StF,
            "patch_icpp(1)%length_y": (ye - yb) * StF,
            "patch_icpp(1)%vel(1)": 0.0e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": p0a,
            "patch_icpp(1)%alpha_rho(1)": 0.0e00,
            "patch_icpp(1)%alpha_rho(2)": rho0a,
            "patch_icpp(1)%alpha(1)": 0.0e00,
            "patch_icpp(1)%alpha(2)": 1.0e00,
            # Patch 2: Shocked state (AIR - 2)
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -ISD - (xe - xb) / 2 * StF,
            "patch_icpp(2)%y_centroid": (yb + ye) / 2 * StF,
            "patch_icpp(2)%length_x": (xe - xb) * StF,
            "patch_icpp(2)%length_y": (ye - yb) * StF,
            "patch_icpp(2)%vel(1)": vel,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%pres": ps,
            "patch_icpp(2)%alpha_rho(1)": 0.0e00,
            "patch_icpp(2)%alpha_rho(2)": rhos,
            "patch_icpp(2)%alpha(1)": 0.0e00,
            "patch_icpp(2)%alpha(2)": 1.0e00,
            # Patch 3: Droplet (WATER - 1)
            "patch_icpp(3)%geometry": 2,
            "patch_icpp(3)%x_centroid": 0.0e00,
            "patch_icpp(3)%y_centroid": 0.0e00,
            "patch_icpp(3)%radius": D0 / 2,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%vel(1)": 0.0e00,
            "patch_icpp(3)%vel(2)": 0.0e00,
            "patch_icpp(3)%pres": p0w,
            "patch_icpp(3)%alpha_rho(1)": rho0w,
            "patch_icpp(3)%alpha_rho(2)": 0.0e00,
            "patch_icpp(3)%alpha(1)": 1.0e00,
            "patch_icpp(3)%alpha(2)": 0.0e00,
            # Patch 4: Cavity (AIR - 2)
            "patch_icpp(4)%geometry": 2,
            "patch_icpp(4)%x_centroid": ecc * D0 / 2,
            "patch_icpp(4)%y_centroid": 0.0e00,
            "patch_icpp(4)%radius": CtD * D0 / 2,
            "patch_icpp(4)%alter_patch(3)": "T",
            "patch_icpp(4)%vel(1)": 0.0e00,
            "patch_icpp(4)%vel(2)": 0.0e00,
            "patch_icpp(4)%pres": p0a,
            "patch_icpp(4)%alpha_rho(1)": 0.0e00,
            "patch_icpp(4)%alpha_rho(2)": rho0a,
            "patch_icpp(4)%alpha(1)": 0.0e00,
            "patch_icpp(4)%alpha(2)": 1.0e00,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gamw - 1),
            "fluid_pp(1)%pi_inf": gamw * piw / (gamw - 1),
            "fluid_pp(2)%gamma": 1.0e00 / (gama - 1),
            "fluid_pp(2)%pi_inf": gama * pia / (gama - 1),
        }
    )
)
