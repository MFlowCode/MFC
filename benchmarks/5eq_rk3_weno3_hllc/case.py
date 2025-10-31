#!/usr/bin/env python3
# Benchmark model_equations_2_time_stepper_3_weno_order_3_riemann_solver_2
# Additional Benchmarked Features
# - model_equations : 2
# - time_stepper : 3
# - weno_order : 3
# - riemann_solver : 2

import json, math, argparse

parser = argparse.ArgumentParser(prog="Benchmarking Case 1", description="This MFC case was created for the purposes of benchmarking MFC.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--gbpp", type=int, metavar="MEM", default=16, help="Adjusts the problem size per rank to fit into [MEM] GB of GPU memory per GPU.")

ARGS = vars(parser.parse_args())
DICT = ARGS["mfc"]

size = 1 if DICT["gpu"] else 0

ppg = 8000000 / 16.0
procs = DICT["nodes"] * DICT["tasks_per_node"]
ncells = math.floor(ppg * procs * ARGS["gbpp"])
s = math.floor((ncells / 2.0) ** (1 / 3))
Nx, Ny, Nz = 2 * s, s, s

# athmospheric pressure - Pa (used as reference value)
patm = 101325

# Initial Droplet Diameter / Reference length - m
D0 = 1.0e-3

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
Min = 2.4

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
xb = -8.4707 * D0
xe = 9.6226 * D0

# xb = -10 * D0
# xe = 10 * D0

# y direction
yb = 0 * D0
ye = 10 * D0

# y direction
zb = 0 * D0
ze = 10 * D0

# Stretching factor, to make sure the domaing is sufficiently large after the mesh stretch
StF = 4.0

# grid delta x if mesh were uniform in x direction - m. Note that I do not need a measure for dy
dx = (xe - xb) / Nx

# I calculate tend twice; first is an estimate, second is
# the actual value used. This is because I am getting errors in the
# post process part every time I approximate the actual Nt by an integer
# number (think of a smarter way).

# dimensionless time
ttilde = 1.92

# auxiliary simulation physical time - s. This is not YET the total simulation time, as it will be corrected so as to avoid
# mismatches in simulation and post_process parts. Note that I wrote it this way so I have better control over the # of autosaves
tendA = ttilde * D0 / vel

cfl = 0.1

# time-step - s
dt = dx * cfl / ss

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
            "run_time_info": "F",
            # Computational Domain Parameters
            "x_domain%beg": xb,
            "x_domain%end": xe,
            "y_domain%beg": yb,
            "y_domain%end": ye,
            "z_domain%beg": zb,
            "z_domain%end": ze,
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "cyl_coord": "F",
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": int(20 * (5 * size + 5)),
            "t_step_save": int(20 * (5 * size + 5)),
            # Simulation Algorithm Parameters
            "num_patches": 3,
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
            "bc_y%end": -3,
            "bc_z%beg": -2,
            "bc_z%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # I will use 1 for WATER properties, and 2 for AIR properties
            # Patch 1: Background (AIR - 2)
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": (xb + xe) / 2 * StF,
            "patch_icpp(1)%y_centroid": (yb + ye) / 2 * StF,
            "patch_icpp(1)%z_centroid": (yb + ye) / 2 * StF,
            "patch_icpp(1)%length_x": (xe - xb) * StF,
            "patch_icpp(1)%length_y": (ye - yb) * StF,
            "patch_icpp(1)%length_z": (ze - zb) * StF,
            "patch_icpp(1)%vel(1)": 0.0e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%vel(3)": 0.0e00,
            "patch_icpp(1)%pres": p0a,
            "patch_icpp(1)%alpha_rho(1)": 0.0e00,
            "patch_icpp(1)%alpha_rho(2)": rho0a,
            "patch_icpp(1)%alpha(1)": 0.0e00,
            "patch_icpp(1)%alpha(2)": 1.0e00,
            # Patch 2: Shocked state (AIR - 2)
            "patch_icpp(2)%geometry": 9,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -ISD - (xe - xb) / 2 * StF,
            "patch_icpp(2)%y_centroid": (yb + ye) / 2 * StF,
            "patch_icpp(2)%z_centroid": (zb + ze) / 2 * StF,
            "patch_icpp(2)%length_x": (xe - xb) * StF,
            "patch_icpp(2)%length_y": (ye - yb) * StF,
            "patch_icpp(2)%length_z": (ze - zb) * StF,
            "patch_icpp(2)%vel(1)": vel,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%vel(3)": 0.0e00,
            "patch_icpp(2)%pres": ps,
            "patch_icpp(2)%alpha_rho(1)": 0.0e00,
            "patch_icpp(2)%alpha_rho(2)": rhos,
            "patch_icpp(2)%alpha(1)": 0.0e00,
            "patch_icpp(2)%alpha(2)": 1.0e00,
            # Patch 3: Droplet (WATER - 1)
            "patch_icpp(3)%geometry": 8,
            "patch_icpp(3)%x_centroid": 0.0e00,
            "patch_icpp(3)%y_centroid": 0.0e00,
            "patch_icpp(3)%z_centroid": 0.0e00,
            "patch_icpp(3)%radius": D0 / 2,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%vel(1)": 0.0e00,
            "patch_icpp(3)%vel(2)": 0.0e00,
            "patch_icpp(3)%vel(3)": 0.0e00,
            "patch_icpp(3)%pres": p0w,
            "patch_icpp(3)%alpha_rho(1)": rho0w,
            "patch_icpp(3)%alpha_rho(2)": 0.0e00,
            "patch_icpp(3)%alpha(1)": 1.0e00,
            "patch_icpp(3)%alpha(2)": 0.0e00,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gamw - 1),
            "fluid_pp(1)%pi_inf": gamw * piw / (gamw - 1),
            "fluid_pp(2)%gamma": 1.0e00 / (gama - 1),
            "fluid_pp(2)%pi_inf": gama * pia / (gama - 1),
        }
    )
)
