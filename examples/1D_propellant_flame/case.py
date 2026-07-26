#!/usr/bin/env python3
# 1D premixed H2/O2/Ar deflagration: a self-propagating premixed flame, used
# as a stand-in for a solid-propellant burn front. A hot equilibrium-product
# kernel at the left end ignites a flame that eats into the fresh premix.
import argparse
import json

import cantera as ct

parser = argparse.ArgumentParser(prog="1D_propellant_flame")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--length", type=float, default=0.008, help="Domain length [m].")
parser.add_argument("--m", type=int, default=200, help="Number of grid cells.")
parser.add_argument("--kernel", type=float, default=1e-3, help="Ignition-kernel width at the left end [m].")
parser.add_argument("--ar", type=float, default=3.0, help="Ar dilution (H2:2,O2:1,AR:<ar>); higher = calmer/slower flame.")
parser.add_argument("--tend", type=float, default=1.2e-4, help="Physical end time [s].")
parser.add_argument("--cfl", type=float, default=0.35, help="Target acoustic CFL number.")
parser.add_argument("--frames", type=int, default=10, help="Number of saved output frames (raise for a smooth movie).")
args = parser.parse_args()

ctfile = "h2o2.yaml"
X = f"H2:2,O2:1,AR:{args.ar}"
T0, P0 = 300.0, 101325.0

# Fresh premixed reactants (fill the domain).
fresh = ct.Solution(ctfile)
fresh.TPX = T0, P0, X

# Ignition kernel: constant-pressure (isobaric) equilibrium products of the
# same mixture -> a thermodynamically self-consistent hot patch (T, rho, Y
# all from Cantera) that ignites a flame without an initial pressure jump.
burned = ct.Solution(ctfile)
burned.TPX = T0, P0, X
burned.equilibrate("HP")

L = args.length
Nx = args.m
dx = L / Nx

c_max = max(fresh.sound_speed, burned.sound_speed)
dt = args.cfl * dx / c_max

NT = int(args.tend / dt)
NS = max(1, NT // args.frames)

case = {
    "run_time_info": "T",
    # Domain
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "m": Nx,
    "n": 0,
    "p": 0,
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": NS,
    "parallel_io": "T",
    # Algorithm
    "model_eqns": 2,
    "num_fluids": 1,
    "num_patches": 2,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "weno_avg": "F",
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 2,
    "avg_state": 1,
    # BCs: left reflective wall, right outflow
    "bc_x%beg": -2,
    "bc_x%end": -3,
    # Flame is diffusion-controlled: viscous + species/thermal diffusion on
    "viscous": "T",
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "T",
    "chem_params%transport_model": 2,
    "cantera_file": ctfile,
    # Output
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "chem_wrt_T": "T",
    # Patch 1: fresh premix, whole domain
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": L / 2,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": fresh.P,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha_rho(1)": fresh.density,
    # Patch 2: hot equilibrium-product ignition kernel, left end
    "patch_icpp(2)%geometry": 1,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": args.kernel / 2,
    "patch_icpp(2)%length_x": args.kernel,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%pres": burned.P,
    "patch_icpp(2)%alpha(1)": 1.0,
    "patch_icpp(2)%alpha_rho(1)": burned.density,
    # Fluid EOS (ideal-gas closure is bypassed by chemistry, but gamma/Re must be set)
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / fresh.viscosity,
}

# Species mass fractions per patch + OH output (flame marker)
for i in range(fresh.n_species):
    case[f"patch_icpp(1)%Y({i + 1})"] = float(fresh.Y[i])
    case[f"patch_icpp(2)%Y({i + 1})"] = float(burned.Y[i])
case[f"chem_wrt_Y({fresh.species_index('OH') + 1})"] = "T"

if __name__ == "__main__":
    print(json.dumps(case))
