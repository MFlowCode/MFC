#!/usr/bin/env python3
"""
A fixed carbon cylinder (2D) or sphere (--ndim 3) gasifying in quiescent hot, high-pressure moist air with a CO flame off the surface,
for the burn rate and flame against a 1D radially symmetric transient reference (ref_cylinder.txt, ref_sphere.txt; ref1d.py).

    ./mfc.sh run examples/nD_ibm_carbon_burn/case.py -n 1 -t pre_process simulation -- --cells_per_D 24 --tend 6e-5
    python3 examples/nD_ibm_carbon_burn/compare_burn_rate.py examples/nD_ibm_carbon_burn

Production values: D = 100 um, 10 atm, T_w = T_inf = 2000 K, 12 to 48 cells per D (cylinder) or 12 to 24 (sphere), 60 us. The body appears
at t = 0 in a box of half-width 6 R with extrapolation faces. The surface chemistry is solved with ib_flux_wall: the wall faces carry the
surface solve's species fluxes. Mechanisms and gas are those of examples/2D_ibm_reacting_surface.
"""

import argparse
import json
import os

import cantera as ct
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MECH_DIR = os.path.join(HERE, "..", "2D_ibm_reacting_surface")
GAS = os.path.normpath(os.path.join(MECH_DIR, "carbon_gasphase_reduced_gri11.yaml"))
SURF = os.path.normpath(os.path.join(MECH_DIR, "carbon_surface_bradley_11species.yaml"))
AIR = {"N2": 0.760557, "O2": 0.233, "H2O": 0.0058, "CO2": 0.000643}

parser = argparse.ArgumentParser(prog="nD_ibm_carbon_burn")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--ndim", type=int, default=2, choices=(2, 3), help="2: cylinder, 3: sphere.")
parser.add_argument("--D", type=float, default=100e-6, help="Body diameter [m].")
parser.add_argument("--p_atm", type=float, default=10.0, help="Pressure [atm].")
parser.add_argument("--T", type=float, default=2000.0, help="Gas and wall temperature [K].")
parser.add_argument("--H_over_R", type=float, default=6.0, help="Box half-width over the body radius.")
parser.add_argument("--cells_per_D", type=int, default=12, help="Cells across the diameter.")
parser.add_argument("--tend", type=float, default=6.0e-5, help="End time [s].")
parser.add_argument("--frames", type=int, default=30, help="Saves.")
parser.add_argument("--cfl", type=float, default=0.5, help="Target CFL (adaptive time step).")
parser.add_argument("--set", action="append", default=[], metavar="KEY=VAL", help="Override a case parameter.")
args, _ = parser.parse_known_args()

p = args.p_atm * ct.one_atm
gas = ct.Solution(GAS)
gas.TPY = args.T, p, AIR
R = 0.5 * args.D
H = args.H_over_R * R
dx = args.D / args.cells_per_D
n = int(round(2 * H / dx))  # cells per direction
ND = args.ndim
case = {
    "run_time_info": "T",
    "cyl_coord": "F",
    "cfl_adap_dt": "T",
    "cfl_target": args.cfl,
    "n_start": 0,
    "t_save": args.tend / args.frames,
    "t_stop": args.tend,
    "num_fluids": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": 2,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "viscous": "T",
    "fd_order": 4,
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "T",
    "chem_params%transport_model": 1,
    "chem_params%reaction_substeps": 2,
    "chem_params%adap_substeps": "T",
    "chem_params%reaction_substeps_max": 16,
    "chem_params%transport_dt": "T",
    "cantera_file": GAS,
    "surface_cantera_file": SURF,
    "surface_phase": "carbon_surface",
    "chem_wrt_T": "T",
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / gas.viscosity,
    "parallel_io": "T" if args.mfc.get("mpi", True) else "F",
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "num_patches": 1,
    "patch_icpp(1)%geometry": 3 if ND == 2 else 9,
    "patch_icpp(1)%pres": p,
    "patch_icpp(1)%alpha_rho(1)": gas.density,
    "patch_icpp(1)%alpha(1)": 1.0,
    # the body, held at T_w with the surface reaction and the flux wall
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 2 if ND == 2 else 8,
    "patch_ib(1)%radius": R,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%thermal_bc": 1,
    "patch_ib(1)%Twall": args.T,
    "patch_ib(1)%surface_reaction": 1,
    "ib_flux_wall": "T",
    "ib_state_wrt": "T",
    "ib_heat_wrt": "T",
}
for i, y in enumerate(gas.Y):
    case[f"patch_icpp(1)%Y({i + 1})"] = float(y)
    case[f"chem_wrt_Y({i + 1})"] = "T"
for k, ax in enumerate("xyz"[:ND]):
    case.update(
        {
            f"{ax}_domain%beg": -H,
            f"{ax}_domain%end": H,
            "mnp"[k]: n - 1,
            f"bc_{ax}%beg": -3,
            f"bc_{ax}%end": -3,
            f"patch_icpp(1)%{ax}_centroid": 0.0,
            f"patch_icpp(1)%length_{ax}": 2 * H,
            f"patch_icpp(1)%vel({k + 1})": 0.0,
            f"patch_ib(1)%{ax}_centroid": 0.0,
        }
    )
if ND == 2:
    case["p"] = 0
for kv in args.set:
    k, v = kv.split("=", 1)
    try:
        case[k] = json.loads(v)
    except json.JSONDecodeError:
        case[k] = v

if __name__ == "__main__":
    print(json.dumps(case))
