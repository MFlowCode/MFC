#!/usr/bin/env python3
# References:
# + https://doi.org/10.1016/j.compfluid.2013.10.014: 4.3. Multi-component inert shock tube

import json
import argparse

import cantera as ct

parser = argparse.ArgumentParser(prog="nD_inert_shocktube", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    "--mfc",
    type=json.loads,
    default="{}",
    metavar="DICT",
    help="MFC's toolchain's internal state.",
)
parser.add_argument(
    "--no-chem",
    dest="chemistry",
    default=True,
    action="store_false",
    help="Disable chemistry.",
)

args = parser.parse_args()

ctfile = "h2o2.yaml"
sol_L = ct.Solution(ctfile)
sol_L.TPX = 400, 8000, "H2:2,O2:1,AR:7"
sol_R = ct.Solution(ctfile)
sol_R.TPX = 1200, 80000, "H2:2,O2:1,AR:7"

L = 0.10
Nx = 400
dx = L / Nx
dt = 20e-8
Tend = 40e-6

NT = int(Tend / dt)
SAVE_COUNT = 200
NS = NT // SAVE_COUNT

case = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    "x_domain%beg": -L / 2,
    "x_domain%end": +L / 2,
    "m": Nx,
    "n": 0,
    "p": 0,
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": NS,
    "parallel_io": "F",
    # Simulation Algorithm Parameters
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
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -2,
    "bc_x%end": -3,
    # Chemistry
    "chemistry": "F" if not args.chemistry else "T",
    "chem_params%diffusion": "F",
    "chem_params%reactions": "T",
    # Formatted Database Files Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": -L / 4,
    "patch_icpp(1)%length_x": L / 2,
    "patch_icpp(1)%vel(1)": 0,
    "patch_icpp(1)%pres": sol_L.P,
    "patch_icpp(1)%alpha(1)": 1,
    "patch_icpp(1)%alpha_rho(1)": sol_L.density,
    "patch_icpp(2)%geometry": 1,
    "patch_icpp(2)%x_centroid": L / 4,
    "patch_icpp(2)%length_x": L / 2,
    "patch_icpp(2)%vel(1)": 0,
    "patch_icpp(2)%pres": sol_R.P,
    "patch_icpp(2)%alpha(1)": 1,
    "patch_icpp(2)%alpha_rho(1)": sol_R.density,
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (1.55e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0,
    # Chemistry
    "cantera_file": ctfile,
}

if args.chemistry:
    for i in range(len(sol_L.Y)):
        case[f"patch_icpp(1)%Y({i+1})"] = sol_L.Y[i]
        case[f"patch_icpp(2)%Y({i+1})"] = sol_R.Y[i]

if __name__ == "__main__":
    print(json.dumps(case))
