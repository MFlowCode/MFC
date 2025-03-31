#!/usr/bin/env python3
# Reference:
# + https://doi.org/10.1063/1.1696266
import json, argparse
import cantera as ct

from mfc.case_utils import *

parser = argparse.ArgumentParser(prog="nD_perfect_reactor", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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
parser.add_argument("--scale", type=float, default=1, help="Scale.")
parser.add_argument("--ndim", type=int, default=1, help="Number of dimensions.")

args = parser.parse_args()

ctfile = "h2o2.yaml"
sol = ct.Solution(ctfile)

sol.TPX = 1_600, ct.one_atm, "H2:0.04, O2:0.02, AR:0.94"

Nx = int(25 * args.scale)
Tend = 1e-4
s = 1e-2
dt = 1e-7

NT = int(Tend / dt)
SAVE_COUNT = 20
NS = NT // SAVE_COUNT

case = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    "x_domain%beg": -s / 2,
    "x_domain%end": +s / 2,
    "y_domain%beg": -s / 2,
    "y_domain%end": +s / 2,
    "z_domain%beg": -s / 2,
    "z_domain%end": +s / 2,
    "m": Nx,
    "n": Nx,
    "p": Nx,
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": NS,
    "parallel_io": "T" if args.ndim > 1 and args.mfc.get("mpi", True) else "F",
    # Simulation Algorithm Parameters
    "model_eqns": 2,
    "num_fluids": 1,
    "num_patches": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "weno_avg": "F",
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 1,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "bc_z%beg": -1,
    "bc_z%end": -1,
    # Formatted Database Files Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "chem_wrt_T": "T",
    # Patch 1
    "patch_icpp(1)%geometry": 3 ** (args.ndim - 1),
    "patch_icpp(1)%x_centroid": 0,
    "patch_icpp(1)%y_centroid": 0,
    "patch_icpp(1)%z_centroid": 0,
    "patch_icpp(1)%length_x": s,
    "patch_icpp(1)%length_y": s,
    "patch_icpp(1)%length_z": s,
    "patch_icpp(1)%vel(1)": 0,
    "patch_icpp(1)%vel(2)": 0,
    "patch_icpp(1)%vel(3)": 0,
    "patch_icpp(1)%pres": sol.P,
    "patch_icpp(1)%alpha(1)": 1,
    "patch_icpp(1)%alpha_rho(1)": sol.density,
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0,
}

if args.chemistry:
    case.update(
        {
            # Chemistry
            "chemistry": "T",
            "chem_params%diffusion": "F",
            "chem_params%reactions": "T",
            "cantera_file": ctfile,
        }
    )

    for i in range(len(sol.Y)):
        case[f"chem_wrt_Y({i + 1})"] = "T"
        case[f"patch_icpp(1)%Y({i+1})"] = sol.Y[i]

case = remove_higher_dimensional_keys(case, args.ndim)

if __name__ == "__main__":
    print(json.dumps(case))
