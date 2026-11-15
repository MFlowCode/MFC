#!/usr/bin/env python3
"""
Gas-mass budget of a reacting immersed boundary in a closed (periodic) box.

A hot carbon cylinder (patch_ib%surface_reaction = 1) sits in quiescent air in a fully periodic box, so the gas can
gain or lose mass only through the cylinder's surface. The surface solve gasifies carbon into the gas at the rate
integrated by ib_surface_wrt; mass_budget.py compares the gas the box gains with that gasified mass. The inert species
should not change at all.

    ./mfc.sh run examples/2D_ibm_wall_mass_budget/case.py -n 1 --gpu                  # ghost-cell wall
    ./mfc.sh run examples/2D_ibm_wall_mass_budget/case.py -n 1 --gpu -- --flux_wall   # ib_flux_wall = T
    python3 examples/2D_ibm_wall_mass_budget/mass_budget.py <case directory>

Case parameters are those of examples/2D_ibm_reacting_surface (same mechanisms, same cylinder), in a closed box.
"""

import argparse
import importlib.util
import json
import os

parser = argparse.ArgumentParser(prog="2D_ibm_wall_mass_budget")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--flux_wall", action="store_true", help="Set ib_flux_wall = T (default: the ghost-cell wall).")
parser.add_argument("--cells_per_D", type=int, default=16, help="Cells across the cylinder diameter.")
parser.add_argument("--t_end", type=float, default=3.0e-4, help="End time [s].")
parser.add_argument("--n_save", type=int, default=40, help="Number of saves (the budget is integrated between saves).")
args = parser.parse_args()

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("reacting_surface", os.path.join(here, "..", "2D_ibm_reacting_surface", "case.py"))
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

D = base.D
L = 6.0 * D
m = n = int(round(L / (D / args.cells_per_D))) - 1

case = dict(base.case)
for key in [k for k in case if k.startswith("bc_x%") and k not in ("bc_x%beg", "bc_x%end")]:
    del case[key]
case.update(
    {
        "run_time_info": "F",
        "m": m,
        "n": n,
        "x_domain%end": L,
        "y_domain%end": L,
        "bc_x%beg": -1,
        "bc_x%end": -1,
        "bc_y%beg": -1,
        "bc_y%end": -1,
        "t_stop": args.t_end,
        "t_save": args.t_end / args.n_save,
        "patch_icpp(1)%x_centroid": 0.5 * L,
        "patch_icpp(1)%y_centroid": 0.5 * L,
        "patch_icpp(1)%length_x": L,
        "patch_icpp(1)%length_y": L,
        "patch_icpp(1)%vel(1)": 0.0,
        "patch_ib(1)%x_centroid": 0.5 * L,
        "patch_ib(1)%y_centroid": 0.5 * L + 0.5 * D / args.cells_per_D / 3.0,
        "parallel_io": "T",
        "ib_state_wrt": "F",
        "ib_surface_wrt": "T",
        "chem_wrt_T": "F",
    }
)
for key in [k for k in case if k.startswith("chem_wrt_Y")]:
    case[key] = "F"
if args.flux_wall:
    case["ib_flux_wall"] = "T"

if __name__ == "__main__":
    print(json.dumps(case))
