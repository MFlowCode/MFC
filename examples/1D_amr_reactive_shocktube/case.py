#!/usr/bin/env python3
# 1D shock-ignited H2/O2/Ar reacting front with AMR tagged on temperature (see docs/documentation/amr.md).
#
# The reactive shock tube of examples/1D_reactive_shocktube: an incident shock reflects off the
# wall at x = 0, the mixture behind the reflected shock ignites, and a reacting front runs after
# the shock into the tube. A 2:1 refined level follows the temperature jumps (amr_tag_var = 6) by
# dynamic regridding, so blocks ride on the shock and on the front while the rest of the tube stays
# coarse. Blocks stay a few cells off the domain edges, so the wall strip where ignition starts is coarse.
# --no-amr runs the same case on a uniform grid at the fine resolution, for comparison.
import argparse
import json

import cantera as ct

parser = argparse.ArgumentParser(prog="1D_amr_reactive_shocktube", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--no-amr", dest="amr", default=True, action="store_false", help="Uniform grid at the fine resolution.")
args = parser.parse_args()

ctfile = "h2o2.yaml"
sol_L = ct.Solution(ctfile)
sol_L.DPX = 0.072, 7173, "H2:2,O2:1,AR:7"
sol_R = ct.Solution(ctfile)
sol_R.DPX = 0.18075, 35594, "H2:2,O2:1,AR:7"
u_r = -487.34

L = 0.12
N = 400  # level-0 cells
Nf = 2 * N  # fine-level cells
# AMR has no reaction sub-stepping. Without it this case (AMR or not, even at a quarter of this dt) fails near 170 us,
# so it stops at 155 us, with the reacting front 19 mm behind the reflected shock.
dt = (L / Nf) / abs(u_r) * 0.02
Tend = 155e-6
NT = int(Tend / dt)
NS = NT // 20

case = {
    # Logistics
    "run_time_info": "T",
    # Computational domain
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "m": (N if args.amr else Nf) - 1,
    "n": 0,
    "p": 0,
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": NS,
    "parallel_io": "F",
    # Simulation algorithm
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
    "chemistry": "T",
    "chem_params%diffusion": "F",
    "chem_params%reactions": "T",
    "cantera_file": ctfile,
    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "chem_wrt_T": "T",
    # Patches: gas at rest (left) and gas moving toward the wall (right)
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": L / 4,
    "patch_icpp(1)%length_x": L / 2,
    "patch_icpp(1)%vel(1)": 0,
    "patch_icpp(1)%pres": sol_L.P,
    "patch_icpp(1)%alpha(1)": 1,
    "patch_icpp(1)%alpha_rho(1)": sol_L.density,
    "patch_icpp(2)%geometry": 1,
    "patch_icpp(2)%x_centroid": 3 * L / 4,
    "patch_icpp(2)%length_x": L / 2,
    "patch_icpp(2)%vel(1)": u_r,
    "patch_icpp(2)%pres": sol_R.P,
    "patch_icpp(2)%alpha(1)": 1,
    "patch_icpp(2)%alpha_rho(1)": sol_R.density,
    # Fluid
    "fluid_pp(1)%gamma": 1.0 / (4.4 - 1.0),
    "fluid_pp(1)%eos": 1,
    "fluid_pp(1)%pi_inf": 0,
}

for i in range(len(sol_L.Y)):
    case[f"chem_wrt_Y({i + 1})"] = "T"
    case[f"patch_icpp(1)%Y({i + 1})"] = sol_L.Y[i]
    case[f"patch_icpp(2)%Y({i + 1})"] = sol_R.Y[i]

if args.amr:
    case.update(
        {
            # --- AMR --- one 2:1 level on the temperature jumps, regridded every 20 steps
            "amr": "T",
            "amr_block_beg(1)": N // 2 - 8,
            "amr_block_end(1)": N // 2 + 7,
            "amr_regrid_int": 20,
            "amr_tag_var": 6,
            "amr_tag_eps": 0.02,
            "amr_buf": 4,
            "amr_max_grid_size": 64,
        }
    )

if __name__ == "__main__":
    print(json.dumps(case))
