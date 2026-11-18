#!/usr/bin/env python3
# Benchmark amr
# Additional Benchmarked Features
# - amr : T (one level, ratio 2, dynamic regrid every 4 steps)
# - amr_subcycle : T

import argparse
import json
import math

parser = argparse.ArgumentParser(prog="Benchmarking Case AMR", description="This MFC case was created for the purposes of benchmarking MFC.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--gbpp", type=int, metavar="MEM", default=16, help="Adjusts the problem size per rank to fit into [MEM] GB of GPU memory per GPU.")
parser.add_argument("--steps", type=int, default=None, help="Override t_step_stop/t_step_save.")

ARGS = vars(parser.parse_args())
DICT = ARGS["mfc"]

size = 1 if DICT["gpu"] else 0

# A spherical blast: the refined shell grows with the shock and is regridded every 4 steps, so the step, the fine-coarse
# coupling and the regrid (tagging, clustering, migration) are all timed. The base grid gets a quarter of the cells a
# uniform-grid case would: the fine blocks, their halos and the store take the rest.
ppg = 8000000 / 16.0
procs = DICT["nodes"] * DICT["tasks_per_node"]
ncells = math.floor(ppg * procs * ARGS["gbpp"] / 4.0)
N = max(32, 4 * math.floor(ncells ** (1 / 3) / 4))

dx = 1.0 / N
steps = ARGS["steps"] if ARGS["steps"] is not None else int(4 * (5 * size + 5))  # twice the others: a regrid every 4 steps

# initial block: the cube around the blast, in base-grid cells
lo, hi = math.floor(0.3 * N), math.ceil(0.7 * N) - 1

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "F",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "z_domain%beg": 0.0,
            "z_domain%end": 1.0,
            "m": N - 1,
            "n": N - 1,
            "p": N - 1,
            "dt": 0.1 * dx,
            "t_step_start": 0,
            "t_step_stop": steps,
            "t_step_save": steps,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 1,
            "mpp_lim": "F",
            "time_stepper": 3,
            "recon_type": 1,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            # AMR
            "amr": "T",
            "amr_subcycle": "T",
            "amr_regrid_int": 4,
            "amr_tag_eps": 0.05,
            "amr_buf": 4,
            "amr_max_blocks": 1024,
            "amr_block_beg(1)": lo,
            "amr_block_end(1)": hi,
            "amr_block_beg(2)": lo,
            "amr_block_end(2)": hi,
            "amr_block_beg(3)": lo,
            "amr_block_end(3)": hi,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: ambient gas
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%z_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%length_z": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 0.1,
            "patch_icpp(1)%alpha_rho(1)": 0.125,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: the high-pressure sphere
            "patch_icpp(2)%geometry": 8,
            "patch_icpp(2)%x_centroid": 0.5,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%z_centroid": 0.5,
            "patch_icpp(2)%radius": 0.1,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": 1.0,
            "patch_icpp(2)%alpha_rho(1)": 1.0,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0,
        }
    )
)
