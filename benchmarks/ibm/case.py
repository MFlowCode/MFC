#!/usr/bin/env python3
# Benchmark ibm_T
# Additional Benchmarked Features
# - ibm : T

import json, math, argparse

parser = argparse.ArgumentParser(prog="Benchmarking Case 4", description="This MFC case was created for the purposes of benchmarking MFC.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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

dx = 1.0 / (1.0 * (Nx + 1))

Tend = 1e-4
Nt = 200
mydt = Tend / (1.0 * Nt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "F",
            # Computational Domain Parameters
            "x_domain%beg": 0.0e00,
            "x_domain%end": 1.0e00,
            "y_domain%beg": 0,
            "y_domain%end": 0.5,
            "z_domain%beg": 0.0,
            "z_domain%end": 0.5,
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "dt": mydt,
            "t_step_start": 0,
            "t_step_stop": int(20 * (5 * size + 5)),
            "t_step_save": int(20 * (5 * size + 5)),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "T",
            "weno_avg": "T",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            # Turn on IBM
            "ib": "T",
            "num_ibs": 1,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1 L
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.25,
            "patch_icpp(1)%z_centroid": 0.25,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 0.5,
            "patch_icpp(1)%length_z": 0.5,
            "patch_icpp(1)%vel(1)": 0.1,
            "patch_icpp(1)%vel(2)": 0,
            "patch_icpp(1)%vel(3)": 0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 0.8e00,
            "patch_icpp(1)%alpha(1)": 0.8e00,
            "patch_icpp(1)%alpha_rho(2)": 0.2e00,
            "patch_icpp(1)%alpha(2)": 0.2e00,
            # Patch: Sphere Immersed Boundary
            "patch_ib(1)%geometry": 8,
            "patch_ib(1)%x_centroid": 0.25,
            "patch_ib(1)%y_centroid": 0.25,
            "patch_ib(1)%z_centroid": 0.25,
            "patch_ib(1)%radius": 0.1,
            # Fluids Physical Parameters
            # Specify 2 fluids
            "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(1)%Re(1)": 54000,
            "fluid_pp(2)%gamma": 1.0e00 / (1.4 - 1.0e00),
            "fluid_pp(2)%pi_inf": 0,
            "fluid_pp(2)%Re(1)": 54000,
        }
    )
)
