#!/usr/bin/env python3
# Benchmark igr_T_viscous_T
# Additional Benchmarked Features
# - igr : T
# - viscous : T
# - igr_order : 5

import json, math, argparse

parser = argparse.ArgumentParser(prog="Benchmarking Case 5", description="This MFC case was created for the purposes of benchmarking MFC.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--gbpp", type=int, metavar="MEM", default=16, help="Adjusts the problem size per rank to fit into [MEM] GB of GPU memory per GPU.")

ARGS = vars(parser.parse_args())
DICT = ARGS["mfc"]

size = 1 if DICT["gpu"] else 0

ppg = 8000000 / 16.0
procs = DICT["nodes"] * DICT["tasks_per_node"]
ncells = math.floor(ppg * procs * ARGS["gbpp"])
s = math.floor((ncells) ** (1 / 3))
Nx, Ny, Nz = s, s, s

Re = 1600
L = 1
P0 = 101325
rho0 = 1
C0 = math.sqrt(1.4 * P0)
V0 = 0.1 * C0
mu = V0 * L / Re

cfl = 0.5
dx = 2 * math.pi * L / (Nx + 1)

dt = cfl * dx / (C0)

tC = L / V0
tEnd = 20 * tC

Nt = int(tEnd / dt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -math.pi * L,
            "x_domain%end": math.pi * L,
            "y_domain%beg": -math.pi * L,
            "y_domain%end": math.pi * L,
            "z_domain%beg": -math.pi * L,
            "z_domain%end": math.pi * L,
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "cyl_coord": "F",
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": int(20 * (5 * size + 5)),
            "t_step_save": int(20 * (5 * size + 5)),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "bc_z%beg": -1,
            "bc_z%end": -1,
            "igr": "T",
            "igr_order": 5,
            "igr_iter_solver": 1,
            "num_igr_iters": 3,
            "num_igr_warm_start_iters": 3,
            "alf_factor": 10,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "omega_wrt(1)": "T",
            "omega_wrt(2)": "T",
            "omega_wrt(3)": "T",
            "qm_wrt": "T",
            "fd_order": 4,
            "parallel_io": "T",
            # Patch 1: Background (AIR - 2)
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0,
            "patch_icpp(1)%y_centroid": 0,
            "patch_icpp(1)%z_centroid": 0,
            "patch_icpp(1)%length_x": 2 * math.pi * L,
            "patch_icpp(1)%length_y": 2 * math.pi * L,
            "patch_icpp(1)%length_z": 2 * math.pi * L,
            "patch_icpp(1)%vel(1)": f"{V0}*sin(x/{L})*cos(y/{L})*sin(z/{L})",
            "patch_icpp(1)%vel(2)": f"-{V0}*cos(x/{L})*sin(y/{L})*sin(z/{L})",
            "patch_icpp(1)%vel(3)": 0,
            "patch_icpp(1)%pres": f"{P0} + ({rho0}*{V0}**2/16)*(cos(2*x/{L}) + cos(2*y/{L}))*(cos(2*z/{L}) + 2)",
            "patch_icpp(1)%alpha_rho(1)": 1,
            "patch_icpp(1)%alpha(1)": 1,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1),
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(1)%Re(1)": 1 / mu,
        }
    )
)
