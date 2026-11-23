#!/usr/bin/env python3
"""
Heat a cold isothermal cylinder absorbs from quiescent hot gas (a particle meeting hot post-shock gas): wall conduction column.

A cylinder of radius R = 50 um at Tw = 300 K (immersed boundary, thermal_bc = 1) in frozen air at 1500 K and 1 atm, in a box of
half width --half R (the heat reaches 0.3 R by the end) with zero-gradient faces. The gas starts in a layer of thickness --layer R (a near-step; patch_ib%thermal_layer).
The run lasts --tend R^2/alpha_inf (alpha_inf of the hot gas). --wall_column 1 turns on patch_ib%wall_column, the sub-grid
conduction column that carries the flux of a layer thinner than a cell; 0 is the default ghost-cell wall. D/ib_heat.dat gets the
heat the body takes every step.

    ./mfc.sh run examples/2D_ibm_cold_wall_transient/case.py -t pre_process simulation -n 1 -- --cells_per_R 40 --tend 0.1 --wall_column 1
    python3 examples/2D_ibm_cold_wall_transient/check_heat.py examples/2D_ibm_cold_wall_transient

Defaults are a CI-sized run (12 cells per R, half width 3 R, t_end = 0.03 R^2/alpha).
"""

import argparse
import json

import cantera as ct

parser = argparse.ArgumentParser(prog="2D_ibm_cold_wall_transient", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--cells_per_R", type=int, default=12, help="Cells per cylinder radius.")
parser.add_argument("--half", type=float, default=3.0, help="Half width of the box [R].")
parser.add_argument("--tend", type=float, default=0.03, help="End time in R^2/alpha_inf.")
parser.add_argument("--tw", type=float, default=300.0, help="Wall temperature [K].")
parser.add_argument("--tinf", type=float, default=1500.0, help="Gas temperature [K].")
parser.add_argument("--layer", type=float, default=0.004, help="Initial conductive layer thickness [R].")
parser.add_argument("--wall_column", type=int, default=1, choices=(0, 1), help="patch_ib%%wall_column.")
parser.add_argument("--frames", type=int, default=3, help="Saves of the run.")
args = parser.parse_args()

R, P0, half = 50.0e-6, 101325.0, args.half
MECH, X = "h2o2.yaml", "O2:0.21,N2:0.79"
gas = ct.Solution(MECH)
gas.TPX = args.tinf, P0, X
alpha = gas.thermal_conductivity / (gas.density * gas.cp_mass)
t_end = args.tend * R**2 / alpha
dx = R / args.cells_per_R
nh = round(half * args.cells_per_R)

case = {
    "run_time_info": "T",
    "x_domain%beg": -nh * dx,
    "x_domain%end": nh * dx,
    "y_domain%beg": -nh * dx,
    "y_domain%end": nh * dx,
    "m": 2 * nh - 1,
    "n": 2 * nh - 1,
    "p": 0,
    "cfl_adap_dt": "T",
    "cfl_target": 0.5,
    "n_start": 0,
    "t_stop": t_end,
    "t_save": t_end / args.frames,
    "num_patches": 1,
    "num_fluids": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "viscous": "F",
    "fd_order": 2,
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "F",
    "chem_params%transport_model": 1,
    "chem_params%transport_dt": "T",
    "cantera_file": MECH,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": 2 * nh * dx,
    "patch_icpp(1)%length_y": 2 * nh * dx,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": P0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha_rho(1)": gas.density,
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": 0.0,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%radius": R,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%thermal_bc": 1,
    "patch_ib(1)%Twall": args.tw,
    "patch_ib(1)%surface_reaction": 0,
    "patch_ib(1)%thermal_layer": args.layer * R,
    "patch_ib(1)%wall_column": "T" if args.wall_column else "F",
    "ib_heat_wrt": "T",
    "ib_force_stride": 1,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "F",
    "parallel_io": "F",
}
for i, y in enumerate(gas.Y):
    case[f"patch_icpp(1)%Y({i + 1})"] = float(y)

if __name__ == "__main__":
    print(json.dumps(case))
