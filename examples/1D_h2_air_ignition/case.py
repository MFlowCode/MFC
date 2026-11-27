#!/usr/bin/env python3
# Constant-volume ignition of stoichiometric H2/air (h2o2.yaml, 1 atm) in a uniform periodic box, checked against a Cantera
# IdealGasReactor by compare_cantera.py. The box is uniform, so there are no fluxes and the run isolates the kinetics, the
# thermodynamics and the reaction integrator.
#
# Reaction integrators (chem_params%reaction_substeps):
#   2: operator-split alpha-QSS sub-stepping, second order in dt (default; stable at the default dt = 1e-8 s)
#   0: reaction source inside the RK stages; exact to the sampling floor but stable only for dt < 2.51/|lambda_max|
#      (about 6e-10 s for the burnt H2/air state), so use e.g. --substeps 0 --dt 5e-10
# Production values of the measured table: T0 = 1000, 1100, 1200 K, tend = 0.8 ms, dt = 1e-8 s.
#
# Reference: Cantera IdealGasReactor, constant volume, adiabatic (computed by compare_cantera.py).
import argparse
import json

import cantera as ct

parser = argparse.ArgumentParser(prog="1D_h2_air_ignition", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--T0", type=float, default=1100.0, help="Initial temperature [K].")
parser.add_argument("--dt", type=float, default=1e-8, help="Time step [s].")
parser.add_argument("--tend", type=float, default=8e-4, help="End time [s].")
parser.add_argument("--nsaves", type=int, default=200, help="Number of saved times.")
parser.add_argument("--substeps", type=int, default=2, help="chem_params%%reaction_substeps (0 = in-stage reactions).")
args = parser.parse_args()

ctfile = "h2o2.yaml"
sol = ct.Solution(ctfile)
sol.TPX = args.T0, ct.one_atm, "H2:2, O2:1, N2:3.76"

NT = int(round(args.tend / args.dt))
NS = max(1, NT // args.nsaves)

case = {
    # Logistics
    "run_time_info": "F",
    # Computational Domain Parameters
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "m": 25,
    "n": 0,
    "p": 0,
    "dt": float(args.dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": NS,
    "parallel_io": "F",
    # Simulation Algorithm Parameters
    "model_eqns": "5eq",
    "num_fluids": 1,
    "num_patches": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1e-16,
    "weno_avg": "F",
    "mapped_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -1,
    "bc_x%end": -1,
    # Formatted Database Files Structure Parameters
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "chem_wrt_T": "T",
    # Patch 1: uniform state at rest
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": sol.P,
    "patch_icpp(1)%alpha(1)": 1,
    "patch_icpp(1)%alpha_rho(1)": sol.density,
    # Fluids Physical Parameters (overridden by the mixture thermodynamics of the mechanism)
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0,
    # Chemistry
    "chemistry": "T",
    "chem_params%diffusion": "F",
    "chem_params%reactions": "T",
    "chem_params%reaction_substeps": args.substeps,
    "cantera_file": ctfile,
}

for i in range(sol.n_species):
    case[f"chem_wrt_Y({i + 1})"] = "T"
    case[f"patch_icpp(1)%Y({i + 1})"] = float(sol.Y[i])

if __name__ == "__main__":
    print(json.dumps(case))
