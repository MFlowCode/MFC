#!/usr/bin/env python3
# 1D reactive shock tube (Fedkiw H2/O2/Ar), driven by the *fully automated*
# stiff-chemistry pipeline:
#
#   * cfl_const_dt  -- the flow timestep is computed once from the hydrodynamic
#                      CFL condition (dt = cfl_target * min dx/(|u|+c)); no hand-
#                      tuned dt. Because the reaction source is operator-split out
#                      of the flow RHS (see below), the CFL sees a clean sound
#                      speed and picks the largest *hydro*-stable step.
#   * adap_substeps -- the reaction integrator (operator-split alpha-QSS) chooses
#                      its sub-step count per flow step, per rank, from the local
#                      chemical stiffness: it sits at reaction_substeps (floor) in
#                      inert/burned gas and ramps toward reaction_substeps_max
#                      (ceiling) only across the ignition front.
#
# The net effect on this case: the flow takes the large hydro-limited step while
# the (cheap) chemistry is sub-cycled only where it is actually stiff -- reaching
# the same detonation structure as a fine-dt reference at a fraction of the cost,
# with no manual dt or sub-step tuning. Without operator splitting, that large a
# step integrated with the reaction in the RK3 RHS goes unstable (NaN).
#
# References:
# + https://doi.org/10.1016/j.ijhydene.2023.03.190:  Verification of numerical method
# + https://doi.org/10.1016/j.compfluid.2013.10.014: 4.7. Multi-species reactive shock tube

import argparse
import json

import cantera as ct

parser = argparse.ArgumentParser(prog="1D_reactive_shocktube_adaptive", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--no-chem", dest="chemistry", default=True, action="store_false", help="Disable chemistry.")
parser.add_argument("--scale", type=float, default=1, help="Grid scale (multiplies the 400-cell base mesh).")
parser.add_argument("--cfl-target", type=float, default=0.5, help="Target hydrodynamic CFL; sets the (constant) flow timestep.")
parser.add_argument("--substeps", type=int, default=2, help="alpha-QSS sub-step floor (used where the chemistry is not stiff).")
parser.add_argument("--substeps-max", type=int, default=16, help="alpha-QSS sub-step ceiling (used across the ignition front).")
args = parser.parse_args()

ctfile = "h2o2.yaml"
sol_L = ct.Solution(ctfile)
sol_L.DPX = 0.072, 7173, "H2:2,O2:1,AR:7"

sol_R = ct.Solution(ctfile)
sol_R.DPX = 0.18075, 35594, "H2:2,O2:1,AR:7"

u_l = 0
u_r = -487.34

L = 0.12
Nx = int(400 * args.scale)
Tend = 230e-6

case = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "m": Nx,
    "n": 0,
    "p": 0,
    # Automated flow timestep: computed once from the CFL target (no fixed dt)
    "cfl_const_dt": "T",
    "cfl_target": args.cfl_target,
    "n_start": 0,
    "t_save": Tend / 100.0,
    "t_stop": Tend,
    "parallel_io": "F",
    # Simulation Algorithm Parameters
    "model_eqns": "5eq",
    "num_fluids": 1,
    "num_patches": 2,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1e-16,
    "weno_avg": "F",
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -2,
    "bc_x%end": -3,
    # Chemistry: operator-split alpha-QSS reaction integrator with per-rank adaptive sub-stepping
    "chemistry": "F" if not args.chemistry else "T",
    "chem_params%diffusion": "F",
    "chem_params%reactions": "T",
    "chem_params%reaction_substeps": args.substeps,
    "chem_params%reaction_substeps_max": args.substeps_max,
    "chem_params%adap_substeps": "T",
    "cantera_file": ctfile,
    # Formatted Database Files Structure Parameters
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "chem_wrt_T": "T",
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": L / 4,
    "patch_icpp(1)%length_x": L / 2,
    "patch_icpp(1)%vel(1)": u_l,
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
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0,
}

if args.chemistry:
    for i in range(len(sol_L.Y)):
        case[f"chem_wrt_Y({i + 1})"] = "T"
        case[f"patch_icpp(1)%Y({i + 1})"] = sol_L.Y[i]
        case[f"patch_icpp(2)%Y({i + 1})"] = sol_R.Y[i]

if __name__ == "__main__":
    print(json.dumps(case))
