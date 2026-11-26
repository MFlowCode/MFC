#!/usr/bin/env python3
# Near-vacuum collapse: a disc of gas at 1/1000 of the ambient density and pressure (same temperature) is filled by the
# surrounding gas, which converges on the disc centre and rebounds as a strong shock. Unlimited WENO5 drives face pressures
# negative at the disc edge within a few steps; face_bound = 2 (Zhang-Shu bound on the face states) carries it through the
# collapse. The disc centre is off the grid's symmetry lines on purpose. See docs/documentation/validation_positivity.md.
import argparse
import json

parser = argparse.ArgumentParser(description="2D near-vacuum collapse")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=128, help="Cells per direction (default: 128)")
parser.add_argument("--face-bound", type=int, default=2, help="Face-state positivity bound: 0=off, 1, 2 (default: 2)")
parser.add_argument("--fluids", type=int, choices=[1, 2], default=1, help="1: one ideal gas; 2: the disc holds a second gas")
parser.add_argument("--model", choices=["5eq", "6eq"], default="5eq", help="Multiphase model (default: 5eq)")
parser.add_argument("--t-end", type=float, default=0.12, help="End time (default: 0.12, after the collapse at t ~ 0.07)")
parser.add_argument("--steps", type=int, default=0, help="Run this many steps instead of to --t-end")
parser.add_argument("--saves", type=int, default=6, help="Number of outputs (default: 6)")
args = parser.parse_args()

dx = 1.0 / args.N
dt = 0.25 * dx / 8.0  # CFL 0.25 on |u| + c <= 8: gas escaping into vacuum reaches 2c/(gamma - 1) = 5.9
nt = args.steps or int(round(args.t_end / dt))
eps = 1.0e-6  # trace volume fraction of the absent gas
two = args.fluids == 2

case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "y_domain%beg": 0.0,
    "y_domain%end": 1.0,
    "m": args.N - 1,
    "n": args.N - 1,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": nt,
    "t_step_save": max(1, nt // args.saves),
    "num_patches": 2,
    "model_eqns": args.model,
    "num_fluids": args.fluids,
    "mpp_lim": "T" if two else "F",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    # Patch 1: ambient gas at rest
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%y_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%length_y": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 1.0,
    # Patch 2: near-vacuum disc, off-centre
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0.48,
    "patch_icpp(2)%y_centroid": 0.53,
    "patch_icpp(2)%radius": 0.3,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": 1.0e-3,
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "ideal_gas",
}
if args.face_bound:
    case["face_bound"] = args.face_bound
if two:
    # disc gas: gamma = 5/3, its density is 1/1000 of the ambient
    case.update(
        {
            "patch_icpp(1)%alpha_rho(1)": 1.0 - eps,
            "patch_icpp(1)%alpha_rho(2)": 1.0e-3 * eps,
            "patch_icpp(1)%alpha(1)": 1.0 - eps,
            "patch_icpp(1)%alpha(2)": eps,
            "patch_icpp(2)%alpha_rho(1)": eps,
            "patch_icpp(2)%alpha_rho(2)": 1.0e-3 * (1.0 - eps),
            "patch_icpp(2)%alpha(1)": eps,
            "patch_icpp(2)%alpha(2)": 1.0 - eps,
            "fluid_pp(2)%gamma": 1.0 / (5.0 / 3.0 - 1.0),
            "fluid_pp(2)%eos": "ideal_gas",
        }
    )
else:
    case.update(
        {
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(2)%alpha_rho(1)": 1.0e-3,
            "patch_icpp(2)%alpha(1)": 1.0,
        }
    )

print(json.dumps(case))
