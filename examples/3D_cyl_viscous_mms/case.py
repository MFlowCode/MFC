#!/usr/bin/env python3
# Manufactured-solution check of the cylindrical viscous operator (3D cylindrical and 2D axisymmetric).
# The patch is uniform at rest; run_mms.py overwrites the initial state with the manufactured velocity
# between pre_process and simulation (no analytic patch, so no per-case compile), runs one step with and
# without viscosity, and compares (q_visc - q_inviscid)/dt with the exact rates in exact_rates.py.
# Coordinates: x axial, y = r (axis at 0), z = theta in radians.
import argparse
import json
import math

parser = argparse.ArgumentParser(prog="3D_cyl_viscous_mms", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--case", default="A3", choices=["A3", "B3", "C3", "E2"], help="manufactured field; trailing digit = dimension")
parser.add_argument("-N", type=int, default=32, help="cells in r (and theta in 3D)")
parser.add_argument("--viscous", default="T", choices=["T", "F"])
args = parser.parse_args()

d3 = args.case.endswith("3")
NX, R = 26, 2.0  # fields do not depend on x; WENO5 needs m + 1 >= 26
case = {
    "run_time_info": "F",
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "y_domain%beg": 0.0,
    "y_domain%end": R,
    "cyl_coord": "T",
    "m": NX - 1,
    "n": args.N - 1,
    "p": args.N - 1 if d3 else 0,
    "dt": 1.0e-8,
    "t_step_start": 0,
    "t_step_stop": 1,
    "t_step_save": 1,
    "num_patches": 1,
    "model_eqns": 2,
    "num_fluids": 1,
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "weno_Re_flux": "F",
    "weno_avg": "F",
    "mapped_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "viscous": args.viscous,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -14 if d3 else -2,
    "bc_y%end": -3,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 1.0,
    "patch_icpp(1)%alpha_rho(1)": 1.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / 0.4,
    "fluid_pp(1)%pi_inf": 0.0,
    "parallel_io": "T",
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "F",
}
if args.viscous == "T":
    case["fluid_pp(1)%Re(1)"] = 1.0  # mu = 1
    if args.case in ("C3", "E2"):
        case["fluid_pp(1)%Re(2)"] = 1.0  # bulk viscosity 1
if d3:
    case.update(
        {
            "z_domain%beg": 0.0,
            "z_domain%end": 2.0 * math.pi,
            "bc_z%beg": -1,
            "bc_z%end": -1,
            "patch_icpp(1)%geometry": 10,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%radius": 2.0 * R,
            "patch_icpp(1)%vel(3)": 0.0,
        }
    )
else:
    case.update({"patch_icpp(1)%geometry": 3, "patch_icpp(1)%x_centroid": 0.5, "patch_icpp(1)%y_centroid": 0.5 * R, "patch_icpp(1)%length_x": 1.0, "patch_icpp(1)%length_y": 2.0 * R})

print(json.dumps(case))
