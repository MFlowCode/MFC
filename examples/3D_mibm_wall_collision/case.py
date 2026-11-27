#!/usr/bin/env python3
"""
Soft-sphere collision of a heavy sphere with a wall: normal restitution and Coulomb friction.

A sphere (R = 5 m, mass 1e6, heavy enough that the gas force is negligible, in air at 1e5 Pa) moves at 3.9 m/s toward the wall y = 0
(bc_y%beg = -15) at an angle --angle from the wall normal and rebounds under collision_model = 1. The rigid-body result for a
sphere on a fixed wall is v_n' = -e v_n, and for the tangential speed, with friction coefficient mu,
    v_t' = v_t - mu (1 + e) |v_n|   if tan(angle) >= 3.5 mu (1 + e)  (sliding throughout the contact),
    v_t' = (5/7) v_t                otherwise                           (sticking before the end of the contact).

    ./mfc.sh run examples/3D_mibm_wall_collision/case.py -n 1 --gpu
    ./mfc.sh run examples/3D_mibm_wall_collision/case.py -n 1 --gpu -- --angle 30 --mu 0.092
    python3 examples/3D_mibm_wall_collision/check_collision.py <case directory> --angle 30 --e 0.98 --mu 0.092

Defaults are a CI-sized run (16^3 cells, 30.5 lead steps, normal impact, e = 0.8, 20 steps per collision). The rebound error of e is first order
in the steps per collision; --cells 40, --lead-steps 200 and --steps-per-collision 80 reproduce the table in the validation page to within the onset-phase effect noted there.
The setup follows the verification cases of Di Maio and Di Renzo (linear frictional-elastic collision model).
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(prog="3D_mibm_wall_collision")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--angle", type=float, default=0.0, help="Impact angle from the wall normal [deg].")
parser.add_argument("--e", type=float, default=0.8, help="coefficient_of_restitution.")
parser.add_argument("--mu", type=float, default=0.0, help="coefficient_of_friction.")
parser.add_argument("--speed", type=float, default=3.9, help="Impact speed [m/s].")
parser.add_argument("--cells", type=int, default=16, help="Cells per direction.")
parser.add_argument("--steps-per-collision", type=float, default=20.0, help="collision_time / dt.")
parser.add_argument("--lead-steps", type=float, default=30.5, help="Free-flight steps before (and after) the impact; a half step puts the contact onset mid-step.")
parser.add_argument("--stride", type=int, default=10, help="ib_force_stride.")
parser.add_argument("--weno_order", type=int, default=1, help="Reconstruction order: 1 is well conditioned (default), 5 is WENO5.")
args = parser.parse_args()

gam = 1.4
dt = 5.0e-6
radius = 5.0
theta = math.radians(args.angle)
domain = 4.0 * radius
collision_time = args.steps_per_collision * dt
lead_distance = args.lead_steps * dt * args.speed

# Free flight to the wall, the contact, and the same flight away
num_steps = round(2.0 * args.lead_steps + args.steps_per_collision)

case = {
    "run_time_info": "F",
    "x_domain%beg": -0.5 * domain,
    "x_domain%end": 0.5 * domain,
    "y_domain%beg": 0.0,
    "y_domain%end": domain,
    "z_domain%beg": -0.5 * domain,
    "z_domain%end": 0.5 * domain,
    "m": args.cells,
    "n": args.cells,
    "p": args.cells,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": num_steps,
    "t_step_save": num_steps,
    "num_patches": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": args.weno_order,
    "avg_state": 2,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -15,
    "bc_y%end": -15,
    "bc_z%beg": -3,
    "bc_z%end": -3,
    "ib": "T",
    "num_ibs": 1,
    "fd_order": 2,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "F",
    "parallel_io": "T",
    "ib_state_wrt": "T",
    "ib_force_wrt": "T",
    "ib_force_stride": args.stride,
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.5 * domain,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": domain,
    "patch_icpp(1)%length_y": domain,
    "patch_icpp(1)%length_z": domain,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": 1.0e5,
    "patch_icpp(1)%alpha_rho(1)": 1.2,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_ib(1)%geometry": 8,
    "patch_ib(1)%x_centroid": -lead_distance * math.sin(theta),
    "patch_ib(1)%y_centroid": radius + lead_distance * math.cos(theta),
    "patch_ib(1)%z_centroid": 0.0,
    "patch_ib(1)%radius": radius,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%mass": 1.0e6,
    "patch_ib(1)%vel(1)": args.speed * math.sin(theta),
    "patch_ib(1)%vel(2)": -args.speed * math.cos(theta),
    "patch_ib(1)%moving_ibm": 2,
    "collision_model": 1,
    "coefficient_of_friction": args.mu,
    "collision_time": collision_time,
    "coefficient_of_restitution": args.e,
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%eos": "ideal_gas",
}

# Above first order, WENO with the usual smoothness weights
if args.weno_order > 1:
    case.update({"weno_eps": 1.0e-16, "mapped_weno": "T", "mp_weno": "T"})

print(json.dumps(case))
