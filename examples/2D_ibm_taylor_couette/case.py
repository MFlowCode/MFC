#!/usr/bin/env python3
"""
2D Taylor-Couette flow between two no-slip immersed walls: accuracy of the wall-normal ghost velocity (ib_wall_order).

The inner circle (radius R = 1, prescribed rotation, moving_ibm = 1) turns at Omega; the outer wall is an STL ring whose inner
edge is at r2 = 2 R and is fixed. Both are no slip. Non-dimensional ideal gas: rho = 1, sound speed 1, Re = Omega R^2 / nu with
constant mu, tip Mach number M = Omega R / c. The gas starts from rest; the steady incompressible solution is
u_theta = A r + B / r (compressibility corrections are O(M^2)); see check_couette.py.

    ./mfc.sh run examples/2D_ibm_taylor_couette/case.py -t pre_process simulation -n 1 -- --cells_per_R 20 --ntau 12 --wall_order 2
    python3 examples/2D_ibm_taylor_couette/check_couette.py examples/2D_ibm_taylor_couette --cells_per_R 20

Defaults are a CI-sized run (12 cells per radius, 64 x 64 cells, 300 fixed steps, far from steady state). Studied values:
M = 0.05, Re = 10, 12 slowest viscous times, 5 to 80 cells per radius (the page validation_wall_order_couette).
Output is serial text (parallel_io = F), so run on one rank.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(prog="2D_ibm_taylor_couette")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--cells_per_R", type=int, default=12, help="Cells per inner radius.")
parser.add_argument("--wall_order", type=int, default=2, help="ib_wall_order (1 mirror, 2 two-probe quadratic).")
parser.add_argument("--ip_interp", type=int, default=1, help="ib_ip_interp (1 inverse distance, 2 TSC).")
parser.add_argument("--visc_compact", type=int, default=1, help="ib_visc_compact: 1 compact viscous gradients next to a body (the solver default), 0 the WENO gradient.")
parser.add_argument("--mach", type=float, default=0.05, help="Tip Mach number Omega R / c.")
parser.add_argument("--re", type=float, default=10.0, help="Re = Omega R^2 / nu.")
parser.add_argument("--sign", type=float, default=1.0, help="Sign of the inner-wall rotation.")
parser.add_argument("--shift", type=float, default=0.0, help="Grid shift in x and y, in cells (0 puts a vertex at the centre).")
parser.add_argument("--steps", type=int, default=300, help="Number of fixed time steps.")
parser.add_argument("--ntau", type=float, default=0.0, help="If > 0, run this many slowest viscous times of the gap instead of --steps.")
parser.add_argument("--n_save", type=int, default=10, help="Number of saves (the checker reads the last two).")
parser.add_argument("--stl", type=str, default="ring.stl", help="Ring STL (inner edge at 2 R); ring_stl.py writes others.")
args = parser.parse_args()

gam, rho, R, r2 = 1.4, 1.0, 1.0, 2.0
p0 = rho / gam  # sound speed 1
omega, nu = args.sign * args.mach / R, args.mach * R / args.re
dx = R / args.cells_per_R
dt = 0.3 * dx / (1.0 + args.mach)  # acoustic CFL 0.3
tau = ((r2 - R) / math.pi) ** 2 / nu  # slowest viscous diffusion time of the gap
steps = math.ceil(args.ntau * tau / dt) if args.ntau > 0 else args.steps
save = max(steps // args.n_save, 1)
steps = save * math.ceil(steps / save)
half = math.ceil((r2 + max(0.3 * R, 8 * dx)) / dx)  # cells from the centre to the edge
lo, hi = (-half + args.shift) * dx, (half + args.shift) * dx

case = {
    "run_time_info": "F",
    "x_domain%beg": lo,
    "x_domain%end": hi,
    "y_domain%beg": lo,
    "y_domain%end": hi,
    "m": 2 * half - 1,
    "n": 2 * half - 1,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": steps,
    "t_step_save": save,
    "num_patches": 1,
    "num_fluids": 1,
    "model_eqns": 2,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-10,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": 2,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "viscous": "T",
    "fd_order": 4,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * (lo + hi),
    "patch_icpp(1)%y_centroid": 0.5 * (lo + hi),
    "patch_icpp(1)%length_x": hi - lo,
    "patch_icpp(1)%length_y": hi - lo,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%eos": "ideal_gas",
    "fluid_pp(1)%Re(1)": 1.0 / (rho * nu),
    "ib": "T",
    "num_ibs": 2,
    "ib_wall_order": args.wall_order,
    "ib_state_wrt": "T",
    "ib_force_wrt": "T",
    "ib_force_stride": max(save // 10, 1),
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": 0.0,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%radius": R,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%moving_ibm": 1,
    "patch_ib(1)%vel(1)": 0.0,
    "patch_ib(1)%vel(2)": 0.0,
    "patch_ib(1)%angular_vel(3)": omega,
    "num_stl_models": 1,
    "stl_models(1)%model_filepath": args.stl,
    "stl_models(1)%model_threshold": 0.5,
    "patch_ib(2)%geometry": 5,
    "patch_ib(2)%model_id": 1,
    "patch_ib(2)%x_centroid": 0.0,
    "patch_ib(2)%y_centroid": 0.0,
    "patch_ib(2)%slip": "F",
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "F",
    "parallel_io": "F",
}
case["ib_ip_interp"] = args.ip_interp  # always written: the solver default is 2
case["ib_visc_compact"] = "T" if args.visc_compact else "F"  # always written

print(json.dumps(case))
