#!/usr/bin/env python3
"""
2D plane Poiseuille flow between two flat no-slip immersed walls: accuracy of the wall-normal ghost velocity (ib_wall_order).

The walls are two static rectangle immersed boundaries at y = 0 and y = H = 1, periodic in x, driven by the body acceleration g_x.
Non-dimensional ideal gas: rho = 1, sound speed 1, constant mu. The exact steady incompressible solution is
u = g y (H - y) / (2 nu), with peak u_max = g H^2 / (8 nu) (equal to the Mach number) and wall shear rho g H / 2.

    ./mfc.sh run examples/2D_ibm_plane_poiseuille/case.py -t pre_process simulation -n 1 -- --cells 16 --ntau 12 --wall_order 2
    python3 examples/2D_ibm_plane_poiseuille/check_poiseuille.py examples/2D_ibm_plane_poiseuille --cells 16

Defaults are a CI-sized run (8 cells across the gap, 32 x 28 cells, 300 fixed steps, far from steady state). Studied values:
u_max = 0.05, Re = u_max H / nu = 1, 12 slowest viscous times, 4 to 64 cells across the gap, wall on a cell face or 0.37 cell
inside it (the page validation_wall_order_couette). Output is serial text (parallel_io = F), so run on one rank.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(prog="2D_ibm_plane_poiseuille")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--cells", type=int, default=8, help="Cells across the gap H = 1.")
parser.add_argument("--wall_order", type=int, default=2, help="ib_wall_order (1 mirror, 2 two-probe quadratic).")
parser.add_argument("--mach", type=float, default=0.05, help="Peak velocity over the sound speed.")
parser.add_argument("--re", type=float, default=1.0, help="Re = u_max H / nu.")
parser.add_argument("--offset", type=float, default=0.0, help="Wall offset from a cell face, in cells.")
parser.add_argument("--nx", type=int, default=32, help="Cells in x (periodic; at least 25).")
parser.add_argument("--steps", type=int, default=300, help="Number of fixed time steps.")
parser.add_argument("--ntau", type=float, default=0.0, help="If > 0, run this many slowest viscous times of the gap instead of --steps.")
parser.add_argument("--n_save", type=int, default=10, help="Number of saves (the checker reads the last two).")
args = parser.parse_args()

gam, rho, H = 1.4, 1.0, 1.0
p0 = rho / gam  # sound speed 1
u_max = args.mach
nu = u_max * H / args.re
accel = 8.0 * nu * u_max / H**2
dx = H / args.cells
dt = 0.3 * dx / (1.0 + u_max)  # acoustic CFL 0.3
tau = (H / math.pi) ** 2 / nu  # slowest viscous diffusion time of the gap
steps = math.ceil(args.ntau * tau / dt) if args.ntau > 0 else args.steps
save = max(steps // args.n_save, 1)
steps = save * math.ceil(steps / save)
margin = max(6, (26 - args.cells) // 2 + 1)  # cells behind each wall: at least 25 cells in y (WENO5 stencil rule)
lo = -(margin + args.offset) * dx
ny = args.cells + 2 * margin
hi = lo + ny * dx
Lx = args.nx * dx
slab = 1.5 * (margin + args.offset) * dx  # wall thickness; the centroid stays inside the domain

case = {
    "run_time_info": "F",
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": lo,
    "y_domain%end": hi,
    "m": args.nx - 1,
    "n": ny - 1,
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
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -16,
    "bc_y%end": -16,
    "bf_x": "T",
    "g_x": accel,
    "k_x": 0.0,
    "w_x": 0.0,
    "p_x": 0.0,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * Lx,
    "patch_icpp(1)%y_centroid": 0.5 * (lo + hi),
    "patch_icpp(1)%length_x": Lx,
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
    "patch_ib(1)%geometry": 3,
    "patch_ib(1)%x_centroid": 0.5 * Lx,
    "patch_ib(1)%y_centroid": -0.5 * slab,
    "patch_ib(1)%length_x": 2.0 * Lx,
    "patch_ib(1)%length_y": slab,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%mass": 0.0,
    "patch_ib(2)%geometry": 3,
    "patch_ib(2)%x_centroid": 0.5 * Lx,
    "patch_ib(2)%y_centroid": H + 0.5 * slab,
    "patch_ib(2)%length_x": 2.0 * Lx,
    "patch_ib(2)%length_y": slab,
    "patch_ib(2)%slip": "F",
    "patch_ib(2)%mass": 0.0,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "F",
    "parallel_io": "F",
}

print(json.dumps(case))
