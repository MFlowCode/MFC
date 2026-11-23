#!/usr/bin/env python3
"""
Bow-shock stand-off of a cylinder (--dim 2) or a sphere (--dim 3) in supersonic flow.

Ideal gas (gamma 1.4, rho = 1, sound speed 1, D = 1) at Mach --mach streaming over a fixed no-slip body (immersed boundary), an
impulsive start with extrapolation on every face, inviscid unless --Re is given. standoff.py measures the distance of the
detached shock from the body on the stagnation line and compares it with Billig's correlation.

    ./mfc.sh run examples/2D_ibm_cylinder_standoff/case.py -t pre_process simulation -n 1 --gpu -- --dim 3 --cells_per_D 30 --steps 0 --tstar 30
    python3 examples/2D_ibm_cylinder_standoff/standoff.py examples/2D_ibm_cylinder_standoff

Defaults are a CI-sized run: the cylinder at Mach 2 with 16 cells per diameter in a 4 D box and 150 fixed steps, before the shock
has settled. The table runs use --steps 0 (a time step from --cfl and the run length --tstar in R / u) and save --frames times.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(prog="2D_ibm_cylinder_standoff")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--dim", type=int, default=2, choices=(2, 3), help="2: cylinder, 3: sphere.")
parser.add_argument("--mach", type=float, default=2.0, help="Free-stream Mach number.")
parser.add_argument("--cells_per_D", type=int, default=16, help="Cells across the diameter.")
parser.add_argument("--xu", type=float, default=1.5, help="Upstream extent [D].")
parser.add_argument("--xd", type=float, default=2.5, help="Downstream extent [D].")
parser.add_argument("--hy", type=float, default=2.0, help="Lateral half width [D].")
parser.add_argument("--Re", type=float, default=0.0, help="Re = rho U D / mu (0: inviscid).")
parser.add_argument("--slip", action="store_true", help="Slip wall.")
parser.add_argument("--ip_interp", type=int, default=1, choices=(1, 2), help="Image-point interpolation (1: inverse distance, 2: TSC).")
parser.add_argument("--steps", type=int, default=150, help="Fixed number of steps; 0 runs --tstar.")
parser.add_argument("--tstar", type=float, default=30.0, help="Run length in R / u when --steps is 0.")
parser.add_argument("--frames", type=int, default=6, help="Saves of the run of --tstar.")
parser.add_argument("--cfl", type=float, default=0.4, help="Acoustic CFL number.")
parser.add_argument("--serial_io", action="store_true", help="parallel_io = F (a build without MPI); the checker then reads p_all/.")
args = parser.parse_args()

gam, D = 1.4, 1.0
u = args.mach  # c = 1
dx = D / args.cells_per_D
Lx, Ly = args.xu + args.xd, 2.0 * args.hy
m, n = round(Lx / dx) - 1, round(Ly / dx) - 1
dt = args.cfl * dx / (u + 1.0)
if args.steps > 0:
    steps = args.steps
else:
    steps = args.frames * math.ceil(args.tstar * 0.5 * D / u / dt / args.frames)
    dt = args.tstar * 0.5 * D / u / steps
three = args.dim == 3

case = {
    "run_time_info": "F",
    "x_domain%beg": -args.xu,
    "x_domain%end": args.xd,
    "y_domain%beg": -args.hy,
    "y_domain%end": args.hy,
    "m": m,
    "n": n,
    "p": n if three else 0,
    "cyl_coord": "F",
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": steps,
    "t_step_save": steps if args.steps > 0 else steps // args.frames,
    "num_patches": 1,
    "num_fluids": 1,
    "model_eqns": 2,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "avg_state": 2,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "patch_icpp(1)%geometry": 9 if three else 3,
    "patch_icpp(1)%x_centroid": 0.5 * (args.xd - args.xu),
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": u,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 1.0 / gam,
    "patch_icpp(1)%alpha_rho(1)": 1.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 8 if three else 2,
    "patch_ib(1)%x_centroid": 0.0,
    "patch_ib(1)%y_centroid": 0.5 * dx,  # on a cell-centre row: the stagnation line is a grid line
    "patch_ib(1)%radius": 0.5 * D,
    "patch_ib(1)%slip": "T" if args.slip else "F",
    "ib_ip_interp": args.ip_interp,
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    "fd_order": 2,
    "parallel_io": "F" if args.serial_io else "T",
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "F",
}
if three:
    case.update(
        {
            "z_domain%beg": -args.hy,
            "z_domain%end": args.hy,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_z": Ly,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_ib(1)%z_centroid": 0.5 * dx,
        }
    )
if args.Re > 0:
    case.update({"viscous": "T", "weno_Re_flux": "F", "fluid_pp(1)%Re(1)": args.Re / (u * D)})

if __name__ == "__main__":
    print(json.dumps(case))
