#!/usr/bin/env python3
"""
Galilean twin of a circle in a doubly periodic box: the gas-mass budget, heat and force of a moving immersed boundary.

Non-dimensional ideal gas (rho = 1, p = 1/1.4, sound speed 1, D = 1), no chemistry. --frame moving: the circle translates at
Mach --mach through gas at rest (moving_ibm = 1). --frame fixed: the circle is held and the gas moves at -U, the same relative
flow. In the closed box the gas mass is constant, so any change of it is created by the body: check_budget.py reads it from the
saves. --Re > 0 adds viscosity and, with --Twall, an isothermal wall (conduction heat, compared between the twins).
Without --relax the solver default applies (ib_mass_relax resolves to 30 for a moving body and 0 for a fixed one).

    ./mfc.sh run examples/2D_ibm_moving_mass_budget/case.py -t pre_process simulation -n 1 --gpu -- \
        --cells_per_D 24 --box 8 --travel 6 --nsave 12 --relax 0
    python3 examples/2D_ibm_moving_mass_budget/check_budget.py mass examples/2D_ibm_moving_mass_budget

Defaults are a CI-sized run: 24 cells per diameter in a box of 5 D, 100 fixed steps (a quarter of a diameter of travel).
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(prog="2D_ibm_moving_mass_budget")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--frame", choices=("moving", "fixed"), default="moving", help="Body moves through gas at rest, or is held in gas moving at -U.")
parser.add_argument("--cells_per_D", type=int, default=24, help="Cells across the diameter.")
parser.add_argument("--box", type=float, default=5.0, help="Box side [D].")
parser.add_argument("--mach", type=float, default=0.2, help="Relative Mach number U / c.")
parser.add_argument("--Re", type=float, default=0.0, help="Re = U D / nu (0: inviscid).")
parser.add_argument("--Pr", type=float, default=0.7, help="Prandtl number (viscous runs).")
parser.add_argument("--Twall", type=float, default=0.0, help="Isothermal wall, T_wall / T_gas (viscous runs; 0: none).")
parser.add_argument("--steps", type=int, default=100, help="Number of fixed steps.")
parser.add_argument("--travel", type=float, default=0.0, help="Distance travelled relative to the gas [D]; > 0 sets the number of steps.")
parser.add_argument("--nsave", type=int, default=1, help="Number of saves.")
parser.add_argument("--cfl", type=float, default=0.4, help="Acoustic CFL number of the fixed step.")
parser.add_argument("--relax", type=int, default=None, help="ib_mass_relax (omitted: the solver default).")
parser.add_argument("--ip_interp", type=int, default=1, choices=(1, 2), help="Image-point interpolation (1: inverse distance, 2: TSC).")
parser.add_argument("--wall_order", type=int, default=2, choices=(1, 2), help="ib_wall_order (1 wall value, 2 quadratic, the solver default).")
parser.add_argument("--visc_compact", type=int, default=1, choices=(0, 1), help="ib_visc_compact (1 compact viscous gradient next to the body, the solver default; 0 the WENO gradient).")
parser.add_argument("--slip", action="store_true", help="Slip wall.")
parser.add_argument("--stride", type=int, default=10, help="ib_force_stride: steps between force rows.")
parser.add_argument("--serial_io", action="store_true", help="parallel_io = F (a build without MPI); the checker then reads p_all/.")
args = parser.parse_args()

gam, rho, p0, c, D = 1.4, 1.0, 1.0 / 1.4, 1.0, 1.0
L = args.box * D
N = round(args.box * args.cells_per_D)
dx = L / N
U = args.mach * c
T_gas = p0 / ((gam - 1.0) * rho)  # cv = 1
dt = args.cfl * dx / (c + U)
steps = math.ceil(args.travel * D / (U * dt)) if args.travel > 0 else args.steps
steps = args.nsave * math.ceil(steps / args.nsave)
moving = args.frame == "moving"
# Off the grid lines by fractions of dx, so the body crosses cell faces at different phases in x and y. Not dx / 3 in x: at the
# default Mach and CFL a step moves the body dx / 15, so every Runge-Kutta stage would sit on a multiple of dx / 30 and some on a
# cell-centre line, where round-off picks the image-point stencil.
xc, yc = 0.5 * L + 0.1371 * dx, 0.5 * L + dx / 5.0

case = {
    "run_time_info": "F",
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "y_domain%beg": 0.0,
    "y_domain%end": L,
    "m": N - 1,
    "n": N - 1,
    "p": 0,
    "cyl_coord": "F",
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": steps,
    "t_step_save": steps // args.nsave,
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
    "null_weights": "F",
    "mp_weno": "F",
    "avg_state": 2,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * L,
    "patch_icpp(1)%y_centroid": 0.5 * L,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%length_y": L,
    "patch_icpp(1)%vel(1)": 0.0 if moving else -U,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0,
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": xc,
    "patch_ib(1)%y_centroid": yc,
    "patch_ib(1)%radius": 0.5 * D,
    "patch_ib(1)%slip": "T" if args.slip else "F",
    "ib_ip_interp": args.ip_interp,
    "ib_wall_order": args.wall_order,
    "ib_state_wrt": "T",
    "ib_force_wrt": "T",
    "ib_force_stride": args.stride,
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    "fd_order": 2,
    "parallel_io": "F" if args.serial_io else "T",
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "F",
}
if moving:
    case.update({"patch_ib(1)%moving_ibm": 1, "patch_ib(1)%vel(1)": U, "patch_ib(1)%vel(2)": 0.0})
if args.Re > 0:
    mu = rho * U * D / args.Re
    case.update({"viscous": "T", "weno_Re_flux": "T", "weno_avg": "T", "fd_order": 4, "fluid_pp(1)%Re(1)": 1.0 / mu})
    case["ib_visc_compact"] = "T" if args.visc_compact else "F"
    if args.Twall > 0:
        case.update(
            {
                "fluid_pp(1)%cv": 1.0,
                "fluid_pp(1)%k_therm": gam * mu / args.Pr,
                "patch_ib(1)%thermal_bc": 1,
                "patch_ib(1)%Twall": args.Twall * T_gas,
            }
        )
if args.relax is not None:
    case["ib_mass_relax"] = args.relax

if __name__ == "__main__":
    print(json.dumps(case))
