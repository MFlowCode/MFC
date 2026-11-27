#!/usr/bin/env python3
"""
Planar shock through a fixed curtain of no-slip cylinders (2D) or spheres (3D), periodic across the flow.

A Mach 1.66 shock in air (p1 = 82.7 kPa, T1 = 291.9 K) hits a random curtain of diameter D = 115 um bodies at solid fraction
phi = 0.2 (the configuration of Sridhar, Fox and Capecelatro 2025, Case 5, which is also the geometry of a shock-tube curtain
experiment). The gas is a single-fluid ideal gas, inviscid by default: an inviscid no-slip wall (slip F, viscous F) reproduces the
viscous mean flow through a dense curtain, the inviscid slip wall (--slip) does not. --re_ps makes the run viscous with a constant
mu from Re_ps = rho2 u2 D / mu. The shock sits half a diameter upstream of the curtain at t = 0. Fixed bodies by default;
--moving frees them (no slip, mass from --rho_p, soft-sphere collisions). The mean gas velocity in the curtain is the
quantity checked: check_curtain.py compares the Favre-averaged profile with the published one.

    ./mfc.sh run examples/2D_ibm_shock_curtain/case.py -n 1 --gpu                          # CI-sized 2D run
    ./mfc.sh run examples/2D_ibm_shock_curtain/case.py -n 1 --gpu -- --dim 3 --single --mach 1.22 --xu 3 --xd 7 \
        --width 6 --cells_per_D 12 --steps 0 --tau 8   # one sphere: unsteady drag (check_curtain.py --drag)
    ./mfc.sh run examples/2D_ibm_shock_curtain/case.py -t pre_process simulation post_process -n 1 --gpu -- \
        --dim 3 --h 17.4 --width 12 --xu 7 --xd 23 --cells_per_D 12 --tau 2 --silo
    python3 examples/2D_ibm_shock_curtain/check_curtain.py <case directory>

Defaults are the CI-sized run: 12 cylinders over a curtain 6 D thick and 8 D wide, 12 cells per diameter, 240 x 96 cells, 300
fixed steps (the shock is inside the curtain at the end). The 3D flags above are the reference configuration (30 D x 12 D x 12 D,
957 spheres, 7.5 M cells at 12 cells per diameter).
"""

import argparse
import importlib.util
import json
import math
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("periodic_array", os.path.join(here, "..", "2D_ibm_periodic_array", "case.py"))
array = importlib.util.module_from_spec(spec)
spec.loader.exec_module(array)

parser = argparse.ArgumentParser(prog="2D_ibm_shock_curtain")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--dim", type=int, default=2, choices=(2, 3), help="Space dimension.")
parser.add_argument("--mach", type=float, default=1.66, help="Shock Mach number.")
parser.add_argument("--phi", type=float, default=0.2, help="Solid (area in 2D) fraction in the curtain.")
parser.add_argument("--h", type=float, default=6.0, help="Curtain thickness [D].")
parser.add_argument("--width", type=float, default=8.0, help="Periodic cross-stream width [D].")
parser.add_argument("--xu", type=float, default=6.0, help="Domain extent upstream of the curtain [D].")
parser.add_argument("--xd", type=float, default=14.0, help="Domain extent downstream of the curtain's upstream face [D], the curtain included.")
parser.add_argument("--gap", type=float, default=0.2, help="Minimum initial surface gap [D].")
parser.add_argument("--seed", type=int, default=1)
parser.add_argument("--cells_per_D", type=float, default=12.0)
parser.add_argument("--steps", type=int, default=300, help="Fixed number of steps (CI-sized run); 0 selects the adaptive step with --tau.")
parser.add_argument("--tau", type=float, default=2.0, help="Run length t / tau_L, tau_L = h D / W_s (D / W_s with --single), when --steps is 0.")
parser.add_argument("--frames", type=int, default=10, help="Saves of the adaptive run (a fixed-step run saves once, at the end).")
parser.add_argument("--cfl", type=float, default=0.4)
parser.add_argument("--slip", action="store_true", help="Slip walls (default no slip).")
parser.add_argument("--wall_order", type=int, default=2, choices=(1, 2), help="ib_wall_order (1 wall value, 2 quadratic, the solver default).")
parser.add_argument("--re_ps", type=float, default=None, help="Viscous: constant mu from Re_ps = rho2 u2 D / mu (bulk viscosity 0.6 mu).")
parser.add_argument("--moving", action="store_true", help="Free bodies with collisions (default fixed).")
parser.add_argument("--rho_p", type=float, default=2520.0, help="Body density [kg/m^3] (--moving).")
parser.add_argument("--single", action="store_true", help="One body at the shock's arrival point instead of a curtain (unsteady drag, check_curtain.py --drag).")
parser.add_argument("--stride", type=int, default=20, help="ib_force_stride.")
parser.add_argument("--silo", action="store_true", help="Write density, pressure and x-velocity for post_process (check_curtain.py).")
args = parser.parse_known_args()[0]

D, g, R = 115.0e-6, 1.4, 287.05
p1, T1 = 82.7e3, 291.9
a1, rho1 = math.sqrt(g * R * T1), p1 / (R * T1)
M2 = args.mach**2
p2 = p1 * (1.0 + 2.0 * g / (g + 1.0) * (M2 - 1.0))
rho2 = rho1 * (g + 1.0) * M2 / ((g - 1.0) * M2 + 2.0)
Ws = args.mach * a1
u2 = Ws * (1.0 - rho1 / rho2)
dim, dx = args.dim, D / args.cells_per_D
xlo, xhi, Wd = -args.xu * D, args.xd * D, args.width * D
cells = lambda length: round(length / dx)  # noqa: E731
vol = math.pi / 4.0 * D**2 if dim == 2 else math.pi / 6.0 * D**3
n = round(args.phi * args.h * D * Wd ** (dim - 1) / vol)
d_min = (1.0 + args.gap) * D
lo = [0.5 * D] + [0.0] * (dim - 1)
hi = [args.h * D - 0.5 * D] + [Wd] * (dim - 1)
if args.single:
    n, centers, x_shock = 1, [[0.5 * D] + [0.5 * Wd] * (dim - 1)], 0.0  # the shock reaches the nose at t = 0
else:
    centers, x_shock = array.pack(n, lo, hi, range(1, dim), d_min, args.seed), -0.5 * D
tau_L = (1.0 if args.single else args.h) * D / Ws  # h D / W_s for a curtain, D / W_s for one body


def box(idx, a, b, p, rho, u):
    out = {f"patch_icpp({idx})%geometry": 3 if dim == 2 else 9, f"patch_icpp({idx})%pres": p, f"patch_icpp({idx})%alpha_rho(1)": rho, f"patch_icpp({idx})%alpha(1)": 1.0}
    for k, name in enumerate("xyz"[:dim]):
        e0, e1 = (a, b) if k == 0 else (0.0, Wd)
        out.update({f"patch_icpp({idx})%{name}_centroid": 0.5 * (e0 + e1), f"patch_icpp({idx})%length_{name}": e1 - e0, f"patch_icpp({idx})%vel({k + 1})": u if k == 0 else 0.0})
    return out


case = {
    "run_time_info": "T",
    "cyl_coord": "F",
    "x_domain%beg": xlo,
    "x_domain%end": xhi,
    "m": cells(xhi - xlo) - 1,
    "n": cells(Wd) - 1,
    "p": cells(Wd) - 1 if dim == 3 else 0,
    "num_patches": 2,
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
    "weno_avg": "T",
    "avg_state": 2,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "viscous": "F" if args.re_ps is None else "T",
    "fd_order": 4,
    "fluid_pp(1)%gamma": 1.0 / (g - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    # Post-shock inflow held at the left face, extrapolation outflow, periodic across the flow
    "bc_x%beg": -17,
    "bc_x%end": -3,
    "ib": "T",
    "ib_wall_order": args.wall_order,
    "num_ibs": n,
    "many_ib_patch_parallelism": "T",
    "ib_state_wrt": "T",
    "ib_force_wrt": "T",
    "ib_force_stride": args.stride,
    "parallel_io": "F",
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "F",
}
for name in "yz"[: dim - 1]:
    case.update({f"{name}_domain%beg": 0.0, f"{name}_domain%end": Wd, f"bc_{name}%beg": -1, f"bc_{name}%end": -1})
case.update(box(1, xlo, x_shock, p2, rho2, u2))
case.update(box(2, x_shock, xhi, p1, rho1, 0.0))
unit = math.pi / 4.0 * D**2 if dim == 2 else math.pi / 6.0 * D**3
for i, c in enumerate(centers):
    j = i + 1
    case.update({f"patch_ib({j})%geometry": 2 if dim == 2 else 8, f"patch_ib({j})%radius": 0.5 * D, f"patch_ib({j})%slip": "T" if args.slip else "F"})
    case[f"patch_ib({j})%moving_ibm"] = 2 if args.moving else 0
    if args.moving:
        case[f"patch_ib({j})%mass"] = args.rho_p * unit
    for name, x in zip("xyz", c):
        case[f"patch_ib({j})%{name}_centroid"] = x
if args.moving:
    case.update({"collision_model": 1, "coefficient_of_restitution": 0.9, "collision_time": 20.0 * args.cfl * dx / (u2 + a1), "ib_coefficient_of_friction": 0.1})
if args.re_ps is not None:
    mu = rho2 * u2 * D / args.re_ps
    case.update({"fluid_pp(1)%Re(1)": 1.0 / mu, "fluid_pp(1)%Re(2)": 1.0 / (0.6 * mu)})
if args.silo:
    case.update({"rho_wrt": "T", "pres_wrt": "T", "vel_wrt(1)": "T"})
if args.steps > 0:
    case.update({"dt": args.cfl * dx / (u2 + math.sqrt(g * p2 / rho2)), "t_step_start": 0, "t_step_stop": args.steps, "t_step_save": args.steps})
else:
    t_stop = -x_shock / Ws + args.tau * tau_L
    case.update({"cfl_adap_dt": "T", "cfl_target": args.cfl, "n_start": 0, "t_save": t_stop / args.frames, "t_stop": t_stop})

if __name__ == "__main__":
    print(f"{dim}D  Ms {args.mach:g}  N {n}  Ws {Ws:.1f} m/s  u2 {u2:.1f} m/s  p2/p1 {p2 / p1:.3f}  tau_L {tau_L:.3e} s  grid {case['m'] + 1} x {case['n'] + 1} x {case['p'] + 1}", file=sys.stderr)
    print(json.dumps(case))
