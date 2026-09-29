#!/usr/bin/env python3
"""
3D Taylor-Green vortex at Re = 1600 and a chosen Mach number, for timing the all-Mach pressure
projection against the explicit solver.

Hardcoded IC 380 (as in examples/3D_TaylorGreenVortex) sets a peak velocity U0 = 37.66 and the
incompressible pressure field p0 + rho*U0^2/16 (...), which does not depend on the sound speed. The
Mach number is set through the equation of state instead: a stiffened gas with pi_inf chosen so that
c = U0/M (M = 0.1 is the ideal gas of the original case). Both solvers take a fixed step,
cfl*dx/U0 for the projection and cfl*dx/(U0 + c) for --explicit, so the number of steps per unit of
physical time is known exactly. The convective time is tC = L/U0. See run_sweep.sh for the timing
suite and analyze.py for the kinetic energy history.
"""

import argparse
import json
import math
import sys

parser = argparse.ArgumentParser(description="3D Taylor-Green vortex, all-Mach pressure projection")
parser.add_argument("--mach", type=float, default=0.01, help="Mach number U0/c, at most 0.1 (default: %(default)s)")
parser.add_argument("--N", type=int, default=64, help="cells per direction (default: %(default)s)")
parser.add_argument("--copies", type=int, default=1, help="periods of the vortex along x, for weak scaling (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.5, help="advective CFL, or acoustic with --explicit (default: %(default)s)")
parser.add_argument("--tend", type=float, default=1.0, help="final time in convective times tC (default: %(default)s)")
parser.add_argument("--steps", type=int, default=0, help="run this many steps instead of --tend (default: %(default)s)")
parser.add_argument("--saves", type=int, default=1, help="number of restart/output saves (default: %(default)s)")
parser.add_argument("--info", action="store_true", help="write run_time.inf (costs a reduction and file write per step)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead of the projection")
parser.add_argument("--rdma", action="store_true", help="GPU-aware MPI (rdma_mpi) instead of staging halos through the host")
parser.add_argument("--low-mach", type=int, default=0, choices=[0, 1, 2], help="HLLC low-Mach correction with --explicit (default: %(default)s)")
args, _ = parser.parse_known_args()
if not 0.0 < args.mach <= 0.1:
    parser.error("--mach must be in (0, 0.1]: M = 0.1 is already the ideal gas")

L, Re, rho0, P0, gamma = 1.0, 1600.0, 1.0, 101325.0, 1.4
U0 = 0.1 * math.sqrt(gamma * P0 / rho0)  # the velocity IC 380 sets
c = U0 / args.mach
pi_inf = rho0 * c**2 / gamma - P0  # stiffened gas giving sound speed c at P0
mu = rho0 * U0 * L / Re
tC = L / U0

dx = 2 * math.pi * L / args.N
dt = args.cfl * dx / (U0 + c if args.explicit else U0)
if args.steps > 0:
    Nt = args.steps
else:
    Nt = int(math.ceil(args.tend * tC / dt))
    dt = args.tend * tC / Nt
print(f"dt = {dt:.4e}  steps per tC = {tC / dt:.1f}  acoustic CFL = {dt * (U0 + c) / dx:.1f}  steps = {Nt}", file=sys.stderr)

print(
    json.dumps(
        {
            "run_time_info": "T" if args.info else "F",
            "rdma_mpi": "T" if args.rdma else "F",
            "x_domain%beg": -math.pi * L,
            "x_domain%end": (2 * args.copies - 1) * math.pi * L,
            "y_domain%beg": -math.pi * L,
            "y_domain%end": math.pi * L,
            "z_domain%beg": -math.pi * L,
            "z_domain%end": math.pi * L,
            "m": args.copies * args.N - 1,
            "n": args.N - 1,
            "p": args.N - 1,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": max(Nt // args.saves, 1),
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "bc_z%beg": -1,
            "bc_z%end": -1,
            "viscous": "T",
            "proj_method": "F" if args.explicit else "T",
            "low_Mach": args.low_mach,
            "format": 1,
            "precision": 2,
            "parallel_io": "T",
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%hcid": 380,
            "patch_icpp(1)%x_centroid": (args.copies - 1) * math.pi * L,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_x": 2 * math.pi * L * args.copies,
            "patch_icpp(1)%length_y": 2 * math.pi * L,
            "patch_icpp(1)%length_z": 2 * math.pi * L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": P0,
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": gamma * pi_inf / (gamma - 1.0),
            "fluid_pp(1)%Re(1)": 1.0 / mu,
        }
    )
)
