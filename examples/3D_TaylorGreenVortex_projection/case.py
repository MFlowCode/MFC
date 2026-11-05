#!/usr/bin/env python3
"""
3D Taylor-Green vortex at Re = 1600 and a chosen Mach number, for timing the all-Mach pressure
projection against the explicit solver.

Hardcoded IC 380 (as in examples/3D_TaylorGreenVortex) sets a peak velocity U0 = 37.66 and the
incompressible pressure field p0 + rho*U0^2/16 (...), which does not depend on the sound speed. The
Mach number is set through the equation of state instead: a stiffened gas with pi_inf chosen so that
c = U0/M (M = 0.1 is the ideal gas of the original case); M = 0.01 here. Both solvers take a fixed
step, cfl*dx/U0 for the projection and cfl*dx/(U0 + c) for --explicit, so the number of steps per unit
of physical time is known exactly. The convective time is tC = L/U0.
"""

import argparse
import json
import math
import sys

parser = argparse.ArgumentParser(description="3D Taylor-Green vortex, all-Mach pressure projection")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC (low_Mach = 2) at the acoustic limit instead of the projection")
args, _ = parser.parse_known_args()

mach, N, cfl, tend, saves = 0.01, 64, 0.5, 10.0, 20  # Mach number, cells per direction, CFL, end time in tC, outputs
L, Re, rho0, P0, gamma = 1.0, 1600.0, 1.0, 101325.0, 1.4
U0 = 0.1 * math.sqrt(gamma * P0 / rho0)  # the velocity IC 380 sets
c = U0 / mach
pi_inf = rho0 * c**2 / gamma - P0  # stiffened gas giving sound speed c at P0
mu = rho0 * U0 * L / Re
tC = L / U0

dx = 2 * math.pi * L / N
dt = cfl * dx / (U0 + c if args.explicit else U0)
Nt = int(math.ceil(tend * tC / dt))
dt = tend * tC / Nt
print(f"dt = {dt:.4e}  steps per tC = {tC / dt:.1f}  acoustic CFL = {dt * (U0 + c) / dx:.1f}  steps = {Nt}", file=sys.stderr)


print(
    json.dumps(
        {
            "run_time_info": "T",
            **{f"{ax}_domain%{e}": v * math.pi * L for ax in "xyz" for e, v in (("beg", -1), ("end", 1))},
            **{ext: N - 1 for ext in "mnp"},
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": max(Nt // saves, 1),
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
            **({"low_Mach": 2} if args.explicit else {"proj_method": "T"}),
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%hcid": 380,
            **{f"patch_icpp(1)%{ax}_centroid": 0.0 for ax in "xyz"},
            **{f"patch_icpp(1)%length_{ax}": 2 * math.pi * L for ax in "xyz"},
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
