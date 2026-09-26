#!/usr/bin/env python3
"""
Static 2D water drop in air with surface tension and the all-Mach pressure projection.

A quarter drop of radius R = 0.15 sits at the symmetric corner of a 0.375 x 0.375 box
(the setup of examples/2D_laplace_pressure_jump, with the five-equation model). At rest
the pressure jump is sigma/R; any velocity is a spurious current. See analyze.py.
--explicit runs HLLC at the air acoustic limit as a control.
"""

import argparse
import json
import math
import sys

parser = argparse.ArgumentParser(description="Static drop, all-Mach pressure projection with surface tension")
parser.add_argument("--N", type=int, default=50, help="cells per direction (default: %(default)s)")
parser.add_argument("--sigma", type=float, default=8.0, help="surface tension coefficient (default: %(default)s)")
parser.add_argument("--dt", type=float, default=1.0e-3, help="projection time step (default: %(default)s)")
parser.add_argument("--tend", type=float, default=0.5, help="final time [s] (default: %(default)s)")
parser.add_argument("--saves", type=int, default=5, help="number of output snapshots (default: %(default)s)")
parser.add_argument("--model", default="well_balanced", choices=["conservative", "well_balanced"], help="surface_tension_model (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead, as a control")
args, _ = parser.parse_known_args()

L, R, eps, p0 = 0.375, 0.15, 1.0e-9, 1.0e5
gamma_w, p_inf_w, gamma_a = 2.1, 1.0e6, 1.4
dx = L / args.N
dt = 0.25 * dx / math.sqrt(gamma_a * p0) if args.explicit else args.dt
Nt = int(round(args.tend / dt))
print(f"dt = {dt:.3e} s, {Nt} steps, Laplace jump sigma/R = {args.sigma / R:.3f}", file=sys.stderr)


def patch(i, water):
    a1 = 1 - eps if water else eps
    return {
        f"patch_icpp({i})%x_centroid": 0.0,
        f"patch_icpp({i})%y_centroid": 0.0,
        f"patch_icpp({i})%vel(1)": 0.0,
        f"patch_icpp({i})%vel(2)": 0.0,
        f"patch_icpp({i})%pres": p0,
        f"patch_icpp({i})%alpha_rho(1)": a1 * 1000.0,
        f"patch_icpp({i})%alpha_rho(2)": (1 - a1) * 1.0,
        f"patch_icpp({i})%alpha(1)": a1,
        f"patch_icpp({i})%alpha(2)": 1 - a1,
        f"patch_icpp({i})%cf_val": 1 if water else 0,
    }


print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "y_domain%beg": 0.0,
            "y_domain%end": L,
            "m": args.N - 1,
            "n": args.N - 1,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": max(Nt // args.saves, 1),
            "model_eqns": 2,
            "num_fluids": 2,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -2,
            "bc_x%end": -3,
            "bc_y%beg": -2,
            "bc_y%end": -3,
            "num_patches": 2,
            "surface_tension": "T",
            "sigma": args.sigma,
            "surface_tension_model": "conservative" if args.explicit else args.model,
            "proj_method": "F" if args.explicit else "T",
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            **patch(1, False),
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%length_x": 2.0,
            "patch_icpp(1)%length_y": 2.0,
            **patch(2, True),
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%radius": R,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%smoothen": "T",
            "patch_icpp(2)%smooth_patch_id": 1,
            "patch_icpp(2)%smooth_coeff": 0.95,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma_w - 1.0),
            "fluid_pp(1)%pi_inf": gamma_w * p_inf_w / (gamma_w - 1.0),
            "fluid_pp(2)%eos": "ideal_gas",
            "fluid_pp(2)%gamma": 1.0 / (gamma_a - 1.0),
        }
    )
)
