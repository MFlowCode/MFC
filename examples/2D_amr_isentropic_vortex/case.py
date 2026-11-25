#!/usr/bin/env python3
"""
Isentropic vortex convected through AMR blocks (see README.md).

Shu's vortex (eps = 5, gamma = 1.4) in a uniform stream u = (1, 0) on a periodic [-10, 10]^2, from x = -3 to x = 3 at t = 6.
--mode static: one 2:1 block over [-1.25, 1.25] x [-2.5, 2.5], which the vortex enters and leaves. --mode regrid: blocks
follow the vortex. --mode uniform: no AMR, at --N cells or, with --fine, 2N. Every run at the same --N takes the same fine-level
step (dt for a 2N grid at CFL 0.4), so the runs differ only in space.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="2D isentropic vortex through AMR blocks")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("--N", type=int, default=128, help="Base cells per direction (default: 128)")
parser.add_argument("--mode", choices=["uniform", "static", "regrid"], default="static")
parser.add_argument("--fine", action="store_true", help="uniform: run at 2N cells")
parser.add_argument("--prolong", type=int, default=1, help="amr_prolong: 1 linear, 2 quartic (default: 1)")
parser.add_argument("--subcycle", action="store_true", help="Subcycle the refined level")
parser.add_argument("--t-end", type=float, default=6.0, help="End time (default: 6, vortex at x = 3)")
parser.add_argument("--nsave", type=int, default=1, help="Saves over the run (default: 1)")
args = parser.parse_args()

L, gam, eps, x0 = 20.0, 1.4, 5.0, -3.0
N = args.N
cells = 2 * N if args.fine else N

# Fine-level step: CFL 0.4 on a 2N grid at |u| + c <= 3; subcycling takes two of them per level-0 step.
nt = 2 * math.ceil(args.t_end / (0.4 * L / (2 * N) / 3.0) / 2)
nt //= 2 if args.subcycle else 1
dt = args.t_end / nt

# Shu's vortex at (x0, 0) as exact cell averages (3x3-point Gauss, sixth order) of the conserved variables, so that every level
# starts from its own averages. The primitives handed to the patch reproduce those averages: u = <rho u>/<rho>, and so on.
G3 = ((-math.sqrt(0.6), 5.0 / 18.0), (0.0, 8.0 / 18.0), (math.sqrt(0.6), 5.0 / 18.0))


def f(v):
    """Exponent form, which the expression translator gives kind wp."""
    return f"{v:.16e}"


def avg(g):
    return "(" + " + ".join(f"{f(wa * wb)} * {g(f'(x + {f(a / 2)} * dx)', f'(y + {f(b / 2)} * dy)')}" for a, wa in G3 for b, wb in G3) + ")"


def fields(X, Y):
    r2 = f"(({X} - ({x0})) ** 2 + {Y} ** 2)"
    T = f"(1.0 - {f((gam - 1.0) * eps**2 / (8.0 * gam * math.pi**2))} * exp(1.0 - {r2}))"
    du = f"({f(eps / (2.0 * math.pi))} * exp(0.5 * (1.0 - {r2})))"
    u, v = f"(1.0 - {Y} * {du})", f"(({X} - ({x0})) * {du})"
    rho = f"{T} ** {f(1.0 / (gam - 1.0))}"
    return {"rho": rho, "mx": f"{rho} * {u}", "my": f"{rho} * {v}", "E": f"({T} ** {f(gam / (gam - 1.0))} / {f(gam - 1.0)} + 0.5 * {rho} * ({u} ** 2 + {v} ** 2))"}


A = {k: avg(lambda X, Y, k=k: fields(X, Y)[k]) for k in ("rho", "mx", "my", "E")}
ic = {
    "patch_icpp(1)%alpha_rho(1)": A["rho"],
    "patch_icpp(1)%vel(1)": f"{A['mx']} / {A['rho']}",
    "patch_icpp(1)%vel(2)": f"{A['my']} / {A['rho']}",
    "patch_icpp(1)%pres": f"{f(gam - 1.0)} * ({A['E']} - 0.5 * ({A['mx']} ** 2 + {A['my']} ** 2) / {A['rho']})",
}

amr = {}
if args.mode != "uniform":
    amr = {
        "amr": "T",
        "amr_block_beg(1)": 7 * N // 16,
        "amr_block_end(1)": 9 * N // 16 - 1,
        "amr_block_beg(2)": 3 * N // 8,
        "amr_block_end(2)": 5 * N // 8 - 1,
        "amr_regrid_int": 0,
        "amr_prolong": args.prolong,
        "amr_subcycle": "T" if args.subcycle else "F",
        "amr_max_grid_size": N // 2,
    }
    if args.mode == "regrid":
        # tag the density jump down to r = 2.5 (dRho/dr there is 0.006), so the refined region is the same at every N
        amr.update({"amr_regrid_int": 10, "amr_tag_eps": 0.006 * L / N, "amr_buf": 3})

print(
    json.dumps(
        {
            "run_time_info": "F",
            "x_domain%beg": -L / 2,
            "x_domain%end": L / 2,
            "y_domain%beg": -L / 2,
            "y_domain%end": L / 2,
            "m": cells - 1,
            "n": cells - 1,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": nt,
            "t_step_save": max(1, nt // args.nsave),
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "mpp_lim": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "F",
            "parallel_io": "T",
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": L,
            **ic,
            "patch_icpp(1)%alpha(1)": 1.0,
            "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
            **amr,
        }
    )
)
