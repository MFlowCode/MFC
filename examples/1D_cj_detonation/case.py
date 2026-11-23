#!/usr/bin/env python3
"""
1D Chapman-Jouguet detonation in argon-diluted H2/O2 (2H2 + O2 + 7Ar, 100 kPa, 298 K) from a closed end, for the front speed and the
post-detonation state against the equilibrium CJ solution (cj_ref.json, written by cj_ref.py).

    ./mfc.sh run examples/1D_cj_detonation/case.py -n 1 -t pre_process simulation post_process -- --ic tz --ppi 5 --L 0.5
    python3 examples/1D_cj_detonation/check_cj.py examples/1D_cj_detonation --case-args "--ic tz --ppi 5 --L 0.5"

x = 0 is a reflecting wall. Starts (--ic):
  tz      the Taylor-Zeldovich solution at t0 = 40 induction lengths / D_CJ (u = 0 plateau, the fan as 7 constant steps, then a slab of
          fresh gas at the von Neumann state --vnw induction lengths thick); the wave then approaches D_CJ from below
  driver  a short slug of equilibrium products at --pdrive times the explosion pressure against the wall (the CI test: front formation)
Numerics follow examples/2D_detonation_cell (5-equation model, mapped WENO5 + mp_weno, HLLC, SSP-RK3, split alpha-QSS reactions with
adaptive substeps), no diffusion or viscosity, with a fixed dt from the frozen vN state: dt = cfl dx / (u_vN + c_vN).
"""

import argparse
import json
import math
import os

import cantera as ct
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MECH = "h2o2.yaml"

parser = argparse.ArgumentParser(prog="1D_cj_detonation")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--ic", default="driver", choices=("driver", "tz"), help="Initial condition.")
parser.add_argument("--ppi", type=float, default=2.5, help="Cells per ZND induction length.")
parser.add_argument("--L", type=float, default=0.10, help="Domain length [m].")
parser.add_argument("--cfl", type=float, default=0.3, help="CFL on u_vN + c_vN.")
parser.add_argument("--nsave", type=int, default=200, help="Saves over the run.")
parser.add_argument("--steps", type=int, default=0, help="Total steps (0: until the front reaches 0.9 L).")
parser.add_argument("--vnw", type=float, default=5.0, help="tz: von Neumann slab thickness [induction lengths].")
parser.add_argument("--pdrive", type=float, default=10.0, help="driver: slug pressure over the explosion pressure.")
parser.add_argument("--driver-len", type=float, default=20.0, help="driver: slug length [induction lengths].")
parser.add_argument("--substeps", type=int, default=2, help="alpha-QSS sub-step floor.")
parser.add_argument("--substeps-max", type=int, default=16, help="alpha-QSS sub-step ceiling.")
args = parser.parse_args()

ref = json.load(open(os.path.join(HERE, "cj_ref.json")))
fresh = ct.Solution(MECH)
fresh.TPX = ref["T1"], ref["p1"], ref["X"]
drv = ct.Solution(MECH)
drv.TPX = ref["T1"], ref["p1"], ref["X"]
drv.equilibrate("UV")
drv.TP = drv.T, args.pdrive * drv.P
drv.equilibrate("TP")

ind = ref["ind_len"]
Nx = int(math.ceil(args.L / (ind / args.ppi)))
dx = args.L / Nx
dt = args.cfl * dx / (ref["u_vN"] + ref["c_vN"])
xf = 40 * ind  # tz: front position at t0
tend = 0.9 * (args.L - (xf if args.ic == "tz" else 0.0)) / ref["D_CJ"]
NS = max(1, int(round(tend / dt / args.nsave))) if not args.steps else max(1, args.steps // args.nsave)
NT = args.steps or NS * args.nsave

case = {
    "run_time_info": "F",
    "x_domain%beg": 0.0,
    "x_domain%end": args.L,
    "m": Nx - 1,
    "n": 0,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": NS,
    "parallel_io": "T" if args.mfc.get("mpi", True) else "F",
    "model_eqns": 2,
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "weno_avg": "F",
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -2,
    "bc_x%end": -3,
    "chemistry": "T",
    "chem_params%diffusion": "F",
    "chem_params%reactions": "T",
    "chem_params%reaction_substeps": args.substeps,
    "chem_params%reaction_substeps_max": args.substeps_max,
    "chem_params%adap_substeps": "T",
    "cantera_file": MECH,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
}


def patch(k, x0, x1, T, p, u, Y, overwrite=True):
    """Constant-state slab [x0, x1] as patch k (later patches overwrite earlier ones)."""
    g = ct.Solution(MECH)
    g.TPY = T, p, Y
    case.update({f"patch_icpp({k})%geometry": 1, f"patch_icpp({k})%x_centroid": 0.5 * (x0 + x1), f"patch_icpp({k})%length_x": x1 - x0})
    case.update({f"patch_icpp({k})%vel(1)": u, f"patch_icpp({k})%pres": p, f"patch_icpp({k})%alpha(1)": 1.0, f"patch_icpp({k})%alpha_rho(1)": g.density})
    if k > 1:
        case[f"patch_icpp({k})%alter_patch(1)"] = "T"
    case.update({f"patch_icpp({k})%Y({i + 1})": float(y) for i, y in enumerate(g.Y)})


patch(1, 0.0, args.L, fresh.T, fresh.P, 0.0, fresh.Y)
if args.ic == "driver":
    patch(2, 0.0, args.driver_len * ind, drv.T, drv.P, 0.0, drv.Y)
    case["num_patches"] = 2
else:
    tz = np.loadtxt(os.path.join(HERE, "taylor_eq.txt"), skiprows=1)[::-1]  # p T rho u c xi, xi ascending
    t0 = xf / ref["D_CJ"]
    eq = ct.Solution(MECH)

    def eq_Y(T, p):
        eq.TPX = T, p, ref["X"]
        eq.equilibrate("TP")
        return eq.Y

    xi = np.linspace(tz[0, 5], ref["D_CJ"], 8)
    patch(2, 0.0, xi[0] * t0, tz[0, 1], tz[0, 0], 0.0, eq_Y(tz[0, 1], tz[0, 0]))
    for j in range(7):
        xm = 0.5 * (xi[j] + xi[j + 1])
        p_, T_, u_ = (np.interp(xm, tz[:, 5], tz[:, c]) for c in (0, 1, 3))
        patch(3 + j, xi[j] * t0, xi[j + 1] * t0, T_, p_, u_, eq_Y(T_, p_))
    patch(10, xf, xf + args.vnw * ind, ref["T_vN"], ref["p_vN"], ref["u_vN"], fresh.Y)
    case["num_patches"] = 10

if __name__ == "__main__":
    print(json.dumps(case))
