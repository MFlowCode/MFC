#!/usr/bin/env python3
"""
Fixed random array of no-slip isothermal cylinders (2D) or spheres (3D) in a periodic box, held at a bulk mass flux.

The box has side L = (N pi D^2 / (4 phi))^(1/2) (2D) or (N pi D^3 / (6 phi))^(1/3) (3D) and holds N bodies at solid volume
fraction phi, placed by seeded random sequential addition with a minimum surface gap. Bodies may straddle the periodic
faces. The flow is driven by the mean pressure gradient controller (mean_mflux, mean_dpdx_gain), which holds the fluid-averaged
mass flux at Q = Re_m mu / ((1 - phi) D), so that Re_m = rho (1 - phi) W D / mu with W the fluid-phase mean velocity. The
force on every body goes to D/ib_forces.dat, the heat it takes to D/ib_heat.dat, and the controller state to mean_dpdx.dat.
Dry air from h2o2.yaml (chemistry on for thermodynamics and mixture-averaged transport, reactions off); the bodies are
isothermal at the free-stream temperature (--tw_ratio 1) or at --tw_ratio times it. --thermostat adds the thermally
developed controller of Sun et al. (2015): mean_temp_wall at the body temperature, mixing-cup temperature held at T_inf.

    ./mfc.sh run examples/2D_ibm_periodic_array/case.py -n 2 --gpu
    python3 examples/2D_ibm_periodic_array/check_array.py <case directory>

Defaults are the CI-sized run (12 cylinders, phi 0.2, Re_m 20, 14 cells per diameter, 96 x 96 cells, 200 fixed steps: the
flow is far from steady). A drag run needs --tstar 40 or more. --dim 3 gives the cubic box of spheres.
Production-like values: Re_m 20 to 100, phi 0.1 to 0.3, 24 cells per diameter, 80 bodies.
"""

import argparse
import json
import math
import random
import sys

import cantera as ct


def make_parser(dim=2):
    parser = argparse.ArgumentParser(prog=f"{dim}D_ibm_periodic_array")
    parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
    parser.add_argument("--dim", type=int, default=dim, choices=(2, 3), help="Space dimension.")
    parser.add_argument("--phi", type=float, default=0.2, help="Solid volume fraction.")
    parser.add_argument("--Re", type=float, default=20.0, help="Re_m = rho (1 - phi) W D / mu.")
    parser.add_argument("--mach", type=float, default=0.2, help="Mach number of W at the free-stream sound speed.")
    parser.add_argument("--n_bodies", type=int, default=12 if dim == 2 else 8, help="Number of bodies.")
    parser.add_argument("--seed", type=int, default=1, help="Seed of the random placement.")
    parser.add_argument("--gap", type=float, default=0.1, help="Minimum surface gap [D].")
    parser.add_argument("--cells_per_D", type=float, default=14.0, help="Cells across a body diameter.")
    parser.add_argument("--steps", type=int, default=200, help="Fixed number of steps.")
    parser.add_argument("--tstar", type=float, default=0.0, help="Run length t W / D (sets the number of steps when positive).")
    parser.add_argument("--saves", type=int, default=1, help="Field saves.")
    parser.add_argument("--cfl", type=float, default=0.3, help="Acoustic CFL number of the fixed step.")
    parser.add_argument("--gain", type=float, default=0.05, help="mean_dpdx_gain of the flux controller.")
    parser.add_argument("--stride", type=int, default=10, help="ib_force_stride: steps between force and heat rows.")
    parser.add_argument("--tw_ratio", type=float, default=1.0, help="Body temperature over T_inf.")
    parser.add_argument("--thermostat", action="store_true", help="mean_temp_wall controller (needs --tw_ratio != 1).")
    parser.add_argument("--ib_ip_interp", type=int, default=1, choices=(1, 2), help="Image-point interpolation (1: inverse distance).")
    parser.add_argument("--wall_order", type=int, default=2, choices=(1, 2), help="ib_wall_order (1 wall value, 2 quadratic, the solver default).")
    parser.add_argument("--visc_compact", type=int, default=1, choices=(0, 1), help="ib_visc_compact (1 compact viscous gradient next to the body, the solver default; 0 the WENO gradient).")
    return parser


def pack(n, lo, hi, periodic, d_min, seed):
    """n centres in the box [lo, hi] (periodic along the axes in `periodic`, which use the minimum-image distance) with
    separation >= d_min, by seeded random sequential addition."""
    rng = random.Random(seed)
    span = [h - l for l, h in zip(lo, hi)]

    def dist(c, p):
        return math.sqrt(sum((min(abs(x - y), s - abs(x - y)) if k in periodic else abs(x - y)) ** 2 for k, (x, y, s) in enumerate(zip(c, p, span))))

    pts = []
    for _ in range(100000 * n):
        c = [l + rng.random() * s for l, s in zip(lo, span)]
        if all(dist(c, p) >= d_min for p in pts):
            pts.append(c)
            if len(pts) == n:
                return pts
    sys.exit(f"pack placed only {len(pts)} of {n} bodies: lower phi or --gap")


def build(args):
    dim, T_inf, p_inf = args.dim, 300.0, 101325.0
    gas = ct.Solution("h2o2.yaml")
    gas.TPX = T_inf, p_inf, "O2:0.21,N2:0.79"
    rho, mu, a, gamma = gas.density, gas.viscosity, gas.sound_speed, gas.cp / gas.cv
    W = args.mach * a
    D = args.Re * mu / (rho * (1.0 - args.phi) * W)  # D follows from Re_m at fixed W
    n = args.n_bodies
    L = (n * math.pi / (4.0 * args.phi)) ** 0.5 * D if dim == 2 else (n * math.pi / (6.0 * args.phi)) ** (1.0 / 3.0) * D
    phi = n * math.pi / 4.0 * (D / L) ** 2 if dim == 2 else n * math.pi / 6.0 * (D / L) ** 3
    cells = round(L / D * args.cells_per_D)
    dx = L / cells
    dt = args.cfl * dx / (W + a)
    steps = round(args.tstar * D / (W * dt)) if args.tstar > 0 else args.steps
    centers = pack(n, [0.0] * dim, [L] * dim, range(dim), (1.0 + args.gap) * D, args.seed)
    Q = args.Re * mu / ((1.0 - phi) * D)
    case = {
        "run_time_info": "T",
        "cyl_coord": "F",
        "m": cells - 1,
        "n": cells - 1,
        "p": cells - 1 if dim == 3 else 0,
        "dt": dt,
        "t_step_start": 0,
        "t_step_stop": steps,
        "t_step_save": max(1, steps // args.saves),
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
        "weno_Re_flux": "T",
        "weno_avg": "T",
        "avg_state": 2,
        "riemann_solver": 2,
        "wave_speeds": 1,
        "viscous": "T",
        "fd_order": 4,
        # Chemistry: thermodynamics and mixture-averaged transport only
        "chemistry": "T",
        "chem_params%diffusion": "T",
        "chem_params%reactions": "F",
        "chem_params%transport_model": 1,
        "cantera_file": "h2o2.yaml",
        "chem_wrt_T": "T",
        "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
        "fluid_pp(1)%eos": "stiffened_gas",
        "fluid_pp(1)%pi_inf": 0.0,
        "fluid_pp(1)%Re(1)": 1.0 / mu,
        # Flow at W everywhere; the controller holds the fluid-averaged mass flux at Q
        "patch_icpp(1)%geometry": 3 if dim == 2 else 9,
        "patch_icpp(1)%pres": p_inf,
        "patch_icpp(1)%alpha_rho(1)": rho,
        "patch_icpp(1)%alpha(1)": 1.0,
        "mean_mflux(1)": Q,
        "mean_dpdx_gain": args.gain,
        # Fixed no-slip isothermal bodies; mass 0 keeps a body-force term out of the force record
        "ib": "T",
        "num_ibs": n,
        "many_ib_patch_parallelism": "T",
        "ib_ip_interp": args.ib_ip_interp,
        "ib_wall_order": args.wall_order,
        "ib_visc_compact": "T" if args.visc_compact else "F",
        "ib_state_wrt": "T",
        "ib_force_wrt": "T",
        "ib_heat_wrt": "T",
        "ib_force_stride": args.stride,
        "parallel_io": "F",
        "format": "silo",
        "precision": "double",
        "prim_vars_wrt": "F",
    }
    for k, name in enumerate("xyz"[:dim]):
        case.update({f"{name}_domain%beg": 0.0, f"{name}_domain%end": L, f"bc_{name}%beg": -1, f"bc_{name}%end": -1})
        case.update({f"patch_icpp(1)%{name}_centroid": 0.5 * L, f"patch_icpp(1)%length_{name}": L, f"patch_icpp(1)%vel({k + 1})": W if k == 0 else 0.0})
    for i, c in enumerate(centers):
        for key, val in (("geometry", 2 if dim == 2 else 8), ("radius", 0.5 * D), ("slip", "F"), ("thermal_bc", 1), ("Twall", args.tw_ratio * T_inf), ("mass", 0.0)):
            case[f"patch_ib({i + 1})%{key}"] = val
        for name, x in zip("xyz", c):
            case[f"patch_ib({i + 1})%{name}_centroid"] = x
    for i, y in enumerate(gas.Y):
        case[f"patch_icpp(1)%Y({i + 1})"] = float(y)
        case[f"chem_wrt_Y({i + 1})"] = "F"
    # Initial pressure gradient: a Stokes-scale estimate of the array drag; the controller replaces it within a few steps
    case["mean_dpdx(1)"] = -32.0 * mu * W * phi / (1.0 - phi) / D**2
    if args.thermostat:
        k_gas = gas.thermal_conductivity
        case.update({"mean_temp": T_inf, "mean_temp_gain": args.gain, "mean_temp_wall": args.tw_ratio * T_inf, "mean_qdot": -k_gas * (args.tw_ratio - 1.0) * T_inf / D**2})
        case["chem_params%diffusion"] = "T"
    info = f"{dim}D  N {n}  phi {phi:.4f}  L {L / D:.3f} D  Re_m {args.Re:g}  Mach {args.mach:g}  D {D:.4e} m  W {W:.2f} m/s  grid {cells}^{dim}  Q {Q:.5f}"
    return case, info


if __name__ == "__main__":
    args = make_parser().parse_known_args()[0]
    case, info = build(args)
    print(info, file=sys.stderr)
    print(json.dumps(case))
