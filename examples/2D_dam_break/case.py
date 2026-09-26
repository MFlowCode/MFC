#!/usr/bin/env python3
"""
2D dam break (Martin & Moyce, Phil. Trans. R. Soc. A 244:312-324, 1952, tables 2 and 6).

A water column a wide and 2a high (a = 2.25 in) collapses under gravity in an air-filled
box 5a wide and 3a high with free-slip walls, as in the n^2 = 2 experiments. The surge
front z and column height eta are compared with the data as Z = z/a against
T = t*sqrt(2g/a) and H = eta/(2a) against t*sqrt(g/a); see analyze.py.

The all-Mach pressure projection steps at the advective limit, a few hundred times the
water acoustic limit (reported at startup); --explicit runs the HLLC solver at the
acoustic limit as a control. --hydrostatic spans the water across the box, a rest state
the discrete scheme must hold.
"""

import argparse
import json
import math
import sys

parser = argparse.ArgumentParser(description="2D dam break, all-Mach pressure projection")
parser.add_argument("--ppa", type=int, default=40, help="cells per column width a (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.25, help="advective CFL based on sqrt(2 g 2a) (default: %(default)s)")
parser.add_argument("--T", type=float, default=3.4, help="final time in units of sqrt(a/2g) (default: %(default)s)")
parser.add_argument("--saves", type=int, default=34, help="number of output snapshots (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead, as a control")
parser.add_argument("--hydrostatic", action="store_true", help="water layer across the whole box: must stay at rest")
args, _ = parser.parse_known_args()

a = 0.05715
g = 9.81
Lx, Ly = 5 * a, 3 * a
gamma_w, p_inf_w, rho_w = 4.4, 6.0e8, 1000.0
gamma_a, rho_a = 1.4, 1.0
p0 = 1.0e5

dx = a / args.ppa
c_w = math.sqrt(gamma_w * (p0 + p_inf_w) / rho_w)
dt = args.cfl * dx / (c_w if args.explicit else math.sqrt(2 * g * 2 * a))
Nt = int(math.ceil(args.T * math.sqrt(a / (2 * g)) / dt))
Ns = max(Nt // args.saves, 1)
print(f"dt = {dt:.3e} s = {dt * c_w / dx:.3g} x water acoustic limit, {Nt} steps", file=sys.stderr)

# Hydrostatic pressure in each phase with the free surface at y = 2a
p_air = f"{p0} + {rho_a * g} * ({Ly} - y)"
p_water = f"{p0} + {rho_a * g * (Ly - 2 * a)} + {rho_w * g} * ({2 * a} - y)"


def patch(i, xc, yc, lx, ly, water, pres):
    return {
        f"patch_icpp({i})%geometry": 3,
        f"patch_icpp({i})%x_centroid": xc,
        f"patch_icpp({i})%y_centroid": yc,
        f"patch_icpp({i})%length_x": lx,
        f"patch_icpp({i})%length_y": ly,
        f"patch_icpp({i})%vel(1)": 0.0,
        f"patch_icpp({i})%vel(2)": 0.0,
        f"patch_icpp({i})%pres": pres,
        f"patch_icpp({i})%alpha_rho(1)": rho_w if water else 0.0,
        f"patch_icpp({i})%alpha_rho(2)": 0.0 if water else rho_a,
        f"patch_icpp({i})%alpha(1)": 1.0 if water else 0.0,
        f"patch_icpp({i})%alpha(2)": 0.0 if water else 1.0,
    }


wx = Lx if args.hydrostatic else a
water = patch(2, 0.5 * wx, a, wx, 2 * a, True, p_water)
water["patch_icpp(2)%alter_patch(1)"] = "T"

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": Lx,
            "y_domain%beg": 0.0,
            "y_domain%end": Ly,
            "m": 5 * args.ppa - 1,
            "n": 3 * args.ppa - 1,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": Ns,
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 2,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -15,
            "bc_x%end": -15,
            "bc_y%beg": -15,
            "bc_y%end": -15,
            "proj_method": "F" if args.explicit else "T",
            "bf_y": "T",
            "g_y": -g,
            "k_y": 0.0,
            "w_y": 0.0,
            "p_y": 0.0,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            **patch(1, 0.5 * Lx, 0.5 * Ly, Lx, Ly, False, p_air),
            **water,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma_w - 1.0),
            "fluid_pp(1)%pi_inf": gamma_w * p_inf_w / (gamma_w - 1.0),
            "fluid_pp(2)%eos": "ideal_gas",
            "fluid_pp(2)%gamma": 1.0 / (gamma_a - 1.0),
        }
    )
)
