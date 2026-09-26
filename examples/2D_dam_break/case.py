#!/usr/bin/env python3
"""
2D dam break (Martin & Moyce, Phil. Trans. R. Soc. A 244:312-324, 1952, tables 2 and 6).

A water column a wide and 2a high (a = 2.25 in) collapses under gravity in an air-filled
box 5a wide and 3a high with free-slip walls, as in the n^2 = 2 experiments. Water and air
carry their viscosities and surface tension, with MTHINC interface compression. The surge
front z and column height eta are compared with the data as Z = z/a against
T = t*sqrt(2g/a) and H = eta/(2a) against t*sqrt(g/a); see analyze.py.

The run ends at T = 5 with adaptive time stepping. The all-Mach pressure projection steps
at the advective limit, a few hundred times the water acoustic one; --explicit runs the
HLLC solver, at the acoustic limit, as a control.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="2D dam break, all-Mach pressure projection")
parser.add_argument("--ppa", type=int, default=40, help="cells per column width a (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.25, help="cfl_target: advective with the projection, acoustic with --explicit (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead, as a control")
parser.add_argument("--st-model", default="well_balanced", choices=["conservative", "well_balanced"], help="surface_tension_model (default: %(default)s)")
args, _ = parser.parse_known_args()

a = 0.05715
g = 9.81
Lx, Ly = 5 * a, 3 * a
gamma_w, p_inf_w, rho_w = 4.4, 6.0e8, 1000.0
gamma_a, rho_a = 1.4, 1.0
p0 = 1.0e5
mu_w, mu_a, sigma = 1.0e-3, 1.8e-5, 0.0728
t_stop = 5.0 * math.sqrt(a / (2 * g))  # T = 5
saves = 100

# Hydrostatic pressure in each phase with the free surface at y = 2a
p_air = f"{p0} + {rho_a * g} * ({Ly} - y)"
p_water = f"{p0} + {rho_a * g * (Ly - 2 * a)} + {rho_w * g} * ({2 * a} - y)"

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
            "cfl_adap_dt": "T",
            "cfl_target": args.cfl,
            "n_start": 0,
            "t_stop": t_stop,
            "t_save": t_stop / saves,
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
            "int_comp": 2,
            "ic_beta": 1.0,
            # The well-balanced model needs the projection's face pressure gradient
            "surface_tension_model": "conservative" if args.explicit else args.st_model,
            "bf_y": "T",
            "g_y": -g,
            "k_y": 0.0,
            "w_y": 0.0,
            "p_y": 0.0,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: air filling the box
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * Lx,
            "patch_icpp(1)%y_centroid": 0.5 * Ly,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%length_y": Ly,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p_air,
            "patch_icpp(1)%alpha_rho(1)": 0.0,
            "patch_icpp(1)%alpha_rho(2)": rho_a,
            "patch_icpp(1)%alpha(1)": 0.0,
            "patch_icpp(1)%alpha(2)": 1.0,
            "patch_icpp(1)%cf_val": 0,
            # Patch 2: the water column in the lower left corner
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.5 * a,
            "patch_icpp(2)%y_centroid": a,
            "patch_icpp(2)%length_x": a,
            "patch_icpp(2)%length_y": 2 * a,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": p_water,
            "patch_icpp(2)%alpha_rho(1)": rho_w,
            "patch_icpp(2)%alpha_rho(2)": 0.0,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%alpha(2)": 0.0,
            "patch_icpp(2)%cf_val": 1,
            # Fluid parameters
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma_w - 1.0),
            "fluid_pp(1)%pi_inf": gamma_w * p_inf_w / (gamma_w - 1.0),
            "fluid_pp(2)%eos": "ideal_gas",
            "fluid_pp(2)%gamma": 1.0 / (gamma_a - 1.0),
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1.0 / mu_w,
            "fluid_pp(2)%Re(1)": 1.0 / mu_a,
            "surface_tension": "T",
            "sigma": sigma,
        }
    )
)
