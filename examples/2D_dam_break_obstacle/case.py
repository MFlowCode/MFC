#!/usr/bin/env python3
"""
2D dam break against an obstacle (Koshizuka, Tamako & Oka, Computational Fluid Dynamics Journal 4:29-46, 1995).

A water column a = 0.146 m wide and 2a high collapses in a 4a x 4a air-filled tank with
free-slip walls and hits a no-slip block 0.024 m wide and 0.048 m high centred on the floor
at x = 2a, which throws the surge up into a jet that falls back onto the far side. The block
is a stationary immersed boundary. Water and air carry their viscosities and surface tension,
with MTHINC interface compression, and the run ends at t = 0.5 s with adaptive stepping.
"""

import argparse
import json

parser = argparse.ArgumentParser(description="2D dam break against an obstacle, all-Mach pressure projection")
parser.add_argument("--ppa", type=int, default=40, help="cells per column width a (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.25, help="cfl_target: advective with the projection, acoustic with --explicit (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead, as a control")
parser.add_argument("--st-model", default="well_balanced", choices=["conservative", "well_balanced"], help="surface_tension_model (default: %(default)s)")
args, _ = parser.parse_known_args()

a = 0.146
g = 9.81
L = 4 * a
w_obs, h_obs = 0.024, 0.048
gamma_w, p_inf_w, rho_w = 4.4, 6.0e8, 1000.0
gamma_a, rho_a = 1.4, 1.0
p0 = 1.0e5
mu_w, mu_a, sigma = 1.0e-3, 1.8e-5, 0.0728
t_stop = 0.6
saves = 120

# Hydrostatic pressure in each phase with the free surface at y = 2a
p_air = f"{p0} + {rho_a * g} * ({L} - y)"
p_water = f"{p0} + {rho_a * g * (L - 2 * a)} + {rho_w * g} * ({2 * a} - y)"

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "y_domain%beg": 0.0,
            "y_domain%end": L,
            "m": 4 * args.ppa - 1,
            "n": 4 * args.ppa - 1,
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
            # MTHINC, softened: sharper profiles (1.0 and THINC's 1.6) let the thin splash jets drain cells past empty at the
            # projection's CFL
            "int_comp": 2,
            "ic_beta": 0.6,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1.0 / mu_w,
            "fluid_pp(2)%Re(1)": 1.0 / mu_a,
            "surface_tension": "T",
            "sigma": sigma,
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
            # Patch 1: air filling the tank
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%y_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p_air,
            "patch_icpp(1)%alpha_rho(1)": 0.0,
            "patch_icpp(1)%alpha_rho(2)": rho_a,
            "patch_icpp(1)%alpha(1)": 0.0,
            "patch_icpp(1)%alpha(2)": 1.0,
            "patch_icpp(1)%cf_val": 0,
            # Patch 2: the water column against the left wall
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
            # The obstacle: a no-slip block on the floor at mid-tank, extended below the floor so it meets the wall cleanly
            "ib": "T",
            "num_ibs": 1,
            "fd_order": 2,
            "patch_ib(1)%geometry": 3,
            "patch_ib(1)%x_centroid": 2 * a,
            "patch_ib(1)%y_centroid": 0.5 * (h_obs - a),
            "patch_ib(1)%length_x": w_obs,
            "patch_ib(1)%length_y": h_obs + a,
            "patch_ib(1)%slip": "F",
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma_w - 1.0),
            "fluid_pp(1)%pi_inf": gamma_w * p_inf_w / (gamma_w - 1.0),
            "fluid_pp(2)%eos": "ideal_gas",
            "fluid_pp(2)%gamma": 1.0 / (gamma_a - 1.0),
        }
    )
)
