#!/usr/bin/env python3
"""
2D rising bubble benchmark (Hysing et al., International Journal for Numerical Methods in Fluids 60:1259-1288, 2009).

A bubble of diameter 0.5 starts at rest at (0.5, 0.5) in a 1 x 2 box of a heavier fluid and
rises under gravity g = 0.98 until t = 3. The top and bottom walls are no-slip, the sides slip.

    test case   rho_1   rho_2   mu_1   mu_2   sigma   (1: surrounding fluid, 2: bubble)
        1       1000    100     10     1      24.5    Re = 35, Eo = 10: the bubble stays compact
        2       1000    1       10     0.1    1.96    Re = 35, Eo = 125: skirt and trailing filaments

The benchmark is incompressible. Both phases are stiffened gases with a water-like sound speed
c = 1500 at the reference pressure, so the flow Mach number is ~2e-4 and compressibility enters at
O(Mach^2). The all-Mach pressure projection steps at the flow's pace regardless; --explicit runs
the HLLC solver, which must resolve that sound speed. Compare with the reference data using
analyze.py.
"""

import argparse
import json

parser = argparse.ArgumentParser(description="2D rising bubble benchmark (Hysing et al. 2009)")
parser.add_argument("--case", type=int, default=1, choices=[1, 2], help="benchmark test case (default: %(default)s)")
parser.add_argument("--ppl", type=int, default=80, help="cells per unit length; the benchmark's h = 1/ppl (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.25, help="cfl_target: advective with the projection, acoustic with --explicit (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead, as a control")
args, _ = parser.parse_known_args()

rho_1, mu_1 = 1000.0, 10.0
rho_2, mu_2, sigma = {1: (100.0, 1.0, 24.5), 2: (1.0, 0.1, 1.96)}[args.case]
g = 0.98
Lx, Ly = 1.0, 2.0
xb, yb, R = 0.5, 0.5, 0.25
t_stop = 3.0
saves = 300

# Stiffened gases sharing the sound speed c0 at the reference pressure p0 (the top wall)
gamma, c0, p0 = 4.4, 1500.0, 1.0e5
pi_inf_1 = c0**2 * rho_1 / gamma - p0
pi_inf_2 = c0**2 * rho_2 / gamma - p0

# Hydrostatic pressure of the surrounding fluid; the bubble adds the Laplace jump sigma/R
p_1 = f"{p0} + {rho_1 * g} * ({Ly} - y)"
p_2 = f"{p0 + sigma / R} + {rho_1 * g} * ({Ly} - y)"
eps = 1.0e-8

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": Lx,
            "y_domain%beg": 0.0,
            "y_domain%end": Ly,
            "m": int(Lx * args.ppl) - 1,
            "n": int(Ly * args.ppl) - 1,
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
            "bc_y%beg": -16,
            "bc_y%end": -16,
            "proj_method": "F" if args.explicit else "T",
            # MTHINC, softened: sharper profiles let thin filaments drain cells past empty at the projection's CFL
            "int_comp": 2,
            "ic_beta": 1.0,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1.0 / mu_1,
            "fluid_pp(2)%Re(1)": 1.0 / mu_2,
            "surface_tension": "T",
            "sigma": sigma,
            # The well-balanced model needs the projection's face pressure gradient
            "surface_tension_model": "conservative" if args.explicit else "well_balanced",
            "bf_y": "T",
            "g_y": -g,
            "k_y": 0.0,
            "w_y": 0.0,
            "p_y": 0.0,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: the surrounding fluid
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * Lx,
            "patch_icpp(1)%y_centroid": 0.5 * Ly,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%length_y": Ly,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p_1,
            "patch_icpp(1)%alpha_rho(1)": (1 - eps) * rho_1,
            "patch_icpp(1)%alpha_rho(2)": eps * rho_2,
            "patch_icpp(1)%alpha(1)": 1 - eps,
            "patch_icpp(1)%alpha(2)": eps,
            "patch_icpp(1)%cf_val": 0,
            # Patch 2: the bubble, smoothed over a few cells
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%smoothen": "T",
            "patch_icpp(2)%smooth_patch_id": 1,
            "patch_icpp(2)%smooth_coeff": 0.95,
            "patch_icpp(2)%x_centroid": xb,
            "patch_icpp(2)%y_centroid": yb,
            "patch_icpp(2)%radius": R,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": p_2,
            "patch_icpp(2)%alpha_rho(1)": eps * rho_1,
            "patch_icpp(2)%alpha_rho(2)": (1 - eps) * rho_2,
            "patch_icpp(2)%alpha(1)": eps,
            "patch_icpp(2)%alpha(2)": 1 - eps,
            "patch_icpp(2)%cf_val": 1,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": gamma * pi_inf_1 / (gamma - 1.0),
            "fluid_pp(2)%eos": "stiffened_gas",
            "fluid_pp(2)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(2)%pi_inf": gamma * pi_inf_2 / (gamma - 1.0),
        }
    )
)
