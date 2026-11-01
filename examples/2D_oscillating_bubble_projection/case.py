#!/usr/bin/env python3
"""
2D oscillating bubble with the all-Mach pressure projection: a validation of surface tension.

A gas bubble of radius R in a liquid starts at rest as an ellipse of semi-axes R(1 +/- eps), the n = 2 shape mode, and oscillates
about its circle under surface tension. For small eps the axis follows a(t) = R + eps R cos(w t) exp(-b t), with Lamb's frequency
of a 2D (cylindrical) interface, w^2 = n (n^2 - 1) sigma/((rho_g + rho_l) R^3), and the viscous decay of the irrotational mode in
the liquid outside, b = 2 n (n + 1) nu_l/R^2 (dissipation over twice the kinetic energy of the potential r^-n cos(n theta)); the
gas, 1000 times lighter and 100 times less viscous, adds little to either. One quadrant is simulated, with symmetry planes along
both axes and slip walls 4R from the center.

rho_l = R = sigma = 1. The liquid is a stiffened gas with sound speed 100 and the bubble an ideal gas (sound speed ~370) against
flow speeds of ~0.1, so the flow is incompressible to O(Mach^2): the projection steps at the capillary limit, while --explicit
(HLLC, conservative surface tension) must resolve the sound speed. The reference pressure p0 = 100 keeps the bubble's breathing
mode, omega_0^2 ~ 2 gamma p0/(rho_l R^2 ln(L/R)), near 14, well above the shape mode's 2.45, so the two do not exchange energy.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="2D oscillating bubble, all-Mach pressure projection")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead of the projection")
args, _ = parser.parse_known_args()

R, eps, L, ppr = 1.0, 0.05, 4.0, 32  # bubble radius, mode amplitude, domain size, cells per radius
rho_l, rho_g, mu_l, mu_g, sigma = 1.0, 1.0e-3, 1.0e-2, 1.0e-4, 1.0
gamma, gamma_g, c0, p0 = 4.4, 1.4, 100.0, 100.0  # liquid stiffened-gas exponent, gas ratio of heats, liquid sound speed
omega = math.sqrt(6 * sigma / ((rho_l + rho_g) * R**3))  # n = 2
periods = 5
t_stop = periods * 2 * math.pi / omega

pi_inf_l = c0**2 * rho_l / gamma - p0
eps_a = 1.0e-8

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "y_domain%beg": 0.0,
            "y_domain%end": L,
            "m": int(L * ppr) - 1,
            "n": int(L * ppr) - 1,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.4,
            "n_start": 0,
            "t_stop": t_stop,
            "t_save": t_stop / (20 * periods),
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
            "bc_x%beg": -2,
            "bc_x%end": -15,
            "bc_y%beg": -2,
            "bc_y%end": -15,
            **({} if args.explicit else {"proj_method": "T", "proj_max_acfl": 100.0}),
            "int_comp": 2,
            "ic_beta": 0.8,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1.0 / mu_l,
            "fluid_pp(2)%Re(1)": 1.0 / mu_g,
            "surface_tension": "T",
            "sigma": sigma,
            # The well-balanced model needs the projection's face pressure gradient
            "surface_tension_model": "conservative" if args.explicit else "well_balanced",
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Liquid
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%y_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": (1 - eps_a) * rho_l,
            "patch_icpp(1)%alpha_rho(2)": eps_a * rho_g,
            "patch_icpp(1)%alpha(1)": 1 - eps_a,
            "patch_icpp(1)%alpha(2)": eps_a,
            "patch_icpp(1)%cf_val": 0,
            # The elliptic bubble, centered on the domain corner, at the circle's Laplace pressure
            "patch_icpp(2)%geometry": 5,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%smoothen": "T",
            "patch_icpp(2)%smooth_patch_id": 1,
            "patch_icpp(2)%smooth_coeff": 0.95,
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%radii(1)": R * (1 + eps),
            "patch_icpp(2)%radii(2)": R * (1 - eps),
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": p0 + sigma / R,
            "patch_icpp(2)%alpha_rho(1)": eps_a * rho_l,
            "patch_icpp(2)%alpha_rho(2)": (1 - eps_a) * rho_g,
            "patch_icpp(2)%alpha(1)": eps_a,
            "patch_icpp(2)%alpha(2)": 1 - eps_a,
            "patch_icpp(2)%cf_val": 1,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": gamma * pi_inf_l / (gamma - 1.0),
            "fluid_pp(2)%eos": "ideal_gas",
            "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
        }
    )
)
