#!/usr/bin/env python3
"""
2D rising bubble that expands as it rises, with the all-Mach pressure projection and the six-equation model.

A 4 mm gas bubble rises from rest through a 30 cm liquid column whose top is a pressure outlet at 1 atm. The liquid is a synthetic
viscous oil under ten-fold gravity: rho = 965 kg/m^3, mu = 0.23 Pa s, sigma = 0.30 N/m, g = 98.1 m/s^2 (Eotvos 5, Morton 0.01).

The six-equation model is solved with the frozen modulus, the default, which stays stable at large acoustic CFL; --wood uses Wood's
modulus instead (proj_woods_modulus). MTHINC (ic_beta 0.8) keeps the interface sharp: without it the bubble smears and trails gas.
"""

import argparse
import json
import sys

parser = argparse.ArgumentParser(description="2D rising bubble that expands as it rises (all-Mach projection, six-equation model)")
parser.add_argument("--ppd", type=int, default=32, help="cells per bubble diameter (default: %(default)s)")
parser.add_argument("--tstop", type=float, default=0.4, help="end time [s] (default: %(default)s)")
parser.add_argument("--saves", type=int, default=50, help="number of outputs (default: %(default)s)")
parser.add_argument("--five", action="store_true", help="five-equation (Allaire) model instead, as a control")
parser.add_argument("--wood", action="store_true", help="six-equation model solved with Wood's modulus instead of the frozen one")
args, _ = parser.parse_known_args()

# Bubble and column
d = 4.0e-3
R = 0.5 * d
W, H = 4 * d, 0.3  # column width and liquid depth
xc, yc = 0.5 * W, 1.5 * W  # bubble centre
dx = d / args.ppd

# Synthetic viscous liquid (water's stiffened-gas EOS) and air
rho_l, mu_l, sigma = 965.0, 0.23, 0.30
gamma_l, pi_inf_l = 6.12, 3.43e8
gamma_g, R_air, T, mu_g = 1.4, 287.0, 293.0, 1.8e-5
g, p_top = 98.1, 101325.0


def num(v):
    return f"{v:.7g}"


eps = 1.0e-8
w = 0.75 * dx  # half-width of the bubble's tanh interface
inside = f"0.5*(1 - tanh((sqrt((x - {num(xc)})**2 + (y - {num(yc)})**2) - {num(R)})/{num(w)}))"  # 1 in the bubble, 0 outside
alpha_g = f"({num(eps)} + {num(1 - 2 * eps)}*{inside})"
p_l = f"({num(p_top)} + {num(rho_l * g)}*({num(H)} - y))"  # hydrostatic pressure of the liquid
p_c = p_top + rho_l * g * (H - yc) + sigma / R  # the bubble's pressure: the liquid's at its centre plus the 2D Laplace jump
p_bubble = f"({p_l} + {inside}*({num(sigma / R)} + {num(rho_l * g)}*(y - {num(yc)})))"
rho_g = p_c / (R_air * T)

print(
    f"bubble d = {d * 1e3:.0f} mm, column {W * 100:.1f} x {H * 100:.0f} cm, {int(round(W / dx))} x {int(round(H / dx))} cells, "
    f"p_bottom/p_top = {(p_top + rho_l * g * H) / p_top:.3f}, Eo = {rho_l * g * d**2 / sigma:.2f}, Mo = {g * mu_l**4 / (rho_l * sigma**3):.3g}",
    file=sys.stderr,
)

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": W,
            "y_domain%beg": 0.0,
            "y_domain%end": H,
            "m": int(round(W / dx)) - 1,
            "n": int(round(H / dx)) - 1,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.25,
            "n_start": 0,
            "t_stop": args.tstop,
            "t_save": args.tstop / args.saves,
            "num_patches": 2,
            "model_eqns": 2 if args.five else 3,
            "num_fluids": 2,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "proj_method": "T",
            **({"proj_woods_modulus": "T"} if args.wood and not args.five else {}),
            # Slip side walls and a no-slip floor; the top is a pressure outlet whose ghost, half a cell above the surface,
            # holds p_top less half a cell of hydrostatic head
            "bc_x%beg": -15,
            "bc_x%end": -15,
            "bc_y%beg": -16,
            "bc_y%end": -3,
            "bc_y%pres_out": p_top - 0.5 * rho_l * g * dx,
            "int_comp": 2,
            "ic_beta": 0.8,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1.0 / mu_l,
            "fluid_pp(2)%Re(1)": 1.0 / mu_g,
            "surface_tension": "T",
            "sigma": sigma,
            "surface_tension_model": "well_balanced",
            "bf_y": "T",
            "g_y": -g,
            "k_y": 0.0,
            "w_y": 0.0,
            "p_y": 0.0,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: the liquid, hydrostatic
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * W,
            "patch_icpp(1)%y_centroid": 0.5 * H,
            "patch_icpp(1)%length_x": W,
            "patch_icpp(1)%length_y": H,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p_l,
            "patch_icpp(1)%alpha_rho(1)": (1 - eps) * rho_l,
            "patch_icpp(1)%alpha_rho(2)": eps * p_top / (R_air * T),
            "patch_icpp(1)%alpha(1)": 1 - eps,
            "patch_icpp(1)%alpha(2)": eps,
            "patch_icpp(1)%cf_val": 0,
            # Patch 2: the bubble, an analytic patch over the window around it
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.5 * W,
            "patch_icpp(2)%y_centroid": yc,
            "patch_icpp(2)%length_x": W,
            "patch_icpp(2)%length_y": W,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": p_bubble,
            "patch_icpp(2)%alpha_rho(1)": f"(1 - {alpha_g})*{num(rho_l)}",
            "patch_icpp(2)%alpha_rho(2)": f"{alpha_g}*{num(rho_g)}",
            "patch_icpp(2)%alpha(1)": f"1 - {alpha_g}",
            "patch_icpp(2)%alpha(2)": alpha_g,
            "patch_icpp(2)%cf_val": inside,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),
            "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),
            "fluid_pp(2)%eos": "stiffened_gas",
            "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
        }
    )
)
