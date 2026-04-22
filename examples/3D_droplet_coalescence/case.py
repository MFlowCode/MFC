#!/usr/bin/env python3

import argparse
import json
import math

# Parameter table: (We, Re, B)
CASES = {
    "a": (0.2, 14.8, 0.20),
    "b": (0.5, 23.6, 0.10),
    "c": (8.6, 105.9, 0.08),
    "d": (15.2, 139.8, 0.08),
    "e": (19.4, 158.0, 0.05),
    "f": (32.8, 210.8, 0.08),
    "g": (37.2, 228.0, 0.01),
    "h": (61.4, 296.5, 0.06),
    "i": (61.3, 295.3, 0.11),
    "j": (56.3, 288.9, 0.13),
    "k": (70.8, 327.7, 0.25),
    "l": (48.1, 270.1, 0.39),
    "m": (60.1, 302.8, 0.55),
    "n": (65.1, 320.3, 0.49),
    "o": (60.8, 313.7, 0.68),
    "p": (64.9, 312.8, 0.71),
    "q": (48.8, 260.3, 0.72),
    "r": (14.5, 149.1, 0.34),
}

parser = argparse.ArgumentParser()
parser.add_argument("--mfc", type=json.loads, default="{}")
parser.add_argument("--case", type=str, default="a", choices=CASES.keys())
args = parser.parse_args()

We, Re, B = CASES[args.case]

# Physical properties (tetradecane / nitrogen) — SI units
rho_l = 763.0  # kg/m^3
rho_g = rho_l / 666  # kg/m^3
sigma = 0.0266  # N/m
gamma_l = 2.35
gamma_g = 1.4

D = 240e-6  # m
R = D / 2.0

# Pressure and EOS parameters scaled from non-dimensional formulation
# (low effective pressure keeps the speed of sound manageable)
p_ref = sigma / (0.72 * D)
p_gas = p_ref
pi_inf_l = 190.0 * p_ref
pi_inf_g = 0.0
p_liq = p_gas + 2 * sigma / R

# Collision parameters from dimensionless numbers
Ur = math.sqrt(We * sigma / (rho_l * D))
Ud = Ur / 2.0

mu_l = rho_l * Ur * D / Re
mu_g = mu_l / 119

# Domain (m)
x0, x1 = -2.25 * D, 2.25 * D
y0, y1 = -1.5 * D, 1.5 * D
z0, z1 = -1.5 * D, 1.5 * D

Nx = 399
Ny = 279
Nz = 279

eps = 1e-9

sep = 0.55 * D
b = B * D

# Simulation time
t_end = 1e-3  # seconds
t_save = t_end / 1000


data = {
    "run_time_info": "T",
    "x_domain%beg": x0,
    "x_domain%end": x1,
    "y_domain%beg": y0,
    "y_domain%end": y1,
    "z_domain%beg": z0,
    "z_domain%end": z1,
    "m": Nx,
    "n": Ny,
    "p": Nz,
    "cyl_coord": "F",
    "n_start": 0,
    # adaptive timestepping (disabled)
    # "cfl_adap_dt": "T",
    # "cfl_target": 0.1,
    # "t_stop": t_end,
    # "t_save": t_save,
    # fixed timestepping
    "dt": 1e-7,
    "t_step_start": 0,
    "t_step_stop": 1,
    "t_step_save": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "low_Mach": 0,
    "mixture_err": "T",
    "mpp_lim": "T",
    "time_stepper": 3,
    "recon_type": 2,
    "muscl_order": 2,
    "muscl_lim": 4,
    "int_comp": "T",
    "ic_beta": 1.6,
    "avg_state": 2,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "viscous": "T",
    "bc_x%beg": -6,
    "bc_x%end": -6,
    "bc_y%beg": -6,
    "bc_y%end": -6,
    "bc_z%beg": -6,
    "bc_z%end": -6,
    "num_patches": 3,
    "num_fluids": 2,
    "surface_tension": "T",
    "format": 1,
    "precision": 2,
    "alpha_wrt(1)": "T",
    "cf_wrt": "T",
    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",
    "vel_wrt(3)": "T",
    "pres_wrt": "T",
    "parallel_io": "T",
    "sigma": sigma,
    "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),
    "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),
    "fluid_pp(1)%Re(1)": 1.0 / mu_l,
    "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(2)%Re(1)": 1.0 / mu_g,
    # Patch 1: Gas background
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": x1 - x0,
    "patch_icpp(1)%length_y": y1 - y0,
    "patch_icpp(1)%length_z": z1 - z0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p_gas,
    "patch_icpp(1)%alpha_rho(1)": eps * rho_l,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_g,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    "patch_icpp(1)%cf_val": 0,
    # Patch 2: Left droplet
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%x_centroid": -sep,
    "patch_icpp(2)%y_centroid": -b / 2.0,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": R,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.95,
    "patch_icpp(2)%vel(1)": Ud,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p_liq,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(2)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "patch_icpp(2)%cf_val": 1,
    # Patch 3: Right droplet
    "patch_icpp(3)%geometry": 8,
    "patch_icpp(3)%x_centroid": sep,
    "patch_icpp(3)%y_centroid": b / 2.0,
    "patch_icpp(3)%z_centroid": 0.0,
    "patch_icpp(3)%radius": R,
    "patch_icpp(3)%alter_patch(1)": "T",
    "patch_icpp(3)%smoothen": "T",
    "patch_icpp(3)%smooth_patch_id": 1,
    "patch_icpp(3)%smooth_coeff": 0.95,
    "patch_icpp(3)%vel(1)": -Ud,
    "patch_icpp(3)%vel(2)": 0.0,
    "patch_icpp(3)%vel(3)": 0.0,
    "patch_icpp(3)%pres": p_liq,
    "patch_icpp(3)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(3)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(3)%alpha(1)": 1.0 - eps,
    "patch_icpp(3)%alpha(2)": eps,
    "patch_icpp(3)%cf_val": 1,
}

print(json.dumps(data))
