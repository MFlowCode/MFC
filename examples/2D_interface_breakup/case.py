#!/usr/bin/env python3
# 2D water-gas interfacial-breakup analogue (multimode Rayleigh-Taylor / Faraday) for MFC.
#
# Runs with the all-Mach pressure projection by default: the acoustics are implicit, so the
# step is set by the flow (and gravity, capillarity, viscosity) rather than the sound speed,
# and water keeps its physical stiffness. --explicit runs the explicit HLLC solver instead,
# with the settings of the original acoustic-CFL case (softened water, mpp_lim, low-Mach fix).
#
# 2D variant of caseWater3d.py. The heavy fluid is liquid WATER (stiffened-gas EOS,
# nearly incompressible); the light fluid is a GAS (air, ideal gas). A single-valued
# interface y = y_mean + eta(x) is seeded with a multimode statistical perturbation
# (perturbation2d.py) and driven by a background + oscillatory acceleration through
# MFC's body force:
#
#     accel_y(t) = g_y + k_y * sin(w_y * t - p_y)
#
# The perturbed interface is injected by hardcoded IC case 209 (2dHardcodedIC.fpp),
# which reads interface_profile.dat (written below; two columns x, y_int) and
# tanh-blends the two fluids' constant states across y = y_int(x).
#
# Physical anchoring (same as caseWater3d.py)
# The Weber number is defined at the scale surface tension actually acts on: the
# SHORTEST SEEDED WAVELENGTH lambda_min, with the deep-water gravity-wave speed at
# that wavelength as the velocity scale:
#
#     c_gw^2 = g lambda_min / (2 pi)                       gravity-wave speed at lambda_min
#     We     = rho_w c_gw^2 lambda_min / sigma
#            = rho_w g lambda_min^2 / (2 pi sigma)
#
# This is exactly the ratio of the gravity to capillary restoring terms in the
# gravity-capillary dispersion relation c^2 = g/k + sigma k / rho at k = 2 pi / lambda_min:
#     We > 1  ->  the smallest seeded modes are gravity waves (capillarity is a perturbation)
#     We < 1  ->  lambda_min sits below the capillary length; those modes are capillary waves
#
# For a REAL fluid pair the dimensionless numbers are NOT independent. The band ratio
# lambda_min / H is fixed (1/16, section 6), so We still sets the length scale:
# lambda_min = sqrt(2 pi We sigma / (rho_w g)), H = lambda_min / (lambda_min/H), and
#
#     U_c = sqrt(g H)                                      reference velocity
#     At  = (rho_w - rho_g) / (rho_w + rho_g)              ~0.998 for water-air
#     Re  = rho_w U_c H / mu_w  = rho_w sqrt(g) H^1.5 / mu_w    ~ We^(3/4)
#
# so fixing the fluid pair + g + We pins Re: We is the control knob, Re comes out DERIVED.
#
# Relation to the old H-based definition (We_H = rho_w g H^2 / sigma):
#     We_H = 2 pi (H / lambda_min)^2 We = 512 pi We     (for lambda_min = H/16)
# so the legacy We_H = 500 case corresponds to We ~= 0.311 (the default).
#
# Nondimensionalization (code units): rho_w = H = g0 = 1
#   length    H   = mean interface height
#   density   rho_w
#   accel     g0  = |background acceleration|
#   time      t_c = sqrt(H / g0)            = 1
#   velocity  U_c = sqrt(g0 H)              = 1
#   pressure  p_c = rho_w U_c^2 = rho_w g0 H
#   We = rho_w g0 lambda_min^2 / (2 pi sigma)  ->  sigma_code = lambda_min^2 / (2 pi We)
#   Re = rho_w U_c   H / mu     ->  fluid Re(i) = rho_w U_c H / mu_i  (per fluid)
#
# Stiffened-gas water + artificial stiffness reduction
# Water is modeled as a stiffened gas (gamma_w ~ 4.4, pi_inf ~ 6e8 Pa) so it is nearly
# incompressible and its sound speed exceeds the gas's -- the physically correct regime.
# However, with an 833:1 density ratio the GAS sound speed is already pinned high
# (its density is tiny, yet p_ref must stay ~100x the body-force pressure swing to keep
# pressure positive), so it -- not water -- would otherwise set the acoustic timestep.
# The *physical* water pi_inf makes c_water ~ 9x the gas sound speed, shrinking dt by
# ~9x for no dynamical gain (water is already at Mach ~3e-4).
#
# We therefore default to matching c_water to c_gas (water_pi_inf_mode="match_gas"):
# water stays a stiffened gas at Mach ~3e-3 (still incompressible -- density varies
# O(Mach^2) ~ 1e-5), but the timestep is as large as pressure-positivity allows. Set
# water_pi_inf_mode="physical" to use the true pi_inf=6e8 Pa (correct, but ~9x costlier).

import argparse
import json
import math
import os
import sys

# perturbation2d.py sits next to this file, and pre_process reads interface_profile.dat from the
# case directory, so both are resolved relative to this file rather than the working directory
case_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, case_dir)

from perturbation2d import (
    cells_per_lambda,
    generate_perturbation_2d,
    print_diagnostics,
    save_diagnostics,
    write_interface_2d,
)

# 0. Command-line arguments (stdout must stay pure JSON for MFC, so argparse
#    help/errors go to stderr as usual and do not pollute the JSON on success)
parser = argparse.ArgumentParser(
    description="Generate a 2D MFC case for water-gas interfacial breakup.",
)
parser.add_argument("--forcing-freq", type=float, default=40.0, help="w_y * t_c: dimensionless angular forcing frequency " "(default: %(default)s)")
parser.add_argument("--forcing-ratio", type=float, default=40.0, help="k_y / g0: oscillation amplitude / background accel; " "0 => pure RT (default: %(default)s)")
parser.add_argument(
    "--We",
    type=float,
    default=0.311,
    help="Weber number rho_w c_gw^2 lambda_min / sigma at the shortest "
    "seeded wavelength (c_gw = gravity-wave speed at lambda_min); "
    "sets the length scale. Legacy H-based We_H = 512 pi We. "
    "(default: %(default)s, equiv. We_H ~= 500)",
)
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead of the pressure projection")
# MFC passes '--mfc' with a JSON payload when invoking case files; accept and
# ignore any unknown args so this stays compatible with the toolchain.
args, _ = parser.parse_known_args()

# 1. Physical fluid properties (SI) and the dimensionless control number
# Heavy fluid: liquid water.  Light fluid: air (gas).
rho_w_SI = 1000.0  # water density            [kg/m^3]
rho_g_SI = 1.2  # gas (air) density        [kg/m^3]
mu_w_SI = 1.0e-3  # water dynamic viscosity  [Pa.s]
mu_g_SI = 1.8e-5  # gas dynamic viscosity    [Pa.s]
sigma_SI = 0.072  # water-air surface tension[N/m]
g_SI = 9.81  # gravitational accel      [m/s^2]
p_atm_SI = 1.013e5  # atmospheric pressure     [Pa]   (reference / reporting only)

gam_w = 6.12  # water stiffened-gas exponent
pi_inf_w_SI = 3.43e8  # water stiffening pressure [Pa]
gam_g = 1.4  # gas ideal-gas exponent (pi_inf = 0)

We = args.We  # Weber number at lambda_min (gravity-wave speed) -> sets the length scale

# Seeded-band ratio: shortest seeded wavelength as a fraction of H. Defined here
# (rather than section 6) because the We definition is anchored to lambda_min,
# so it enters the length-scale inversion below. Section 6 reuses it.
lambda_min_over_H = 1.0 / 16.0

# 2. Reference scales + derived dimensionless groups
# We = rho_w g lambda_min^2 / (2 pi sigma)
#   -> lambda_min = sqrt(2 pi We sigma / (rho_w g)),  H = lambda_min / (lambda_min/H)
lambda_min_SI = math.sqrt(2.0 * math.pi * We * sigma_SI / (rho_w_SI * g_SI))
H_SI = lambda_min_SI / lambda_min_over_H  # mean interface height [m]
U_SI = math.sqrt(g_SI * H_SI)  # reference velocity   [m/s]
t_SI = math.sqrt(H_SI / g_SI)  # reference time       [s]
p_SI = rho_w_SI * U_SI**2  # pressure scale = rho_w g H [Pa]

At = (rho_w_SI - rho_g_SI) / (rho_w_SI + rho_g_SI)  # Atwood number (~0.998)

# Reference scales in code units (by construction = 1)
rho_H = 1.0  # heavy (water) density
H = 1.0  # length
g0 = 1.0  # background accel magnitude
U_c = math.sqrt(g0 * H)

# Derived fluid properties in code units
rho_L = rho_g_SI / rho_w_SI  # light (gas) density       (~1.2e-3)
# sigma from We = rho_w c_gw^2 lambda_min / sigma with c_gw^2 = g0 lambda_min / (2 pi)
lambda_min = lambda_min_over_H * H  # shortest seeded wavelength (code units)
c_gw = math.sqrt(g0 * lambda_min / (2.0 * math.pi))  # gravity-wave speed at lambda_min
sigma = rho_H * c_gw**2 * lambda_min / We  # = g0 lambda_min^2 / (2 pi We)
We_H = rho_H * g0 * H**2 / sigma  # legacy H-based Weber (= 512 pi We), for reporting
# Per-fluid Reynolds numbers Re(i) = rho_w U_c H / mu_i  (global reference density rho_w)
Re_w = rho_w_SI * U_SI * H_SI / mu_w_SI  # water Reynolds (~1.4e4) -> fluid_pp(1)%Re(1)
Re_g = rho_w_SI * U_SI * H_SI / mu_g_SI  # gas   Reynolds          -> fluid_pp(2)%Re(1)

# 3. Forcing:  accel_y(t) = g_y + k_y sin(w_y t - p_y)   (dimensionless; fluid-agnostic)
g_sign = -1.0  # heavy on the bottom; -1 => statically STABLE background
forcing_ratio = args.forcing_ratio  # k_y / g0  (oscillation amplitude / background); 0 => pure RT
forcing_freq = args.forcing_freq  # w_y * t_c  (dimensionless angular frequency)
forcing_phase = 0.0

g_y = g_sign * g0
k_y = forcing_ratio * g0
w_y = forcing_freq / math.sqrt(H / g0)  # = forcing_freq in code units

# Physical (SI) forcing values: the acceleration scale is g_SI (code g0=1 -> g_SI)
# and the time scale is t_SI, so an angular frequency in code units converts as
# w_SI = w_code / t_SI.
g_y_SI = g_y * g_SI  # background accel magnitude/sign [m/s^2]
k_y_SI = k_y * g_SI  # oscillation accel amplitude     [m/s^2]
w_y_SI = w_y / t_SI  # angular forcing frequency       [rad/s]
f_y_SI = w_y_SI / (2.0 * math.pi)  # forcing frequency               [Hz]
p_y = forcing_phase

# Run length expressed in oscillation periods of the forcing
n_periods = 5.0  # number of forcing periods to simulate
frames_per_period = 40.0  # output snapshots saved per period
T_osc = 2.0 * math.pi / w_y
t_stop = n_periods * T_osc
t_save = T_osc / frames_per_period

# 4. Equation of state + reference pressure (with stiffness reduction for water)
# p_ref must exceed the body-force-driven pressure swing across the column so that
# pressure stays positive everywhere.  Swing ~ rho_H * (|g_y| + |k_y|) * Ly.
Ly_over_H = 2.0  # set again in section 5; used here for p_ref
p_swing = rho_H * (abs(g_y) + abs(k_y)) * Ly_over_H
pressure_safety = 2.5
p_ref = pressure_safety * p_swing  # ~110 code units

# Gas (ideal) sound speed at p_ref:  c^2 = gam_g p_ref / rho_g
c_gas = math.sqrt(gam_g * p_ref / rho_L)

# Water pi_inf in code pressure units (the true physical value)
pi_inf_w_phys = pi_inf_w_SI / p_SI  # ~2.3e6 code units

# Stiffness-reduction choice for water:
#   "match_gas" -> lower pi_inf so c_water == c_gas (cost-optimal; water still Mach<<1)
#   "physical"  -> use the true pi_inf (correct, but ~9x smaller dt)
# The softening only buys a larger acoustic step, which the projection does not need
water_pi_inf_mode = "match_gas" if args.explicit else "physical"
water_pi_inf_mode = "physical"

if water_pi_inf_mode == "match_gas":
    # c_water^2 = gam_w (p_ref + pi_inf) / rho_w  ==  c_gas^2
    pi_inf_w = c_gas**2 * rho_H / gam_w - p_ref
elif water_pi_inf_mode == "physical":
    pi_inf_w = pi_inf_w_phys
else:
    raise ValueError(f"unknown water_pi_inf_mode {water_pi_inf_mode!r}")

c_water = math.sqrt(gam_w * (p_ref + pi_inf_w) / rho_H)
softening_factor = pi_inf_w_phys / pi_inf_w  # how much we softened vs reality

# MFC stiffened-gas coefficients:  Gamma = 1/(gam-1),  Pi = gam pi_inf / (gam-1)
Gamma_w = 1.0 / (gam_w - 1.0)
Pi_w = gam_w * pi_inf_w / (gam_w - 1.0)
Gamma_g = 1.0 / (gam_g - 1.0)
Pi_g = 0.0  # ideal gas

c_max = max(c_water, c_gas)
Mach_w = U_c / c_water
Mach_g = U_c / c_gas

# 5. Domain & grid
Lx = 2.0 * H  # periodic width  -> mode wavelengths lambda = Lx / n
Ly = 2.0 * H  # tall enough for the interface to deform and break
x0, y0 = 0.0, 0.0
Ny = 1024
Nx = int(Ny * Lx / Ly)  # aspect-ratio-preserving resolution
dx, dy = Lx / Nx, Ly / Ny

# Acoustic-limited timestep estimate (adaptive CFL; the gas/water sound speed dominates)
dt_est = 0.4 * min(dx, dy) / c_max
n_steps_est = t_stop / dt_est

# 6. Perturbation seed (the statistical quantities)
# Seeded wavelength band in PHYSICAL units (fractions of H), anchored to the
# physics rather than the box. Because modes are selected by physical k, the
# shortest seeded wave -- and hence cells-per-shortest-wave -- is independent of
# Lx, Ly.
lambda_max = 1.0 * H  # longest seeded wavelength
# lambda_min (= lambda_min_over_H * H) is defined in section 1: the We definition
# is anchored to it, so it participates in the length-scale inversion there.
# It sets cells/shortest-wave here.
p_slope = 0.0  # spectral slope a_n^2 ~ k^p (0 = broadband)
seed = 42
eta_rms = 0.03 * lambda_min  # RMS interface amplitude

# The interface profile is linearly interpolated onto the DNS grid by hcid 209,
# so its resolution is DECOUPLED from (Nx, Ny). It only needs enough samples
# to resolve the shortest seeded wave (lambda_min); ~16 points/wave is ample.
# Across the box that is pts_per_lambda_pert * (Lx / lambda_min) samples.
pts_per_lambda_pert = 12
Nx_pert = max(128, math.ceil(pts_per_lambda_pert * Lx / lambda_min))

x, eta, diag = generate_perturbation_2d(
    Lx,
    Nx_pert,
    eta_rms,
    lambda_min=lambda_min,
    lambda_max=lambda_max,
    p=p_slope,
    seed=seed,
    randomize_amp=True,
    x0=x0,
)
write_interface_2d(x, eta, y_mean=(y0 + H), path=os.path.join(case_dir, "interface_profile.dat"))

# Interface smoothing half-thickness (a few cells) handed to hcid 209 via a(4)
delta = 4.0 * dy

gam = gam_w  # (kept for the Mach reporting helpers below; per-fluid EOS set in `data`)

# 7. Diagnostics (stderr; stdout must be pure JSON for MFC)
print_diagnostics(diag, file=sys.stderr)
print(f"  cells/shortest-wave  {cells_per_lambda(diag['lambda_min'], dy):.1f}", file=sys.stderr)
print(f"  total cells {Nx*Ny:.3g}  (Nx={Nx} Ny={Ny})", file=sys.stderr)
print(f"  physical scales: H={H_SI*100:.3g} cm  U_c={U_SI:.3g} m/s  t_c={t_SI*1e3:.3g} ms  p_c={p_SI:.4g} Pa", file=sys.stderr)
print(f"  At={At:.5f}  We={We:.4g} (at lambda_min; legacy We_H={We_H:.4g})  " f"Re_water={Re_w:.4g}  Re_gas={Re_g:.4g}", file=sys.stderr)
print(f"  c_gw(lambda_min)={c_gw:.4g}  lambda_min={lambda_min:.4g} ({lambda_min_SI*1e3:.3g} mm)", file=sys.stderr)
print(f"  rho_H(water)={rho_H:.4g} rho_L(gas)={rho_L:.4g}  sigma={sigma:.4g}", file=sys.stderr)
print(f"  EOS mode='{water_pi_inf_mode}'  pi_inf_water(code)={pi_inf_w:.4g} (physical={pi_inf_w_phys:.4g}, softened {softening_factor:.3g}x)", file=sys.stderr)
print(f"  p_ref={p_ref:.4g} (swing~{p_swing:.3g}, min p~{p_ref - p_swing:.3g})", file=sys.stderr)
print(f"  c_water={c_water:.4g} (Mach {Mach_w:.3g})  c_gas={c_gas:.4g} (Mach {Mach_g:.3g})", file=sys.stderr)
print(f"  dt_est~{dt_est:.3g}  n_steps_est~{n_steps_est:.3g}  (t_stop={t_stop:.4g})", file=sys.stderr)
print(f"  lambda_dom/H={diag['lambda_dom'] / H:.3g}", file=sys.stderr)
print(f"  T_osc={T_osc:.4g}  n_periods={n_periods:.0f}  t_save={t_save:.4g} ({frames_per_period:.0f}/period)", file=sys.stderr)
print(f"  physical forcing: g_y={g_y_SI:.4g} m/s^2  k_y={k_y_SI:.4g} m/s^2  w_y={w_y_SI:.4g} rad/s  f_y={f_y_SI:.4g} Hz  T_osc={T_osc * t_SI * 1e3:.4g} ms", file=sys.stderr)

# Persist diagnostics + case parameters to a file for the run record
save_diagnostics(
    diag,
    path=os.path.join(case_dir, "diagnostics.dat"),
    extra={
        "At": At,
        "We(lambda_min)": We,
        "We_H(legacy)": We_H,
        "c_gw(lambda_min)": c_gw,
        "lambda_min": lambda_min,
        "Re_water": Re_w,
        "Re_gas": Re_g,
        "H_SI[m]": H_SI,
        "lambda_min_SI[m]": lambda_min_SI,
        "U_SI[m/s]": U_SI,
        "t_SI[s]": t_SI,
        "p_SI[Pa]": p_SI,
        "rho_H(water)": rho_H,
        "rho_L(gas)": rho_L,
        "sigma": sigma,
        "water_pi_inf_mode": water_pi_inf_mode,
        "pi_inf_water(code)": pi_inf_w,
        "pi_inf_water_physical(code)": pi_inf_w_phys,
        "softening_factor": softening_factor,
        "p_ref": p_ref,
        "c_water": c_water,
        "c_gas": c_gas,
        "Mach_water": Mach_w,
        "Mach_gas": Mach_g,
        "dt_est": dt_est,
        "n_steps_est": n_steps_est,
        "lambda_dom/H": diag["lambda_dom"] / H,
        "cells/shortest-wave": cells_per_lambda(diag["lambda_min"], dy),
        "total cells": Nx * Ny,
        "Nx": Nx,
        "Ny": Ny,
        "g_y": g_y,
        "k_y": k_y,
        "w_y": w_y,
        "p_y": p_y,
        "g_y[m/s^2]": g_y_SI,
        "k_y[m/s^2]": k_y_SI,
        "w_y[rad/s]": w_y_SI,
        "f_y[Hz]": f_y_SI,
        "T_osc[s]": T_osc * t_SI,
        "T_osc": T_osc,
        "n_periods": n_periods,
        "t_stop": t_stop,
        "t_save": t_save,
        "frames_per_period": frames_per_period,
    },
)

eps = 1.0e-8

data = {
    # Logistics
    "run_time_info": "T",
    "rdma_mpi": "F",
    # Computational domain (2D: p = 0, no z direction)
    "x_domain%beg": x0,
    "x_domain%end": x0 + Lx,
    "y_domain%beg": y0,
    "y_domain%end": y0 + Ly,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": 0,
    "cyl_coord": "F",
    # MTHINC
    "int_comp": 2,
    "ic_beta": 1.5,
    # CDI
    # "int_comp": 4,
    # "ic_delta": 2.0,
    # Time stepping (adaptive CFL; code time units, t_c = 1)
    "cfl_adap_dt": "T",
    "cfl_target": 0.4,
    "t_stop": t_stop,
    "t_save": t_save,
    "n_start": 0,
    # Simulation algorithm
    "model_eqns": 2,
    "num_fluids": 2,
    "num_patches": 1,
    "alt_soundspeed": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 3,
    "weno_eps": 1e-8,
    "weno_avg": "F",
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -1,  # periodic
    "bc_x%end": -1,
    # Database output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    # "cons_vars_wrt": "T",
    # "cf_wrt": "T",
    "parallel_io": "T",
    "file_per_process": "T",
    # Equation of state:
    #   fluid 1 = WATER  -> stiffened gas (gamma_w, pi_inf_w)
    #   fluid 2 = GAS    -> ideal gas     (gamma_g, pi_inf = 0)
    "fluid_pp(1)%gamma": Gamma_w,
    "fluid_pp(1)%pi_inf": Pi_w,
    "fluid_pp(2)%gamma": Gamma_g,
    "fluid_pp(2)%pi_inf": Pi_g,
    # Single full-domain patch; hcid 209 builds the perturbed two-fluid interface.
    # The constant state below is a placeholder (the hcid overwrites alpha/alpha_rho/
    # color everywhere); only pres and vel are kept. Interface parameters are passed
    # through %a():  a(2)=rho_H (water, below), a(3)=rho_L (gas, above), a(4)=delta.
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%hcid": 209,
    "patch_icpp(1)%x_centroid": x0 + 0.5 * Lx,
    "patch_icpp(1)%y_centroid": y0 + 0.5 * Ly,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p_ref,
    "patch_icpp(1)%alpha_rho(1)": (1.0 - eps) * rho_H,
    "patch_icpp(1)%alpha_rho(2)": eps * rho_L,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": 1.0,
    "patch_icpp(1)%a(2)": rho_H,
    "patch_icpp(1)%a(3)": rho_L,
    "patch_icpp(1)%a(4)": delta,
    # Body force: background + oscillatory acceleration in y
    "bf_y": "T",
    "g_y": g_y,
    "k_y": k_y,
    "w_y": w_y,
    "p_y": p_y,
}

if args.explicit:
    # Acoustic-CFL settings of the original case; the projection rejects mpp_lim and has no Riemann solver to correct
    data.update({"mpp_lim": "T", "low_Mach": 2})
else:
    data.update({"proj_method": "T"})

# if (We < 100):
data.update(
    {
        # Surface tension (Weber number)
        "surface_tension": "T",
        "sigma": sigma,
        # The well-balanced model needs the projection's face pressure gradient
        "surface_tension_model": "conservative" if args.explicit else "well_balanced",
    }
)

# if (Re_w < 1e5) or (Re_g < 1e5):
data.update(
    {
        # Viscosity (per-fluid Reynolds numbers: water and gas differ by ~mu_w/mu_g ~ 56x)
        "viscous": "T",
        "weno_Re_flux": "F",
        "fluid_pp(1)%Re(1)": Re_w,
        "fluid_pp(2)%Re(1)": Re_g,
        "bc_y%beg": -16,  # No slip wall
        "bc_y%end": -16,
    }
)
# else:
# data.update({
# "bc_y%beg": -15,   # Slip wall
# "bc_y%end": -15,
# })

print(json.dumps(data, indent=4))
