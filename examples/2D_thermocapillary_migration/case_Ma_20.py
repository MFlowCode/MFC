#!/usr/bin/env python3
# Thermocapillary migration of a 2D drop at finite Marangoni number -- Samareh, Mostaghimi & Moreau,
# Int. J. Heat Mass Transfer 73 (2014) 616-626, Sec. 4.1.2 / Fig. 7 (test case of Nas & Tryggvason).
# A drop of diameter D = 1 starts 1D above the cold floor of a 2D x 4D slip-wall box. With Ma finite
# the energy equation is coupled: the drop's motion distorts the temperature field, a thermal
# boundary layer forms at the interface, and sigma responds to the evolving local T. Drop and bulk
# are genuinely DIFFERENT fluids (all material properties at ratio 0.5). Non-dimensional targets:
# Re = 5, Ma = 20, Ca = 0.01666; Samareh's migration velocity U* = U/U_r peaks near 0.13 at t*~5.
#
# TEMPERATURE WITH DISTINCT FLUIDS (compressible code, incompressible reference).
# MFC is compressible: at the low migration Mach here, pressure is ~uniform and the EOS locks
# T = (p+p_inf)/((gam-1)*rho*cv), so a temperature gradient is a density gradient. With two distinct
# fluids we impose T(y) by a PER-FLUID density proxy: each pure-fluid density absorbs the field,
#   rho_i(x,y) = (p + p_inf_i)/((gam-1)*cv_i*T(y)).
# The fluid-to-fluid density RATIO is then height-independent (T(y) cancels), so the prescribed
# Nas-Tryggvason ratio 0.5 is preserved by tuning the per-fluid stiffening p_inf_i; and the MIXTURE
# EOS recovers exactly T_mix(y) = T(y) across the interface (the volume fractions cancel) -- no jump.
# Bulk conduction + isothermal gradient walls then sustain the field against the drop's advection.
# Honest caveat: each fluid's ABSOLUTE density stratifies as ~1/T (a compressibility artifact absent
# in Samareh's incompressible solver; magnitude ~ dT/T, shrunk by the T_base offset). So the peak
# matches Fig 7 only qualitatively; the overshoot-then-settle shape is the validation target.

import json

# Non-dimensional targets (Samareh Fig. 7)
Re = 5.0
Ma = 20.0
Ca = 0.01666
prop_ratio = 0.5  # droplet / bulk material-property ratio (rho, mu, cv, k all at 0.5)

# Geometry: D = 1 drop in a 2D wide x 4D tall box; drop 1D above the cold floor, gradient axis = y
D = 1.0
r = D / 2.0
Wx = 2.0 * D
Hy = 4.0 * D
y_bottom = -Hy / 2.0
y_drop = y_bottom + 1.0 * D

Nx = 64  # cells across the box width (Samareh used 64, 128)
dx = Wx / Nx
Ny = round(Hy / dx)

# Bulk reference scales (arbitrary; only Re, Ma, Ca are physical)
rho_b = 1.0
mu_b = 0.02
gam = 2.0

# Invert the non-dimensional numbers for the physical surface-tension and conduction properties
G = Re * mu_b**2 / (rho_b * r**2)  # |sigma_T*gradT| = 0.008
gradT = 1.0 / Hy  # T runs 1 (cold) .. 2 (hot) across the box height = 0.25
sigma_T = -G / gradT  # dsigma/dT < 0 = -0.032
sigma0 = G * r / Ca  # surface tension at T_ref ~ 0.240
alpha_b = G * r**2 / (mu_b * Ma)  # bulk thermal diffusivity = 0.005

# Imposed field T(y) = T_base + gradT*(y - y_bottom). T_base > 0 keeps rho positive and the
# isothermal-wall validator happy; a larger offset reduces the ~1/T compressibility stratification.
T_base = 1.0


def T_of_y(y):
    return T_base + gradT * (y - y_bottom)


T_ref = T_of_y(y_drop)  # sigma = sigma0 at the drop's initial position; T_ref = 1.25

# Stiffened-gas EOS, bulk fluid. Pick a reference sound speed (low migration Mach ~ U_r/c) and
# background pressure (well above the Laplace jump sigma0/r); the stiffening and cv then close the
# EOS so the bulk density is rho_b at (p0, T_ref).
c_ref = 15.0
p0 = 5.0
p_inf_b = c_ref**2 * rho_b / gam - p0  # bulk stiffening
cv_b = (p0 + p_inf_b) / ((gam - 1.0) * rho_b * T_ref)
cp_b = gam * cv_b
k_b = alpha_b * rho_b * cp_b  # bulk conductivity

# Droplet = prop_ratio * bulk for every material property. The EOS stiffening p_inf_d follows so the
# drop density is rho_d at (p0, T_ref): p0+p_inf_d = prop_ratio^2*(p0+p_inf_b) -> ratio rho_d/rho_b = 0.5.
rho_d = prop_ratio * rho_b
mu_d = prop_ratio * mu_b
cv_d = prop_ratio * cv_b
k_d = prop_ratio * k_b
p_inf_d = (gam - 1.0) * rho_d * cv_d * T_ref - p0

U_r = G * r / mu_b  # Marangoni velocity scale = 0.2
t_r = mu_b / G  # capillary-thermal time = 2.5

# Time stepping: min(acoustic CFL, explicit-conduction limit). Density is lowest (sound speed
# highest) at the hot wall, where rho ~ rho_ref*T_ref/T_hot.
T_hot = T_of_y(Hy / 2.0)
rho_hot_min = prop_ratio * rho_b * T_ref / T_hot  # drop fluid at the hot wall: global density min
c_max = (gam * (p0 + p_inf_d) / rho_hot_min) ** 0.5
alpha_d = k_d / (rho_d * cp_b * prop_ratio)  # = k_d/(rho_d*cp_d), the fastest-diffusing phase
alpha_max = max(alpha_b, alpha_d)
mydt = 0.35 * dx / c_max
mydt = min(mydt, 0.35 * dx**2 / (4.0 * alpha_max))
t_step_stop = round(15.0 * t_r / mydt)  # 15 capillary-thermal times
t_step_save = max(1, t_step_stop // 100)

# One analytic patch over the whole box (a second patch's analytic density would leak globally). The
# smooth circle eta(x,y) ~ 1 in the drop / 0 outside drives the color, the volume-fraction split, the
# Laplace pressure jump, and BOTH per-fluid densities together. Hardcode the drop center/radius as
# decimals (`r` and `eps` are reserved IC-expression names).
xc_d, yc_d, r_d = 0.0, y_drop, r
w_if = 0.75 * dx  # interface half-width (~3-cell transition)
laplace = sigma0 / r  # Laplace pressure jump sigma/r
dist = f"sqrt((x - ({xc_d:.9f}))**2 + (y - ({yc_d:.9f}))**2)"
eta = f"0.5*(1.0 - tanh(({dist} - {r_d:.9f})/{w_if:.9f}))"
Texpr = f"({T_base:.9f} + {gradT:.9f}*(y - ({y_bottom:.9f})))"
pres_expr = f"({p0:.9f} + {laplace:.9f}*({eta}))"
# volume fractions: drop = eta, bulk = 1 - eta, both floored at a trace eps_vf so no phase is ever
# exactly empty (pressure relaxation rejects a cell holding volume fraction but no mass)
eps_vf = 1.0e-9
alpha2_expr = f"({eps_vf!r} + (1.0 - 2.0*{eps_vf!r})*({eta}))"
alpha1_expr = f"(1.0 - {alpha2_expr})"
# partial densities: alpha_i * rho_i, with rho_i = (p + p_inf_i)/((gam-1)*cv_i*T(y))
cb = (gam - 1.0) * cv_b
cd = (gam - 1.0) * cv_d
arho1_expr = f"{alpha1_expr}*({pres_expr} + {p_inf_b:.9f})/({cb:.9f}*{Texpr})"
arho2_expr = f"{alpha2_expr}*({pres_expr} + {p_inf_d:.9f})/({cd:.9f}*{Texpr})"
cf_expr = f"({eta})"

# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": -Wx / 2,
    "x_domain%end": Wx / 2,
    "y_domain%beg": -Hy / 2,
    "y_domain%end": Hy / 2,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": 0,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    # Simulation Algorithm
    "model_eqns": 3,
    "alt_soundspeed": "F",
    "mixture_err": "T",
    "mpp_lim": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "weno_avg": "T",
    "weno_Re_flux": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    # Closed free-slip box (Samareh / Nas & Tryggvason): reflective side walls, slip walls on y
    # (isothermal walls require -15/-16)
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -15,
    "bc_y%end": -15,
    "num_patches": 1,
    "num_fluids": 2,
    # Physics: viscosity + temperature-dependent surface tension sigma(T); k_therm > 0 enables conduction
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": sigma_T,
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    # Fluid 1 (bulk)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf_b / (gam - 1.0),
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu_b,
    "fluid_pp(1)%k_therm": k_b,
    # Fluid 2 (droplet), all material properties at prop_ratio * bulk; stiffening keeps rho_d/rho_b = 0.5
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf_d / (gam - 1.0),
    "fluid_pp(2)%cv": cv_d,
    "fluid_pp(2)%Re(1)": 1.0 / mu_d,
    "fluid_pp(2)%k_therm": k_d,
    # Single analytic patch: the drop lives entirely in eta(x,y). Color, volume fractions, pressure,
    # and both per-fluid densities all share it, so the mixture EOS recovers the linear T(y) exactly
    # (drop included) while the distinct fluid-to-fluid density ratio (0.5) is preserved.
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": Wx,
    "patch_icpp(1)%length_y": Hy,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": pres_expr,
    "patch_icpp(1)%alpha_rho(1)": arho1_expr,
    "patch_icpp(1)%alpha_rho(2)": arho2_expr,
    "patch_icpp(1)%alpha(1)": alpha1_expr,
    "patch_icpp(1)%alpha(2)": alpha2_expr,
    "patch_icpp(1)%cf_val": cf_expr,
    # Isothermal Dirichlet gradient walls pin T to the imposed gradient (cold floor / hot ceiling)
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T_of_y(-Hy / 2),
    "bc_y%Twall_out": T_of_y(Hy / 2),
}

print(json.dumps(data))
