#!/usr/bin/env python3
# Phase 1 drainage validation case — acoustic_demulsification project
#   Geometry:  2D axisymmetric, head-on collision of two equal droplets
#   Fluids  :  tetradecane drops in nitrogen ambient
#   Target  :  resolve gas-film drainage; compare h_min(t) to Klaseboer 2000
#   Stop    :  when h_min reaches ~100–200 nm (before numerical rupture)
# Roadmap reference: +Projects/acoustic_demulsification/simulation_roadmap.md
# (Phase 1 — Drainage validation)
# Run:
#   ./mfc.sh run case.py -n <N>          # single node
#   ./mfc.sh run case.py --gpu -n <N>    # GPU build
# Knobs at top of file:
#   SCHEME         — 'wenoz' or 'muscl_thinc'
#   DX_FINE        — cell size in drainage film (m); start 10 nm
#   DX_BULK        — cell size in bulk (m); start 1 micron
#   FILM_HALF      — axial half-width of the uniform fine zone (m)
#   WE             — target Weber number (default 30)

import json
import math

# User-facing knobs

SCHEME = "wenoz"  # 'wenoz' | 'muscl_thinc'
DX_FINE = 10e-9  # 10 nm — film cell size target
DX_BULK = 1e-6  #  1 um — bulk cell size target
FILM_HALF = 2e-6  # ± half-width of the uniform fine zone on the axis (m)
WE = 30.0  # Weber number (coalescence regime)

# For quick debugging, set coarser values, e.g. DX_FINE=40e-9, DX_BULK=4e-6.

# Physical constants (tetradecane / nitrogen at ~25 C)

# Tetradecane (liquid)
rho_l = 762.0  # kg/m^3
mu_l = 2.1e-3  # Pa*s
sigma = 0.027  # N/m (≈ 27 mN/m)

# Stiffened-gas fit for tetradecane: gamma = 4.4, p_inf = 3.0e8 Pa
# => c_l = sqrt(gamma*(p+p_inf)/rho) ≈ 1316 m/s at 1 atm — matches bulk data.
gamma_l = 4.4
pi_inf_l = 3.0e8

# Nitrogen (ideal gas)
rho_g = 1.165  # kg/m^3 at 1 atm, 20 C
gamma_g = 1.4
pi_inf_g = 0.0

p_amb = 1.01325e5  # 1 atm

# Droplet geometry and kinematics

D = 300e-6  # droplet diameter, 300 um
R = D / 2.0  # radius, 150 um

# Weber number: We = rho_l * U_rel^2 * D / sigma,  U_rel = 2U per Qian & Law.
# U_rel^2 = We*sigma/(rho_l*D)
#         = 30 * 0.027 / (762 * 3e-4)
#         = 0.810   / 0.2286
#         = 3.5433 m^2/s^2
# U_rel   = 1.8823 m/s  ->  U = 0.9412 m/s per droplet
U_rel = math.sqrt(WE * sigma / (rho_l * D))
U = 0.5 * U_rel

# Initial axial gap between droplet surfaces (pre-drainage), and centers.
gap0 = 20e-6  # 20 um initial gap — leaves inertial approach + drainage
x_center = R + 0.5 * gap0  # droplet centers at x = +/- x_center

# Domain & grid

# Axial extent — big enough that outgoing waves don't reflect back onto the film.
# Keep a few D each side of the droplets.
Lx_half = 5.0 * R  # = 750 um axial half-domain
Ly_max = 3.0 * R  # = 450 um radial extent

# Uniform fine zone in x: [-FILM_HALF, +FILM_HALF] at DX_FINE.
# Stretched zone: DX_FINE -> DX_BULK, geometric growth handled by MFC stretch_x.
# Cell-count estimate (for sanity; MFC's hyperbolic stretch won't be exactly geometric):
#     N_fine = 2*FILM_HALF / DX_FINE
#     N_stretch (each side) ~ ln(DX_BULK/DX_FINE)/ln(r)   with r = 1.1 -> ~48
N_fine = int(2 * FILM_HALF / DX_FINE)  #  400 for defaults
N_stretch = int(math.ceil(math.log(DX_BULK / DX_FINE) / math.log(1.1)))  # ~48 each side
Nx = N_fine + 2 * N_stretch  # 496 for defaults

# Radial: fine near axis (film rim), coarsens outward.
# Same strategy on y (positive only, because axisymmetric).
N_fine_r = int(4 * FILM_HALF / DX_FINE)  # 800 cells out to 8 um
N_stretch_r = int(math.ceil(math.log(DX_BULK / DX_FINE) / math.log(1.1)))  # ~48
Ny = N_fine_r + N_stretch_r  # 848 for defaults

# NOTE: These are generous baselines for a grid-convergence study.
#       On first debug runs, set DX_FINE=40e-9, DX_BULK=4e-6 to cut Nx,Ny by ~4x.

# Time integration

# Characteristic times:
#   tau_inertia  = D / U_rel  = 3e-4 / 1.8823 ≈ 1.59e-4 s
#   tau_cap      = sqrt(rho_l D^3 / sigma)     ≈ 8.73e-4 s
# Want to simulate until h_min ~ 100–200 nm. Experience says ~1–2 * tau_inertia.
t_stop = 2.0e-4  # 200 us
t_save = 5.0e-6  # 40 frames

# Numerical scheme switches

# WENO-Z 5 is the default from MFC 5.0 paper; MUSCL+THINC is the fallback if
# interface width exceeds ~5 cells or if we see early numerical rupture.
if SCHEME == "wenoz":
    recon = {
        "recon_type": 1,
        "weno_order": 5,
        "weno_eps": 1e-6,
        "mapped_weno": "F",
        "wenoz": "T",
        "mp_weno": "T",
        "teno": "F",
    }
elif SCHEME == "muscl_thinc":
    recon = {
        "recon_type": 2,
        "muscl_order": 2,
        "muscl_lim": 4,  # Van Leer
        "int_comp": "T",  # THINC interface compression
        "ic_eps": 1e-4,
        "ic_beta": 1.6,
    }
else:
    raise ValueError(f"Unknown SCHEME: {SCHEME!r}")

# MFC input deck

eps = 1e-9  # avoid exactly 0/1 volume fractions

case = {
    # grid
    "x_domain%beg": -Lx_half,
    "x_domain%end": +Lx_half,
    "y_domain%beg": 0.0,
    "y_domain%end": +Ly_max,
    "m": Nx - 1,  # MFC uses zero-based cell-count (m = Nx-1)
    "n": Ny - 1,
    "p": 0,  # 2D
    "cyl_coord": "T",  # 2D axisymmetric: x is axial, y is radial
    # Grid stretching — fine uniform zone in [-FILM_HALF, +FILM_HALF], stretched outside.
    "stretch_x": "T",
    "a_x": 4.0,
    "x_a": -FILM_HALF,
    "x_b": +FILM_HALF,
    "loops_x": 2,
    "stretch_y": "T",
    "a_y": 4.0,
    "y_a": 0.0,
    "y_b": 4.0 * FILM_HALF,  # keep fine near the axis / film rim
    "loops_y": 2,
    # boundary conditions
    # Axial (x): outflow/extrapolation on both ends.
    "bc_x%beg": -3,
    "bc_x%end": -3,
    # Radial (y): axisymmetric at y=0, outflow at y_max.
    "bc_y%beg": -2,  # reflective / axis
    "bc_y%end": -3,
    # time integration
    "time_stepper": 3,  # TVD-RK3
    "cfl_adap_dt": "T",
    "cfl_target": 0.3,  # conservative; tighten/loosen after first runs
    "t_stop": t_stop,
    "t_save": t_save,
    # physical model
    "model_eqns": 2,  # 5-equation (Allaire) diffuse interface
    "num_fluids": 2,
    "num_patches": 3,
    "mixture_err": "T",
    "mpp_lim": "T",
    "avg_state": 2,
    "riemann_solver": 2,  # HLLC
    "wave_speeds": 1,
    **recon,
    # surface tension
    "surface_tension": "T",
    "sigma": sigma,
    # fluids (index 1 = tetradecane, 2 = nitrogen)
    "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),  # 1/(4.4-1) = 0.2941
    "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),  # 4.4*3e8/3.4 ≈ 3.882e8
    "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),  # 1/(1.4-1) = 2.5
    "fluid_pp(2)%pi_inf": 0.0,
    # patch 1: nitrogen background filling the whole domain
    "patch_icpp(1)%geometry": 3,  # rectangle (2D)
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": Ly_max / 2.0,
    "patch_icpp(1)%length_x": 2.5 * Lx_half,  # cover w/ margin
    "patch_icpp(1)%length_y": 2.0 * Ly_max,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p_amb,
    "patch_icpp(1)%alpha_rho(1)": eps * rho_l,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_g,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    "patch_icpp(1)%cf_val": 0.0,
    # patch 2: left droplet, moving in +x
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.5,
    "patch_icpp(2)%geometry": 2,  # circle (axisymmetric disk becomes a sphere)
    "patch_icpp(2)%x_centroid": -x_center,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radius": R,
    "patch_icpp(2)%vel(1)": +U,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p_amb,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(2)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "patch_icpp(2)%cf_val": 1.0,
    # patch 3: right droplet, moving in -x
    "patch_icpp(3)%alter_patch(1)": "T",
    "patch_icpp(3)%smoothen": "T",
    "patch_icpp(3)%smooth_patch_id": 1,
    "patch_icpp(3)%smooth_coeff": 0.5,
    "patch_icpp(3)%geometry": 2,
    "patch_icpp(3)%x_centroid": +x_center,
    "patch_icpp(3)%y_centroid": 0.0,
    "patch_icpp(3)%radius": R,
    "patch_icpp(3)%vel(1)": -U,
    "patch_icpp(3)%vel(2)": 0.0,
    "patch_icpp(3)%pres": p_amb,
    "patch_icpp(3)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(3)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(3)%alpha(1)": 1.0 - eps,
    "patch_icpp(3)%alpha(2)": eps,
    "patch_icpp(3)%cf_val": 1.0,
}

# Derived-quantity banner (printed to stderr on --case-optimization runs)
# Printed here as a comment for humans; MFC only reads the JSON below.
#   U_rel  = 1.8823 m/s          (We = 30, head-on: U_rel = 2U)
#   U      = 0.9412 m/s per drop
#   h_c    = (A R / sigma)^(1/3) = 38 nm  (A = 1e-20 J)
#   tau_i  = D/U_rel             ≈ 1.59e-4 s
#   tau_s  = sqrt(rho D^3/sigma) ≈ 8.73e-4 s
#   Re     = rho U_rel D / mu    ≈ 205
#   Nx     = N_fine + 2*N_stretch  (defaults: 400 + 96 = 496)
#   Ny     = N_fine_r + N_stretch_r (defaults: 800 + 48 = 848)

print(json.dumps(case))
