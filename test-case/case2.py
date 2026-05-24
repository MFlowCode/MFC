#!/usr/bin/env python3
"""
2D single-material JWL shock driver interacting with a random particle bed
-------------------------------------------------------------------------

This is a numerical JWL EOS benchmark, not a full explosive-products/air
multi-material detonation case. The current JWL implementation supports one
fluid, so the background is a low-pressure JWL material and the left driver is
a high-pressure JWL region. The pressure jump launches a shock-like compression
wave into the particle bed.
"""

import json
import math

# =============================================================================
# 1. CASE NAME
# =============================================================================

case_name = "jwl_shock_driver_2D_particle_bed"

# =============================================================================
# 2. PARTICLE / BED SETTINGS
# =============================================================================

dp0 = 5.0e-3
D = dp0

Grid_Res = 8

target_num_particles = 10
random_seed = 5000

particle_min_spacing = 0.10 * D

# 0 = fixed
# 2 = fully coupled moving IBM
moving_ibm_mode = 2

# =============================================================================
# 3. DOMAIN
# =============================================================================

bed_x_start_in_D = 2.0
bed_x_len_in_D   = 8.0
bed_y_len_in_D   = 8.0

bed_x0 = bed_x_start_in_D * D
bed_x1 = (bed_x_start_in_D + bed_x_len_in_D) * D

bed_y0 = -0.5 * bed_y_len_in_D * D
bed_y1 = +0.5 * bed_y_len_in_D * D

bed_lx = bed_x1 - bed_x0
bed_ly = bed_y1 - bed_y0

bed_x_centroid = 0.5 * (bed_x0 + bed_x1)
bed_y_centroid = 0.5 * (bed_y0 + bed_y1)

# Padding
x_pad_left_D  = 8.0
x_pad_right_D = 8.0
y_pad_D       = 2.0

L_xbeg = bed_x0 - x_pad_left_D * D
L_xend = bed_x1 + x_pad_right_D * D

L_ybeg = bed_y0 - y_pad_D * D
L_yend = bed_y1 + y_pad_D * D

x_len = L_xend - L_xbeg
y_len = L_yend - L_ybeg

# =============================================================================
# 4. GRID
# =============================================================================

dx = D / Grid_Res
dy = dx

m = max(1, int(round(x_len / dx)))
n = max(1, int(round(y_len / dy)))

p = 0

# =============================================================================
# 5. BACKGROUND JWL MATERIAL
# =============================================================================

# The JWL implementation in this branch is single-fluid. These values are a
# low-pressure JWL background used for a shock benchmark, not physical air.
gam_a = 1.4

p_ambient = 101325.0
rho_ambient = 1.22

# =============================================================================
# 6. JWL PARAMETERS
# =============================================================================

# TNT-like values

jwl_A     = 3.712e11
jwl_B     = 3.231e9
jwl_R1    = 4.15
jwl_R2    = 0.95
jwl_omega = 0.30
jwl_rho0  = 1630.0

# =============================================================================
# 7. SHOCK DRIVER INITIALIZATION
# =============================================================================

# A left-side high-pressure driver creates a clear right-going shock-like wave.
# This is much easier to benchmark than a tiny circular pressure spot.
driver_x_end = 1.5 * D
driver_x_centroid = 0.5 * (L_xbeg + driver_x_end)
driver_length_x = driver_x_end - L_xbeg

driver_pressure = 5.0e9
driver_velocity = 300.0
driver_density = jwl_rho0

# =============================================================================
# 8. PARTICLE MASS
# =============================================================================

r_ib = 0.5 * dp0

rho_particle = 0.5 * rho_ambient

particle_area = math.pi * r_ib**2

mass_particle = 1e-2   #rho_particle * particle_area

# =============================================================================
# 9. COLLISION MODEL
# =============================================================================

collision_model = 1

coefficient_of_restitution = 0.8
ib_coefficient_of_friction = 0.02

collision_time = 2.0e-4

# =============================================================================
# 10. TIME CONTROL
# =============================================================================

# Conservative JWL CFL
cfl = 0.50

# Estimated shock speed for output timing
shock_speed_est = 4000.0

# Characteristic particle-crossing time
tau = D / shock_speed_est

# Save frequently
t_save = 0.25 * tau

# Allow blast to fully cross bed
t_stop = 80.0 * tau

# =============================================================================
# 11. CASE DICTIONARY
# =============================================================================

case = {

    # -------------------------------------------------------------------------
    # OUTPUT
    # -------------------------------------------------------------------------

    "run_time_info": "T",

    "format": 1,
    "precision": 2,
    "parallel_io": "T",

    "prim_vars_wrt": "T",
    "cons_vars_wrt": "F",

    "rho_wrt": "T",
    "pres_wrt": "T",

    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",

    "E_wrt": "T",

    "c_wrt": "T",

    "ib_state_wrt": "T",

    # -------------------------------------------------------------------------
    # DOMAIN
    # -------------------------------------------------------------------------

    "x_domain%beg": float(L_xbeg),
    "x_domain%end": float(L_xend),

    "y_domain%beg": float(L_ybeg),
    "y_domain%end": float(L_yend),

    "cyl_coord": "F",

    "m": int(m),
    "n": int(n),
    "p": 0,

    # -------------------------------------------------------------------------
    # TIME
    # -------------------------------------------------------------------------

    "t_step_start": 0,

    "n_start": 0,

    "cfl_adap_dt": "T",
    "cfl_const_dt": "F",

    "cfl_target": float(cfl),

    "t_stop": float(t_stop),
    "t_save": float(t_save),

    # -------------------------------------------------------------------------
    # SOLVER
    # -------------------------------------------------------------------------

    "num_patches": 2,

    "model_eqns": 2,

    "alt_soundspeed": "F",

    "num_fluids": 1,

    "mpp_lim": "F",

    "mixture_err": "T",

    "time_stepper": 3,

    # safer initial JWL setup
    "weno_order": 3,

    "weno_eps": 1.0e-16,

    "weno_Re_flux": "T",

    "weno_avg": "F",

    "avg_state": 2,

    "mapped_weno": "F",

    "null_weights": "F",

    "mp_weno": "F",

    "riemann_solver": 2,

    "wave_speeds": 1,

    "viscous": "T",

    # -------------------------------------------------------------------------
    # BOUNDARY CONDITIONS
    # -------------------------------------------------------------------------

    "bc_x%beg": -3,
    "bc_x%end": -3,

    "bc_y%beg": -3,
    "bc_y%end": -3,

    # -------------------------------------------------------------------------
    # IBM SETTINGS
    # -------------------------------------------------------------------------

    "ib": "T",

    "num_ibs": 0,

    "collision_model": int(collision_model),

    "coefficient_of_restitution":
        float(coefficient_of_restitution),

    "ib_coefficient_of_friction":
        float(ib_coefficient_of_friction),

    "collision_time":
        float(collision_time),

    # -------------------------------------------------------------------------
    # INTERNAL PARTICLE BED
    # -------------------------------------------------------------------------

    "num_particle_beds": 1,

    "particle_bed(1)%x_centroid":
        float(bed_x_centroid),

    "particle_bed(1)%y_centroid":
        float(bed_y_centroid),

    "particle_bed(1)%z_centroid":
        0.0,

    "particle_bed(1)%length_x":
        float(bed_lx),

    "particle_bed(1)%length_y":
        float(bed_ly),

    "particle_bed(1)%length_z":
        0.0,

    "particle_bed(1)%num_particles":
        int(target_num_particles),

    "particle_bed(1)%radius":
        float(r_ib),

    "particle_bed(1)%mass":
        float(mass_particle),

    "particle_bed(1)%min_spacing":
        float(particle_min_spacing),

    "particle_bed(1)%moving_ibm":
        int(moving_ibm_mode),

    "particle_bed(1)%seed":
        int(random_seed),

    # -------------------------------------------------------------------------
    # PATCH 1 : LOW-PRESSURE JWL BACKGROUND
    # -------------------------------------------------------------------------

    "patch_icpp(1)%geometry": 3,

    "patch_icpp(1)%x_centroid":
        float(0.5 * (L_xbeg + L_xend)),

    "patch_icpp(1)%y_centroid":
        0.0,

    "patch_icpp(1)%length_x":
        float(x_len),

    "patch_icpp(1)%length_y":
        float(y_len),

    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,

    "patch_icpp(1)%pres":
        float(p_ambient),

    "patch_icpp(1)%alpha_rho(1)":
        float(rho_ambient),

    "patch_icpp(1)%alpha(1)":
        1.0,

    # -------------------------------------------------------------------------
    # PATCH 2 : HIGH-PRESSURE JWL SHOCK DRIVER
    # -------------------------------------------------------------------------

    "patch_icpp(2)%geometry": 3,

    "patch_icpp(2)%alter_patch(1)": "T",

    "patch_icpp(2)%x_centroid":
        float(driver_x_centroid),

    "patch_icpp(2)%y_centroid":
        0.0,

    "patch_icpp(2)%length_x":
        float(driver_length_x),

    "patch_icpp(2)%length_y":
        float(y_len),

    "patch_icpp(2)%vel(1)":
        float(driver_velocity),
    "patch_icpp(2)%vel(2)": 0.0,

    "patch_icpp(2)%pres":
        float(driver_pressure),

    "patch_icpp(2)%alpha_rho(1)":
        float(driver_density),

    "patch_icpp(2)%alpha(1)":
        1.0,

    # -------------------------------------------------------------------------
    # FLUID PROPERTIES
    # -------------------------------------------------------------------------

    "fluid_pp(1)%eos": 2,

    "fluid_pp(1)%Re(1)": 1000.0,

    # Still needed internally
    "fluid_pp(1)%gamma":
        float(1.0 / (gam_a - 1.0)),

    "fluid_pp(1)%pi_inf":
        0.0,

    # JWL constants
    "fluid_pp(1)%jwl_A":
        float(jwl_A),

    "fluid_pp(1)%jwl_B":
        float(jwl_B),

    "fluid_pp(1)%jwl_R1":
        float(jwl_R1),

    "fluid_pp(1)%jwl_R2":
        float(jwl_R2),

    "fluid_pp(1)%jwl_omega":
        float(jwl_omega),

    "fluid_pp(1)%jwl_rho0":
        float(jwl_rho0),
}

# =============================================================================
# 12. PRINT JSON
# =============================================================================

print(json.dumps(case, indent=4))