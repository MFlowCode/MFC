#!/usr/bin/env python3
"""
MFC case generator for a 2D shock interaction with a random packed bed.
Uses MFC's internal particle_bed generator instead of python packing.
"""

import json
import math
import sys

# =============================================================================
# 1. USER INPUTS
# =============================================================================
case_name = "stability_debug_20_particles_2D_mfc_bed"

# Particle / bed controls
dp0 = 5e-3
D = dp0
Grid_Res = 12
target_num_particles = 40
random_seed = 5000
particle_min_spacing = 0.05 * D

# Shock / gas controls
Ms = 3.0
p_ratio_factor = 1.0
Re = 300

# Bed dimensions in units of D
bed_x_start_in_D = 0.0
bed_x_len_in_D = 10.0
bed_y_len_in_D = 10.0

# Collision controls
collision_model = 1
ib_coefficient_of_friction = 0.09
coefficient_of_restitution = 0.7

collision_time_steps = 50.0
collision_time_min_tau_factor = 0.30

# Moving IBM mode
# 0 = fixed particles, 1 = prescribed motion, 2 = two-way coupled motion
moving_ibm_mode = 2

# Flow / numerics
p1 = 101325.0
rho1 = 1.204
u1 = 0.0
gam_a = 1.4
slip_flag = "F"

# Conservative adaptive-CFL target
cfl = 0.5

# Output cadence
points_per_tau = 1.0

# =============================================================================
# 2. DOMAIN SETUP
# =============================================================================
xs = -3.0 * D

bed_x0 = bed_x_start_in_D * D
bed_x1 = (bed_x_start_in_D + bed_x_len_in_D) * D
bed_y0 = -0.5 * bed_y_len_in_D * D
bed_y1 = +0.5 * bed_y_len_in_D * D

bed_lx = bed_x1 - bed_x0
bed_ly = bed_y1 - bed_y0

# Center of the bounding box for the particle bed
bed_x_centroid = 0.5 * (bed_x0 + bed_x1)
bed_y_centroid = 0.5 * (bed_y0 + bed_y1)

x_pad_upstream_D = 4.0
x_pad_downstream_D = 7.0

# Added 1D of padding to the top and bottom so MFC's native particle_bed
# generator doesn't place particles that clip through the boundary.
y_pad_D = 0.0

L_xbeg = bed_x0 - x_pad_upstream_D * D
L_xend = bed_x1 + x_pad_downstream_D * D
L_ybeg = bed_y0 - y_pad_D * D
L_yend = bed_y1 + y_pad_D * D

x_len = L_xend - L_xbeg
y_len = L_yend - L_ybeg

r_ib = 0.5 * dp0
dx = dp0 / Grid_Res
dy = dx

m = max(1, int(round(x_len / dx)))
n = max(1, int(round(y_len / dy)))
p = 0 # 2D Case
total_cells = m * n

# =============================================================================
# 3. SHOCK / POST-SHOCK STATE
# =============================================================================
c1 = math.sqrt(gam_a * p1 / rho1)
us = Ms * c1

base_p_ratio = (2.0 * gam_a * Ms**2 - (gam_a - 1.0)) / (gam_a + 1.0)
rho_ratio = (gam_a + 1.0) * Ms**2 / ((gam_a - 1.0) * Ms**2 + 2.0)
p_ratio = p_ratio_factor * base_p_ratio

p2 = p1 * p_ratio
rho2 = rho1 * rho_ratio
u2 = us * (1.0 - rho1 / rho2)
c2 = math.sqrt(gam_a * p2 / rho2)

t_stop = 1.0 * (L_xend - xs) / us

# =============================================================================
# 4. TIME SCALES & COLLISION TIME
# =============================================================================
tau = D / us
t_save = tau / points_per_tau

max_wave_speed_est = max(c1, abs(u1) + c1, abs(u2) + c2)
dt_cfl_est = cfl * dx / max_wave_speed_est

collision_time_from_dt = collision_time_steps * dt_cfl_est
collision_time_from_tau = collision_time_min_tau_factor * tau
collision_time = max(collision_time_from_dt, collision_time_from_tau)

# =============================================================================
# 5. PARTICLE MASS MODEL (2D AREA)
# =============================================================================
rho_particle = 3 * rho1
particle_area = math.pi * r_ib**2
mass_particle = rho_particle * particle_area


# =============================================================================
# 7. BUILD MFC CASE DICTIONARY
# =============================================================================
case = {
    # Output
    "run_time_info": "T",
    "format": 1,
    "precision": 2,
    "parallel_io": "T",
    "prim_vars_wrt": "F",
    "cons_vars_wrt": "F",
    "E_wrt": "F",
    "ib_state_wrt": "T",
    "rho_wrt": "F",
    "pres_wrt": "T",
    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",
    "alpha_wrt(1)": "F",
    "alpha_rho_wrt(1)": "F",
    "gamma_wrt": "F",
    "pi_inf_wrt": "F",
    "c_wrt": "F",
    
    # Domain
    "x_domain%beg": float(L_xbeg),
    "x_domain%end": float(L_xend),
    "y_domain%beg": float(L_ybeg),
    "y_domain%end": float(L_yend),
    "cyl_coord": "F",
    "m": int(m),
    "n": int(n),
    "p": 0,
    
    # Time Stepping (Adaptive CFL)
    "cfl_adap_dt": "T",
    "cfl_const_dt": "F",
    "cfl_target": float(cfl),
    "n_start": 0,
    "t_stop": float(t_stop),
    "t_save": float(t_save),
    
    # Solver / Numerics
    "num_patches": 2,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "weno_Re_flux": "T",
    "weno_avg": "F",
    "avg_state": 2,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "viscous": "T",
    
    # Boundary Conditions
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -1, # Slip wall condition
    "bc_y%end": -1,
    
    # IBM Global Settings
    "ib": "T",
    "num_ibs": 0, # Set to 0 because we use particle_bed
    "collision_model": int(collision_model),
    "ib_coefficient_of_friction": float(ib_coefficient_of_friction),
    "collision_time": float(collision_time),
    "coefficient_of_restitution": float(coefficient_of_restitution),
    
    # Internal Particle Bed Generator
    "num_particle_beds": 1,
    "particle_bed(1)%x_centroid": float(bed_x_centroid),
    "particle_bed(1)%y_centroid": float(bed_y_centroid),
    "particle_bed(1)%z_centroid": 0.0,
    "particle_bed(1)%length_x": float(bed_lx),
    "particle_bed(1)%length_y": float(bed_ly),
    "particle_bed(1)%length_z": 0.0,
    "particle_bed(1)%num_particles": int(target_num_particles),
    "particle_bed(1)%radius": float(r_ib),
    "particle_bed(1)%mass": float(mass_particle),
    "particle_bed(1)%min_spacing": float(particle_min_spacing),
    "particle_bed(1)%moving_ibm": int(moving_ibm_mode),
    "particle_bed(1)%seed": int(random_seed),
    
    # Patch 1: Post-shock region
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": float((xs + L_xbeg) / 2.0),
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": float(xs - L_xbeg),
    "patch_icpp(1)%length_y": float(y_len),
    "patch_icpp(1)%vel(1)": float(u2),
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": float(p2),
    "patch_icpp(1)%alpha_rho(1)": float(rho2),
    "patch_icpp(1)%alpha(1)": 1.0,
    
    # Patch 2: Pre-shock region
    "patch_icpp(2)%geometry": 3,
    "patch_icpp(2)%x_centroid": float((L_xend + xs) / 2.0),
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%length_x": float(L_xend - xs),
    "patch_icpp(2)%length_y": float(y_len),
    "patch_icpp(2)%vel(1)": float(u1),
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": float(p1),
    "patch_icpp(2)%alpha_rho(1)": float(rho1),
    "patch_icpp(2)%alpha(1)": 1.0,
    
    # Fluid properties
    "fluid_pp(1)%gamma": float(1.0 / (gam_a - 1.0)),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": float(Re),
}

# =============================================================================
# 8. PRINT JSON
# =============================================================================
print(json.dumps(case, indent=4))