#!/usr/bin/env python3
"""
2D Minimal Case Template
------------------------
A minimal 2D case with a circular perturbation.

Usage:
    ./mfc.sh run case.py
"""
import math
import json

# =============================================================================
# SIMULATION PARAMETERS - Modify these for your case
# =============================================================================

# Grid resolution
Nx = 99                     # Cells in x-direction
Ny = 99                     # Cells in y-direction

# Domain size
x_start, x_end = 0.0, 1.0
y_start, y_end = 0.0, 1.0

# Time stepping
dt = 1.0e-6
Nt = 1000

# Background state
rho_bg = 1.0
vel_x_bg = 0.0
vel_y_bg = 0.0
pres_bg = 1.0e5

# Perturbation (circular region)
x_center = 0.5
y_center = 0.5
radius = 0.1
rho_pert = 2.0
pres_pert = 2.0e5

# Fluid properties
gamma = 1.4

# =============================================================================
# CASE DICTIONARY - MFC configuration
# =============================================================================
print(json.dumps({
    # Logistics
    "run_time_info": "T",

    # Computational Domain
    "x_domain%beg": x_start,
    "x_domain%end": x_end,
    "y_domain%beg": y_start,
    "y_domain%end": y_end,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": max(1, Nt // 10),

    # Simulation Algorithm
    "num_patches": 2,
    "model_eqns": 2,
    "num_fluids": 1,
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,

    # Boundary Conditions
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,

    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",

    # Patch 1: Background
    "patch_icpp(1)%geometry": 3,  # Rectangle
    "patch_icpp(1)%x_centroid": (x_start + x_end) / 2,
    "patch_icpp(1)%y_centroid": (y_start + y_end) / 2,
    "patch_icpp(1)%length_x": x_end - x_start,
    "patch_icpp(1)%length_y": y_end - y_start,
    "patch_icpp(1)%vel(1)": vel_x_bg,
    "patch_icpp(1)%vel(2)": vel_y_bg,
    "patch_icpp(1)%pres": pres_bg,
    "patch_icpp(1)%alpha_rho(1)": rho_bg,
    "patch_icpp(1)%alpha(1)": 1.0,

    # Patch 2: Circular perturbation
    "patch_icpp(2)%geometry": 2,  # Circle
    "patch_icpp(2)%x_centroid": x_center,
    "patch_icpp(2)%y_centroid": y_center,
    "patch_icpp(2)%radius": radius,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%vel(1)": vel_x_bg,
    "patch_icpp(2)%vel(2)": vel_y_bg,
    "patch_icpp(2)%pres": pres_pert,
    "patch_icpp(2)%alpha_rho(1)": rho_pert,
    "patch_icpp(2)%alpha(1)": 1.0,

    # Fluid Properties
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}))
