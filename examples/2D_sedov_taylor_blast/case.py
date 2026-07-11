#!/usr/bin/env python3
# Sedov-Taylor cylindrical blast wave -- a quantitative gas-dynamics benchmark.
#
# A concentrated deposit of internal energy in an ideal gas drives a strong blast
# whose shock radius follows the exact self-similar law
#
#     R_shock(t) = xi0 * (E_L / rho0)**(1/4) * t**(1/2)          (cylindrical, nu=2)
#
# with xi0 a known gamma-dependent constant (~1.0 for gamma=1.4). The t**(1/2) power
# law and the strong-shock density ratio (gamma+1)/(gamma-1) = 6 are parameter-free
# consequences of the Euler equations -- this is THE standard check that a compressible
# solver gets blast energetics right, independent of any EOS calibration.
#
# One quadrant [0,L]^2 with reflective axes represents the full line-symmetric blast
# (energy per unit length E_L = 4 * energy in the quadrant deposit). The deposit is a
# small quarter-disk of radius r0 at the origin at high pressure; everything else is
# air at rest. analyze_sedov.py extracts R_shock(t) and compares to the exact solution.
import json
import math

# --- Ideal-gas air (gamma = 1.4) ---
gamma = 1.4
gamma_mfc = 1.0 / (gamma - 1.0)  # = 2.5
rho0 = 1.225
# Low ambient pressure -> high blast Mach number across the whole window, i.e. the
# genuine strong-shock regime where Sedov-Taylor is exact (density ratio -> 6).
p_amb = 1.0e3
cv = 717.5

# --- Energy deposit (point-like blast driver) ---
r0 = 0.03  # deposit radius (kept << measurement radii so Sedov asymptotics hold)
p_hot = 2.0e8
# Excess internal energy per unit length of the full line source (4 quadrants):
#   E_L = 4 * (p_hot - p_amb)/(gamma-1) * (pi r0^2 / 4) = (p_hot - p_amb)/(gamma-1) * pi r0^2
E_L = (p_hot - p_amb) / (gamma - 1.0) * math.pi * r0**2  # ~ 1.41e6 J/m

L = 1.5

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "y_domain%beg": 0.0,
            "y_domain%end": L,
            "m": 299,
            "n": 299,
            "p": 0,
            "dt": 1.5e-7,
            "t_step_start": 0,
            "t_step_stop": 6200,
            "t_step_save": 200,
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            # Reflective symmetry on the two axes; non-reflecting on the outer edges.
            "bc_x%beg": -2,
            "bc_x%end": -3,
            "bc_y%beg": -2,
            "bc_y%end": -3,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: ambient air at rest.
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%y_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p_amb,
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: concentrated energy deposit (quarter-disk at the origin).
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%radius": r0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": p_hot,
            "patch_icpp(2)%alpha_rho(1)": rho0,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Single ideal-gas fluid.
            "fluid_pp(1)%gamma": gamma_mfc,
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%cv": cv,
        }
    )
)
