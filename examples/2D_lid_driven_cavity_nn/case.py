#!/usr/bin/env python3
"""
2D Lid-Driven Cavity Flow with Herschel-Bulkley Non-Newtonian Fluid

Re = 5000, n = 1.5 (shear-thickening) with mesh stretching near walls.

HB model: mu = (tau0/gdot)*(1 - exp(-m*gdot)) + K * gdot^(n-1)
Re_gen = rho * U^(2-n) * L^n / K = 1 * 1^0.5 * 1^1.5 / 0.0002 = 5000

Mesh stretching: cosh-based clustering near all 4 walls (x_a, x_b, y_a, y_b).
"""
import json

eps = 1e-6

# HB model parameters
tau0 = 0.0          # Yield stress (set to 0 for power-law fluid)
K = 0.0002          # Consistency index (Re=5000: K = 1/5000)
nn = 1.5            # Flow behavior index (shear-thickening)
mu_min = 0.00002    # K * gdot_min^(n-1) = 0.0002 * (0.01)^0.5
mu_max = 0.0632     # K * gdot_max^(n-1) = 0.0002 * (1e5)^0.5
hb_m = 1000.0       # Papanastasiou regularization parameter
mu_bulk = 0.0

lid_velocity = 1.0  # Lid velocity (m/s)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "m": 255,
            "n": 255,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.5,
            "n_start": 0,
            "t_stop": 50.0,
            "t_save": 25.0,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1e-16,
            "mapped_weno": "T",
            "weno_Re_flux": "T",
            "mp_weno": "T",
            "weno_avg": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -16,
            "bc_x%end": -16,
            "bc_y%beg": -16,
            "bc_y%end": -16,
            "bc_y%ve1": lid_velocity,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "omega_wrt(3)": "T",
            "fd_order": 4,
            "parallel_io": "T",
            # Patch 1: Base
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1e3,
            "patch_icpp(1)%alpha_rho(1)": 0.5,
            "patch_icpp(1)%alpha(1)": 0.5,
            "patch_icpp(1)%alpha_rho(2)": 0.5,
            "patch_icpp(1)%alpha(2)": 0.5,
            # Fluids Physical Parameters
            # Fluid 1:
            "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": 1.0 / K,
            "fluid_pp(1)%non_newtonian": "T",
            "fluid_pp(1)%tau0": tau0,
            "fluid_pp(1)%K": K,
            "fluid_pp(1)%nn": nn,
            "fluid_pp(1)%mu_max": mu_max,
            "fluid_pp(1)%mu_min": mu_min,
            "fluid_pp(1)%mu_bulk": mu_bulk,
            "fluid_pp(1)%hb_m": hb_m,
            # Fluid 2:
            "fluid_pp(2)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%Re(1)": 1.0 / K,
            "fluid_pp(2)%non_newtonian": "T",
            "fluid_pp(2)%tau0": tau0,
            "fluid_pp(2)%K": K,
            "fluid_pp(2)%nn": nn,
            "fluid_pp(2)%mu_max": mu_max,
            "fluid_pp(2)%mu_min": mu_min,
            "fluid_pp(2)%mu_bulk": mu_bulk,
            "fluid_pp(2)%hb_m": hb_m,
        }
    )
)
