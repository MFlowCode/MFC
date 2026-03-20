#!/usr/bin/env python3
import json
import math

# 2D axisymmetric hypoelastic case with HLLD solver.
# Solid sphere in liquid hit by focused acoustic pulse.
# Adapted from HLL_new/circle_1_interface for cylindrical coordinates.
# x = axial, y = radial (axis of symmetry at y = 0).

config = {
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": 0,
    "x_domain%end": 1.0,
    "y_domain%beg": 0,
    "y_domain%end": 1.0,
    "cyl_coord": "T",
    "m": 49,
    "n": 49,
    "p": 0,
    "dt": 2.0e-6,
    "t_step_start": 0,
    "t_step_stop": 250,
    "t_step_save": 50,
    # Simulation Algorithm
    "num_patches": 2,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 2,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 1,
    "weno_order": 1,
    "weno_eps": 1.0e-20,
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 4,
    # "hypo_hll_interface_rhs": "T",
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "rho_wrt": "T",
    "parallel_io": "T",
    # Hypoelasticity
    "hypoelasticity": "T",
    "fd_order": 4,
    # Patch 1: Liquid background
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%y_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%length_y": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 1e05,
    "patch_icpp(1)%tau_e(1)": 0.0,
    "patch_icpp(1)%alpha_rho(1)": 1000 * (1.0 - 1e-8),
    "patch_icpp(1)%alpha(1)": 1.0 - 1e-8,
    "patch_icpp(1)%alpha_rho(2)": 1000 * 1e-8,
    "patch_icpp(1)%alpha(2)": 1e-8,
    # Patch 2: Solid sphere (circle, offset from axis)
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0.6,
    "patch_icpp(2)%y_centroid": 0.2,
    "patch_icpp(2)%radius": 0.1,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": 1e05,
    "patch_icpp(2)%tau_e(1)": 0.0,
    "patch_icpp(2)%alpha_rho(1)": 1000 * 1e-8,
    "patch_icpp(2)%alpha(1)": 1e-8,
    "patch_icpp(2)%alpha_rho(2)": 1000 * (1.0 - 1e-8),
    "patch_icpp(2)%alpha(2)": 1.0 - 1e-8,
    # Acoustic source (axisymmetric focused)
    "acoustic_source": "T",
    "num_source": 1,
    "acoustic(1)%support": 6,
    "acoustic(1)%loc(1)": 0.1,
    "acoustic(1)%loc(2)": 0.0,
    "acoustic(1)%pulse": 2,
    "acoustic(1)%npulse": 1,
    "acoustic(1)%mag": 1.0,
    "acoustic(1)%foc_length": 0.8,
    "acoustic(1)%aperture": 0.8,
    "acoustic(1)%gauss_sigma_time": 4e-5,
    "acoustic(1)%delay": 2e-4,
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 4.4e00 * 5.57e08 / (4.4e00 - 1.0e00),
    "fluid_pp(1)%G": 0.0,
    "fluid_pp(2)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
    "fluid_pp(2)%pi_inf": 4.4e00 * 5.57e08 / (4.4e00 - 1.0e00),
    "fluid_pp(2)%G": 1e9,
}

print(json.dumps(config, indent=4))
