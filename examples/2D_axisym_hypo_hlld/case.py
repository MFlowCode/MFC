#!/usr/bin/env python3
import json

# 2D axisymmetric hypoelastic case with HLLD solver.
# Epoxy rectangular inclusion in water hit by focused acoustic pulse.
# Water: gamma=4.4, pi_inf=6.0e8, G=0
# Epoxy: gamma=4.4, pi_inf=2.4e9, rho=1180, G=1.5e9

gamma_w, pi_inf_w, rho_w = 4.4, 6.0e8, 1000.0
gamma_e, pi_inf_e, rho_e = 4.4, 2.4e9, 1180.0

config = {
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "y_domain%beg": 0.0,
    "y_domain%end": 1.0,
    "cyl_coord": "T",
    "m": 99,
    "n": 99,
    "p": 0,
    "dt": 1.5e-6,
    "t_step_start": 0,
    "t_step_stop": 400,
    "t_step_save": 100,
    # Simulation Algorithm
    "num_patches": 2,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 2,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 4,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -6,
    "bc_x%end": -6,
    "bc_y%beg": -2,
    "bc_y%end": -6,
    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "rho_wrt": "T",
    "parallel_io": "T",
    # Hypoelasticity
    "hypoelasticity": "T",
    "fd_order": 4,
    # Patch 1: Water background
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%y_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%length_y": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 1e05,
    "patch_icpp(1)%tau_e(1)": 0.0,
    "patch_icpp(1)%alpha_rho(1)": rho_w * (1.0 - 1e-8),
    "patch_icpp(1)%alpha(1)": 1.0 - 1e-8,
    "patch_icpp(1)%alpha_rho(2)": rho_e * 1e-8,
    "patch_icpp(1)%alpha(2)": 1e-8,
    # Patch 2: Epoxy rectangular inclusion
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%geometry": 3,
    "patch_icpp(2)%x_centroid": 0.6,
    "patch_icpp(2)%y_centroid": 0.2,
    "patch_icpp(2)%length_x": 0.2,
    "patch_icpp(2)%length_y": 0.2,
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 2.0,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": 1e05,
    "patch_icpp(2)%tau_e(1)": 0.0,
    "patch_icpp(2)%alpha_rho(1)": rho_w * 1e-8,
    "patch_icpp(2)%alpha(1)": 1e-8,
    "patch_icpp(2)%alpha_rho(2)": rho_e * (1.0 - 1e-8),
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
    "fluid_pp(1)%gamma": 1.0 / (gamma_w - 1.0),
    "fluid_pp(1)%pi_inf": gamma_w * pi_inf_w / (gamma_w - 1.0),
    "fluid_pp(1)%G": 0.0,
    "fluid_pp(2)%gamma": 1.0 / (gamma_e - 1.0),
    "fluid_pp(2)%pi_inf": gamma_e * pi_inf_e / (gamma_e - 1.0),
    "fluid_pp(2)%G": 1.5e9,
}

print(json.dumps(config, indent=4))
