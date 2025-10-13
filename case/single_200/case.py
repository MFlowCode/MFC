#!/usr/bin/env python3
import json
import math

# Calculate 10 evenly spaced probe points on a circle.
center_x = 0.5
center_y = 0.5
radii = [0.45, 0.475]
num_probes = 10
num_layers = 2
probes = {}
for layer in range(num_layers):
    radius = radii[layer]
    for i in range(num_probes):
        angle = 2 * math.pi * i / num_probes
        probes[f"probe({layer*num_probes + i + 1})%x"] = center_x + radius * math.cos(angle)
        probes[f"probe({layer*num_probes + i + 1})%y"] = center_y + radius * math.sin(angle)

config = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    "x_domain%beg": 0,
    "x_domain%end": 1.,
    "y_domain%beg": 0,
    "y_domain%end": 1.,
    "m": 199,
    "n": 199,
    "p": 0,
    "dt": 2e-6,
    "t_step_start": 0,
    # "t_step_stop": 1,
    # "t_step_save": 1,
    "t_step_stop": 1000,
    "t_step_save": 10,
    # Simulation Algorithm Parameters
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
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -6,
    "bc_x%end": -6,
    "bc_y%beg": -6,
    "bc_y%end": -6,
    # Formatted Database Files Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "rho_wrt": "T",
    "parallel_io": "T",
    # Hypoelasticity
    "hypoelasticity": "T",
    "fd_order": 4,
    # Patch 1 Liquid
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%y_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.,
    "patch_icpp(1)%length_y": 1.,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 1e05,
    "patch_icpp(1)%tau_e(1)": 0.0,
    "patch_icpp(1)%alpha_rho(1)": 1000*(1.0 - 1e-8),
    "patch_icpp(1)%alpha(1)": 1.0 - 1e-8,
    "patch_icpp(1)%alpha_rho(2)": 1000*1e-8,
    "patch_icpp(1)%alpha(2)": 1e-8,
    # Patch 1 Solid
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0.6,
    "patch_icpp(2)%y_centroid": 0.61,
    "patch_icpp(2)%radius": 0.15,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": 1e05,
    "patch_icpp(2)%tau_e(1)": 0.0,
    "patch_icpp(2)%alpha_rho(1)": 1000*1e-8,
    "patch_icpp(2)%alpha(1)": 1e-8,
    "patch_icpp(2)%alpha_rho(2)": 1000*(1.0 - 1e-8),
    "patch_icpp(2)%alpha(2)": 1.0 - 1e-8,
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 2.0e00,
    # Acoustic source
    "acoustic_source": "T",
    "num_source": 1,
    "acoustic(1)%support": 5,
    "acoustic(1)%loc(1)": 0.1,
    "acoustic(1)%loc(2)": 0.5,
    "acoustic(1)%pulse": 2,
    "acoustic(1)%npulse": 1,
    "acoustic(1)%dir": 1.0,
    "acoustic(1)%mag": 1.0,
    "acoustic(1)%foc_length": 0.4,
    "acoustic(1)%aperture": 0.75,
    "acoustic(1)%gauss_sigma_time": 5e-5,
    "acoustic(1)%delay": 3e-4,
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 4.4e00 * 5.57e08 / (4.4e00 - 1.0e00),
    "fluid_pp(1)%G": 0.0,
    "fluid_pp(2)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
    "fluid_pp(2)%pi_inf": 4.4e00 * 5.57e08 / (4.4e00 - 1.0e00),
    "fluid_pp(2)%G": 1e9,
    # # Probes
    # "probe_wrt": "T",
    # "num_probes": num_probes*num_layers,
}

# # Add probe coordinates to the config
# config.update(probes)

print(json.dumps(config, indent=4))
