#!/usr/bin/env python3
import json, math

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics ================================================================
            "run_time_info": "T",
            # Computational Domain Parameters ==========================================
            "x_domain%beg": 0,
            "x_domain%end": 0.001,
            "y_domain%beg": 0,
            "y_domain%end": 0.001,
            "m": 199,
            "n": 199,
            "p": 0,
            "dt": 2e-9,
            "t_step_start": 0,
            "t_step_stop": 400,
            "t_step_save": 10,
            # Formatted Database Files Structure Parameters ============================
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1 Liquid ===========================================================
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0005,
            "patch_icpp(1)%y_centroid": 0.0005,
            "patch_icpp(1)%length_x": 0.001,
            "patch_icpp(1)%length_y": 0.001,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1e05,
            "patch_icpp(1)%alpha_rho(1)": 1100,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Acoustic source ==========================================================
            "acoustic_source": "T",
            "num_source": 1,
            "acoustic(1)%support": 9,
            "acoustic(1)%loc(1)": 0.00006,
            "acoustic(1)%loc(2)": 0.0005,
            "acoustic(1)%pulse": 2,
            "acoustic(1)%npulse": 1,
            "acoustic(1)%mag": 1.0,
            "acoustic(1)%gauss_sigma_time": 2e-8,
            "acoustic(1)%foc_length": 0.00054,
            "acoustic(1)%aperture": 0.0008,
            "acoustic(1)%num_elements": 4,
            "acoustic(1)%element_on": 0,
            "acoustic(1)%element_spacing_angle": math.pi / 180 * 5,
            "acoustic(1)%delay": 1e-7,
            # Fluids Physical Parameters ===============================================
            "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 4.4e00 * 5.57e08 / (4.4e00 - 1.0e00),
        }
    )
)
