#!/usr/bin/env python3
import json

rho1 = 1.19
gam_a = 1.4
D = 3

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -4 * D,
            "x_domain%end": 4 * D,
            "y_domain%beg": -2 * D,
            "y_domain%end": 2 * D,
            "cyl_coord": "F",
            "m": 49,
            "n": 24,
            "p": 0,
            "dt": 1.0e-9,
            "t_step_start": 0,
            "t_step_stop": 1000,
            "t_step_save": 10,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "T",
            "avg_state": 2,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "viscous": "F",
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "E_wrt": "T",
            "parallel_io": "T",
            "fd_order": 2,
            # Patch 1: domain filled with air (2D rectangle)
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 100 * D,
            "patch_icpp(1)%length_y": 50 * D,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 100000,
            "patch_icpp(1)%alpha_rho(1)": rho1,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: STL/OBJ model imported as a constant-IC patch (geometry 21),
            # resolved via model_id -> stl_models; denser than air so it shows in the field.
            "patch_icpp(2)%geometry": 21,
            "patch_icpp(2)%model_id": 1,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%smoothen": "F",
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 100000,
            "patch_icpp(2)%alpha_rho(1)": 5.0 * rho1,
            "patch_icpp(2)%alpha(1)": 1.0,
            # STL model referenced by patch 2 (native circle ~0.1 dia -> scale to a centered ~5-unit disk)
            "num_stl_models": 1,
            "stl_models(1)%model_filepath": "Circle_IBM.stl",
            "stl_models(1)%model_scale(1)": 50.0,
            "stl_models(1)%model_scale(2)": 50.0,
            "stl_models(1)%model_scale(3)": 50.0,
            "stl_models(1)%model_translate(1)": 0.0,
            "stl_models(1)%model_translate(2)": 0.0,
            "stl_models(1)%model_translate(3)": 0.0,
            "stl_models(1)%model_threshold": 0.95,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),
            "fluid_pp(1)%pi_inf": 0,
        }
    )
)
