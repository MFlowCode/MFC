import json
import math

Mu = 1.84e-05
gam_a = 1.4
rho1 = 1.19
D = 3

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -4 * D,
            "x_domain%end": 4.0 * D,
            # y direction
            "y_domain%beg": -2 * D,
            "y_domain%end": 2 * D,
            # z direction
            "z_domain%beg": -2 * D,
            "z_domain%end": 2 * D,
            "cyl_coord": "F",
            "m": 99,
            "n": 49,
            "p": 49,
            "dt": 1.0e-9,
            "t_step_start": 0,
            "t_step_stop": 1000,
            "t_step_save": 10,
            # Simulation Algorithm Parameters
            "num_patches": 1,
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
            "bc_z%beg": -3,
            "bc_z%end": -3,
            "ib": "T",
            "num_ibs": 1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "E_wrt": "T",
            "omega_wrt(1)": "T",
            "omega_wrt(2)": "T",
            "omega_wrt(3)": "T",
            "parallel_io": "T",
            "fd_order": 2,
            # Patch: Constant filled with air
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_x": 100 * D,
            "patch_icpp(1)%length_y": 50 * D,
            "patch_icpp(1)%length_z": 50 * D,
            "patch_icpp(1)%vel(1)": 0.001,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%vel(3)": 0.0e00,
            "patch_icpp(1)%pres": 100000,
            "patch_icpp(1)%alpha_rho(1)": (1.0) * rho1,
            "patch_icpp(1)%alpha(1)": 1.0e00,
            # Patch: Model Immersed Boundary
            "patch_ib(1)%geometry": 12,
            "patch_ib(1)%model_filepath": "Cube_IBM.stl",
            "patch_ib(1)%model_translate(1)": 0,
            "patch_ib(1)%model_translate(2)": -0.05,
            "patch_ib(1)%model_translate(3)": -0.05,
            "patch_ib(1)%model_spc": 100,
            "patch_ib(1)%model_threshold": 0.95,
            "patch_ib(1)%slip": "F",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),
            "fluid_pp(1)%pi_inf": 0,
        }
    )
)
