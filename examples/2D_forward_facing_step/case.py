import json
import math

h = 0.2

gam_a = 1.4
p0 = 1
rho0 = 1.4
c0 = math.sqrt(gam_a * p0 / rho0)
v0 = 3 * c0
mu = rho0 * v0 * h / 2e5

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            "x_domain%beg": 0,
            "x_domain%end": 15 * h,
            "y_domain%beg": 0,
            "y_domain%end": 5 * h,
            "cyl_coord": "F",
            "m": 1499,
            "n": 499,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.8,
            "n_start": 0,
            "t_save": 0.04,
            "t_stop": 4,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_avg": "T",
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "avg_state": 2,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -2,
            "bc_y%end": -2,
            "ib": "T",
            "num_ibs": 1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1 Background
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 7.5 * h,
            "patch_icpp(1)%y_centroid": 2.5 * h,
            "patch_icpp(1)%length_x": 15 * h,
            "patch_icpp(1)%length_y": 5 * h,
            "patch_icpp(1)%vel(1)": v0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch: No slip rectangle
            "patch_ib(1)%geometry": 3,
            "patch_ib(1)%x_centroid": 11.5 * h,
            "patch_ib(1)%y_centroid": 0 * h,
            "patch_ib(1)%length_x": 17 * h,
            "patch_ib(1)%length_y": 2 * h,
            "patch_ib(1)%slip": "T",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0 / (gam_a - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1 / mu,
        },
        indent=4,
    )
)
