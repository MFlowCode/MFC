import json
import math

h = 0.2

# Radius as a percentage of height (h)
rc = 0.2

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
            "cfl_target": 0.6,
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
            "igr": "T",
            "igr_pres_lim": "T",
            "igr_order": 3,
            "igr_iter_solver": 1,
            "num_igr_iters": 5,
            "num_igr_warm_start_iters": 50,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -2,
            "bc_y%end": -2,
            "ib": "T",
            "num_ibs": 3,
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
            # Patch: No slip rectangle with rouded corner
            "patch_ib(1)%geometry": 3,
            "patch_ib(1)%x_centroid": 11.5 * h + rc*h,
            "patch_ib(1)%y_centroid": 0 * h,
            "patch_ib(1)%length_x": 17 * h,
            "patch_ib(1)%length_y": 2 * h,
            "patch_ib(1)%slip": "T",
            "patch_ib(2)%geometry": 3,
            "patch_ib(2)%x_centroid": (3 + rc/2)*h,
            "patch_ib(2)%y_centroid": -rc * h,
            "patch_ib(2)%length_x": rc*h,
            "patch_ib(2)%length_y": 2 * h,
            "patch_ib(2)%slip": "T",
            "patch_ib(3)%geometry": 2,
            "patch_ib(3)%x_centroid": (3 + rc)*h,
            "patch_ib(3)%y_centroid": (1 - rc)*h,
            "patch_ib(3)%radius" : rc*h,
            "patch_ib(3)%slip": "T",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0 / (gam_a - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1 / mu,
        },
        indent=4,
    )
)
