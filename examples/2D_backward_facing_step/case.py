import json
import math

h = 1
Re = 5100

gam_a = 1.4
p0 = 1
rho0 = 1
c0 = math.sqrt(gam_a * p0 / rho0)
v0 = 0.5 * c0
mu = v0 * h / Re

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            "x_domain%beg": 0,
            "x_domain%end": 30 * h,
            "y_domain%beg": 0,
            "y_domain%end": 6 * h,
            "cyl_coord": "F",
            "m": 499,
            "n": 99,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.8,
            "n_start": 0,
            "t_save": 150 / (100 * v0),
            "t_stop": 150 / v0,
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
            "weno_avg": "T",
            "avg_state": 2,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "bc_x%beg": -7,
            "bc_x%end": -8,
            "bc_y%beg": -16,
            "bc_y%end": -6,
            "ib": "T",
            "num_ibs": 1,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1 Background
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 15 * h,
            "patch_icpp(1)%y_centroid": 3 * h,
            "patch_icpp(1)%length_x": 30 * h,
            "patch_icpp(1)%length_y": 6 * h,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2 Velocity
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 4.5 * h,
            "patch_icpp(2)%y_centroid": 3 * h,
            "patch_icpp(2)%length_x": 9 * h,
            "patch_icpp(2)%length_y": 6 * h,
            "patch_icpp(2)%vel(1)": v0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": p0,
            "patch_icpp(2)%alpha_rho(1)": rho0,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Patch: Rectangle
            "patch_ib(1)%geometry": 3,
            "patch_ib(1)%x_centroid": 4 * h,
            "patch_ib(1)%y_centroid": 0 * h,
            "patch_ib(1)%length_x": 12 * h,
            "patch_ib(1)%length_y": 2 * h,
            "patch_ib(1)%slip": "F",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0 / (gam_a - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": 1 / mu,
        }
    )
)
