#!/usr/bin/env python3
import math
import json

Ma = 1.4
ps = 238558
rho_post_a = 2.18
rho_a = 1.204
rho_w = 1000
gam_a = 1.4
gam_w = 6.12
pi_w = 3.43e8
vel = 226
rho = 1
c_l = math.sqrt(1.4 * ps / rho)
eps = 1e-9

D = 0.048
Ny = 299.0
Nx = 1199.0
dx = 0.25 / Nx  # 8.3e-6

time_end = 0.005  # 50us
cfl = 0.25

dt = cfl * dx / c_l  # 5.3E-9
Nt = int(time_end / dt)  # 10000

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "F",
            # Computational Domain Parameters
            "x_domain%beg": -4 * D,
            "x_domain%end": 20 * D,
            "y_domain%beg": 0.0,
            "y_domain%end": 6 * D,
            "stretch_y": "T",
            "a_y": 3.67,
            "y_a": -5.7 * D,
            "y_b": 5.7 * D,
            "loops_y": 2,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": math.ceil(Nt / 100),
            # Simulation Algorithm Parameters
            "num_patches": 3,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "recon_type": 2,
            "muscl_order": 2,
            "muscl_lim": 4,
            "int_comp": "T",
            "null_weights": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -6,  # 11,
            "bc_x%end": -6,  # 12
            "bc_y%beg": -2,
            "bc_y%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: Background
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 8 * D,
            "patch_icpp(1)%y_centroid": 6 * D,
            "patch_icpp(1)%length_x": 24 * D,
            "patch_icpp(1)%length_y": 14 * D,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": 101325.0,
            "patch_icpp(1)%alpha_rho(1)": eps * 1000,
            "patch_icpp(1)%alpha_rho(2)": (1 - eps) * 1.17,
            "patch_icpp(1)%alpha(1)": eps,
            "patch_icpp(1)%alpha(2)": 1 - eps,
            # Patch 2: Shocked state
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -2.5 * D,
            "patch_icpp(2)%y_centroid": 6 * D,
            "patch_icpp(2)%length_x": 3 * D,
            "patch_icpp(2)%length_y": 14 * D,
            "patch_icpp(2)%vel(1)": vel,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%pres": ps,
            "patch_icpp(2)%alpha_rho(1)": eps * 1000,
            "patch_icpp(2)%alpha_rho(2)": (1 - eps) * rho_post_a,
            "patch_icpp(2)%alpha(1)": eps,
            "patch_icpp(2)%alpha(2)": 1 - eps,
            # Patch 3: Bubble
            "patch_icpp(3)%geometry": 2,
            "patch_icpp(3)%x_centroid": 0,
            "patch_icpp(3)%y_centroid": 0,
            "patch_icpp(3)%radius": D / 2,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0e00,
            "patch_icpp(3)%pres": 101325.0,
            "patch_icpp(3)%alpha_rho(1)": (1 - eps) * rho_w,
            "patch_icpp(3)%alpha_rho(2)": eps * 1.17,
            "patch_icpp(3)%alpha(1)": 1 - eps,  # 0.95
            "patch_icpp(3)%alpha(2)": eps,  # 0.05,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_w - 1.0e00),
            "fluid_pp(1)%pi_inf": pi_w * gam_w / (gam_w - 1.0e00),
            "fluid_pp(2)%gamma": 1.0e00 / (gam_a - 1.0e00),
            "fluid_pp(2)%pi_inf": 0.0e00,
        }
    )
)
