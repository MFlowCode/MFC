import json
import math
import sys

"""
A case file meant to perform a direct comparison with:
Canuto D, Taira K. Two-dimensional compressible viscous flow around a circular cylinder. Journal of Fluid Mechanics. 2015;785:349-371. doi:10.1017/jfm.2015.635


This performs a comparison to the low reynolds case in Table 4 of that paper
"""

Re_p = 40.0
Ma_inf = 0.1
d_cyl = 1.0
r_cyl = d_cyl / 2.0

gamma = 1.4
U_inf = 1.0
rho_inf = 1.0
P_inf = rho_inf * U_inf**2 / (gamma * Ma_inf**2)

case = {
    # --- Output ---
    "format": 1,
    "precision": 2,
    "run_time_info": "T",
    "parallel_io": "T",
    "prim_vars_wrt": "T",
    "ib_state_wrt": "T",
    # --- Domain ---
    "x_domain%beg": -15.0,
    "x_domain%end": 15.0,
    "y_domain%beg": -15.0,
    "y_domain%end": 15.0,
    "m": 3000,
    "n": 3000,
    "p": 0,
    "cyl_coord": "F",
    # --- Time stepping ---
    "cfl_adap_dt": "T",
    "cfl_target": 0.5,
    "n_start": 0,
    "t_save": 10.0,
    "t_stop": 100.0,
    # --- Numerics ---
    "num_patches": 1,
    "num_fluids": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-10,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": 2,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 2,
    "low_Mach": 2,  # enable low mach disipation setting
    "wave_speeds": 1,
    "viscous": "T",
    "fd_order": 4,
    # --- Initial condition: uniform freestream, 2D rectangle ---
    "patch_icpp(1)%geometry": 3,  # 2D rectangle
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": 30.0,
    "patch_icpp(1)%length_y": 30.0,
    "patch_icpp(1)%vel(1)": U_inf,  # flow in +x
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": P_inf,
    "patch_icpp(1)%alpha_rho(1)": rho_inf,
    "patch_icpp(1)%alpha(1)": 1.0,
    # --- Fluid properties (calorically perfect air) ---
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),  # MFC: 1/(gamma-1)
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": Re_p / (d_cyl * U_inf),
    # --- Boundary conditions ---
    "bc_x%beg": -7,
    "bc_x%grcbc_in": "T",  # grbc for low-mach case
    "bc_x%vel_in(1)": U_inf,
    "bc_x%vel_in(2)": 0.0,
    "bc_x%pres_in": P_inf,
    "bc_x%alpha_rho_in(1)": rho_inf,
    "bc_x%alpha_in(1)": 1.0,
    "bc_x%end": -8,
    "bc_x%grcbc_out": "T",
    "bc_x%pres_out": P_inf,
    "bc_y%beg": -9,
    "bc_y%end": -9,
    # --- IBM: 2D circle (cylinder cross-section, geometry=2) ---
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": 0.0,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%radius": r_cyl,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%moving_ibm": 0,
}

print(json.dumps(case, indent=4))
