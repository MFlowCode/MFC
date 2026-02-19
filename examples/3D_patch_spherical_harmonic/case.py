#!/usr/bin/env python3
"""Minimal 3D acoustic case with geometry 14 (spherical harmonic surface)."""
import math
import json

L = 8.0
N = 50
dx = L / (N + 1)
p_inf = 101325.0
rho_inf = 1.0
gam = 1.4
# Patch 2: same pressure, higher sound speed (different gas) -> lower density
c_ratio = 1.5
rho2 = rho_inf / (c_ratio * c_ratio)
c = math.sqrt(gam * p_inf / rho_inf)
cfl = 0.3
mydt = cfl * dx / c
Tfinal = 10.0 / c
Nt = int(Tfinal / mydt)

# Acoustic source (similar to 3D_acoustic_support7)
t_cross = L / c # domain crossing time
gauss_sigma_time = 0.03 * t_cross
acoustic_delay = 0.15 * t_cross
acoustic_loc_x = -3.5
acoustic_loc_y = 0.0
acoustic_loc_z = 0.0
acoustic_aperture = 6.0
acoustic_foc_length = 3.5

config = {
    "run_time_info": "T",
    "x_domain%beg": -4.0,
    "x_domain%end": 4.0,
    "y_domain%beg": -4.0,
    "y_domain%end": 4.0,
    "z_domain%beg": -4.0,
    "z_domain%end": 4.0,
    "m": N,
    "n": N,
    "p": N,
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": max(1, int(Nt / 10)),
    "num_patches": 2,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 2,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -8,
    "bc_x%end": -8,
    "bc_y%beg": -8,
    "bc_y%end": -8,
    "bc_z%beg": -8,
    "bc_z%end": -8,
    "bc_x%grcbc_in": "F",
    "bc_x%grcbc_out": "T",
    "bc_x%grcbc_vel_out": "F",
    "bc_x%pres_out": p_inf,
    "bc_y%grcbc_in": "F",
    "bc_y%grcbc_out": "T",
    "bc_y%grcbc_vel_out": "F",
    "bc_y%pres_out": p_inf,
    "bc_z%grcbc_in": "F",
    "bc_z%grcbc_out": "T",
    "bc_z%grcbc_vel_out": "F",
    "bc_z%pres_out": p_inf,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%length_y": L,
    "patch_icpp(1)%length_z": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p_inf,
    "patch_icpp(1)%alpha_rho(1)": rho_inf,
    "patch_icpp(1)%alpha_rho(2)": 0.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha(2)": 0.0,
    "patch_icpp(2)%geometry": 14,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": 1.0,
    "patch_icpp(2)%sph_har_coeff(1,0)": 0.5,
    "patch_icpp(2)%sph_har_coeff(1,1)": 0.35,
    "patch_icpp(2)%sph_har_coeff(1,-1)": -0.25,
    "patch_icpp(2)%sph_har_coeff(2,0)": -0.4,
    "patch_icpp(2)%sph_har_coeff(2,1)": 0.35,
    "patch_icpp(2)%sph_har_coeff(2,-1)": -0.2,
    "patch_icpp(2)%sph_har_coeff(2,2)": 0.45,
    "patch_icpp(2)%sph_har_coeff(2,-2)": -0.3,
    "patch_icpp(2)%sph_har_coeff(3,0)": 0.3,
    "patch_icpp(2)%sph_har_coeff(3,1)": -0.25,
    "patch_icpp(2)%sph_har_coeff(3,2)": 0.35,
    "patch_icpp(2)%sph_har_coeff(3,3)": -0.2,
    "patch_icpp(2)%sph_har_coeff(4,0)": -0.2,
    "patch_icpp(2)%sph_har_coeff(4,2)": 0.25,
    "patch_icpp(2)%sph_har_coeff(4,4)": 0.2,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p_inf,
    "patch_icpp(2)%alpha_rho(1)": 0.0,
    "patch_icpp(2)%alpha_rho(2)": rho2,
    "patch_icpp(2)%alpha(1)": 0.0,
    "patch_icpp(2)%alpha(2)": 1.0,
    "patch_icpp(2)%alter_patch(1)": "T",
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    # Acoustic source (similar to 3D_acoustic_support7)
    "acoustic_source": "T",
    "num_source": 1,
    "acoustic(1)%support": 7,
    "acoustic(1)%loc(1)": acoustic_loc_x,
    "acoustic(1)%loc(2)": acoustic_loc_y,
    "acoustic(1)%loc(3)": acoustic_loc_z,
    "acoustic(1)%pulse": 2,
    "acoustic(1)%npulse": 1,
    "acoustic(1)%mag": 1.0,
    "acoustic(1)%gauss_sigma_time": gauss_sigma_time,
    "acoustic(1)%foc_length": acoustic_foc_length,
    "acoustic(1)%aperture": acoustic_aperture,
    "acoustic(1)%delay": acoustic_delay,
}
print(json.dumps(config))
