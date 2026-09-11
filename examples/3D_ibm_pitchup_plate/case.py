#!/usr/bin/env python3
# Canonical pitch-up: a flat plate pitches about its leading edge from 0 to 45 degrees on the smoothed linear
# ramp of the AIAA low-Reynolds-number canonical cases, driven by patch_ib%kin_model = 2. See readme.md.
import json
import math

U, rho, gamma, Ma, Re = 1.0, 1.0, 1.4, 0.2, 300.0
P = rho * U**2 / (gamma * Ma**2)
cs = math.sqrt(gamma * P / rho)

c = 1.0  # chord
semi_span = 2.0 * c  # aspect ratio 4 overall, half of it simulated
th_max = math.radians(45.0)  # the canonical maneuver pitches from 0 to 45 degrees
K = math.pi / 8  # reduced pitch rate: the plate pitches over one chord of travel (case C1)
Omega = 2 * K * U  # nominal pitch rate, rad per c/U
t_p = th_max / Omega  # pitch duration
a_smooth = 21.0  # Eldredge smoothing for C1
t0 = 3.0  # the ramp starts here; before it the plate sits at zero incidence
t_end = t0 + t_p + 4.0

dx = 0.0125 * c
thick = 4 * dx  # at least four cells across the section: see readme.md
x0, x1 = -2.0 * c, 4.0 * c
y0, y1 = 0.0, 2.8 * c  # y = 0 is the symmetry plane at the wing root
z0, z1 = -2.0 * c, 2.0 * c
m, n, p = int((x1 - x0) / dx) - 1, int((y1 - y0) / dx) - 1, int((z1 - z0) / dx) - 1
dt = 0.4 * dx / (U + cs)
nt = int(t_end / dt)

case = {
    "run_time_info": "T",
    "parallel_io": "T",
    "prim_vars_wrt": "T",
    "ib_state_wrt": "T",
    "format": "silo",
    "precision": "double",
    "x_domain%beg": x0,
    "x_domain%end": x1,
    "y_domain%beg": y0,
    "y_domain%end": y1,
    "z_domain%beg": z0,
    "z_domain%end": z1,
    "m": m,
    "n": n,
    "p": p,
    "cyl_coord": "F",
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": nt,
    "t_step_save": max(1, nt // 20),
    "num_patches": 1,
    "num_fluids": 1,
    "model_eqns": "5eq",
    "alt_soundspeed": "F",
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1.0e-10,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": "arithmetic",
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": "hllc",
    "low_Mach": 2,
    "wave_speeds": "direct",
    "viscous": "T",
    "fd_order": 4,
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.5 * (x0 + x1),
    "patch_icpp(1)%y_centroid": 0.5 * (y0 + y1),
    "patch_icpp(1)%z_centroid": 0.5 * (z0 + z1),
    "patch_icpp(1)%length_x": x1 - x0,
    "patch_icpp(1)%length_y": y1 - y0,
    "patch_icpp(1)%length_z": z1 - z0,
    "patch_icpp(1)%vel(1)": U,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": P,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%eos": "ideal_gas",
    "fluid_pp(1)%Re(1)": Re,
    "bc_x%beg": -7,
    "bc_x%grcbc_in": "T",
    "bc_x%vel_in(1)": U,
    "bc_x%vel_in(2)": 0.0,
    "bc_x%vel_in(3)": 0.0,
    "bc_x%pres_in": P,
    "bc_x%alpha_rho_in(1)": rho,
    "bc_x%alpha_in(1)": 1.0,
    "bc_x%end": -8,
    "bc_x%grcbc_out": "T",
    "bc_x%pres_out": P,
    "bc_y%beg": -2,
    "bc_y%end": -9,
    "bc_z%beg": -9,
    "bc_z%end": -9,
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 9,
    # the plate reaches past the symmetry plane so its root is a continuation of the wing, not a wall
    "patch_ib(1)%x_centroid": 0.5 * c,
    "patch_ib(1)%y_centroid": 0.75 * c,
    "patch_ib(1)%z_centroid": 0.0,
    "patch_ib(1)%length_x": c,
    "patch_ib(1)%length_y": semi_span + 0.5 * c,
    "patch_ib(1)%length_z": thick,
    "ib_neighborhood_radius": 3,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%moving_ibm": 1,
    "patch_ib(1)%angles(1)": 0.0,
    "patch_ib(1)%angles(2)": 0.0,
    "patch_ib(1)%angles(3)": 0.0,
    # prescribed kinematics: the Eldredge smoothed linear pitch ramp and hold, about the leading edge
    "patch_ib(1)%kin_model": 2,
    "patch_ib(1)%kin_hinge(1)": 0.0,
    "patch_ib(1)%kin_hinge(2)": 0.0,
    "patch_ib(1)%kin_hinge(3)": 0.0,
    "patch_ib(1)%kin_offset(1)": 0.5 * c,
    "patch_ib(1)%kin_offset(2)": 0.75 * c,
    "patch_ib(1)%kin_offset(3)": 0.0,
    "patch_ib(1)%kin_theta0": th_max,
    "patch_ib(1)%kin_theta_mean": 0.0,
    "patch_ib(1)%kin_pitch_rate": Omega,
    "patch_ib(1)%kin_smooth": a_smooth,
    "patch_ib(1)%kin_t0": t0,
}
print(json.dumps(case, indent=4))
