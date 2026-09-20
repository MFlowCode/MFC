#!/usr/bin/env python3
"""Immersed-boundary force on a thin plate, against a published measurement (issue 1849).

A 2D flat plate pitching about its leading edge, 0 -> 45 degrees on an Eldredge smoothed ramp (a = 21),
K = pi/8 (case C1), Re_c = 300, Ma 0.2. The plate is 2.5 percent of the chord thick, matching the experiment.

    Jantzen, Taira, Granlund & Ol, Phys. Fluids 26, 053606 (2014), Fig. 10, 2D panel, curve C1.

`NCELL` sets how many cells lie across the plate thickness; the physical problem does not change with it, so
the sweep isolates the immersed-boundary resolution requirement from any change of geometry. See README.md
for what the sweep shows and why it matters.

    NCELL=2 python3 case.py      dx = 0.0125 c
    NCELL=4 python3 case.py      dx = 0.00625 c   (default)
    NCELL=8 python3 case.py      dx = 0.003125 c
    NCELL=16 python3 case.py     dx = 0.0015625 c
"""

import json
import math
import os

U, rho, Ma, gamma, Re = 1.0, 1.0, 0.2, 1.4, 300.0
P = rho * U**2 / (gamma * Ma**2)
cs = math.sqrt(gamma * P / rho)
K = math.pi / 8
Omega = 2 * K * U  # rad per c/U; pitch time = 45 deg / Omega = 1 c/U
th_max = math.radians(45.0)
t_p = th_max / Omega
a_smooth = 21.0
THICK = 0.025
t0 = 2.0  # settle at 0 deg before the ramp
t_end = t0 + t_p + 4.0
x0, x1, y0, y1 = -2.0, 5.0, -2.5, 2.5

# L4 added after the first three failed to converge: refining 2 -> 4 -> 8 cells across the thickness moved
# the peak lift 6.41 -> 4.68 -> 4.46 against a reference of 7.00, i.e. away from it and then stalling.
# Two under-resolved answers landing near each other is not convergence. If a few cells across a thin
# body is simply too few for the immersed boundary, 16 should move back toward the reference; if the
# finite thickness is genuinely the difference, it should stay near 4.5.
LEVELS = {"L1": 0.0125, "L2": 0.00625, "L3": 0.003125, "L4": 0.0015625}


NCELL = int(os.environ.get("NCELL", 4))
dx = THICK / NCELL
m, n = int((x1 - x0) / dx) - 1, int((y1 - y0) / dx) - 1
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
    "m": m,
    "n": n,
    "p": 0,
    "cyl_coord": "F",
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": nt,
    "t_step_save": max(1, nt // 40),
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
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * (x0 + x1),
    "patch_icpp(1)%y_centroid": 0.5 * (y0 + y1),
    "patch_icpp(1)%length_x": x1 - x0,
    "patch_icpp(1)%length_y": y1 - y0,
    "patch_icpp(1)%vel(1)": U,
    "patch_icpp(1)%vel(2)": 0.0,
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
    "bc_x%pres_in": P,
    "bc_x%alpha_rho_in(1)": rho,
    "bc_x%alpha_in(1)": 1.0,
    "bc_x%end": -8,
    "bc_x%grcbc_out": "T",
    "bc_x%pres_out": P,
    "bc_y%beg": -9,
    "bc_y%end": -9,
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 3,
    "patch_ib(1)%x_centroid": 0.5,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%length_x": 1.0,
    "patch_ib(1)%length_y": THICK,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%moving_ibm": 1,
    "patch_ib(1)%angles(3)": 0.0,
    "patch_ib(1)%angular_vel(3)": 0.0,
    "omega_wrt(3)": "T",
}
tau = f"(t - {t0})"
th = f"(0.5*{th_max}*(1.0 + (log(cosh({a_smooth}*{tau})) - log(cosh({a_smooth}*({tau} - {t_p}))))/{a_smooth * t_p}))"
thd = f"(0.5*{Omega}*(tanh({a_smooth}*{tau}) - tanh({a_smooth}*({tau} - {t_p}))))"
case["patch_ib(1)%angular_vel(3)"] = thd
case["patch_ib(1)%vel(1)"] = f"-0.5*{thd}*sin({th})"
case["patch_ib(1)%vel(2)"] = f"0.5*{thd}*cos({th})"

if __name__ == "__main__":
    if os.environ.get("SUMMARY"):
        print(f"{NCELL} cells across the {THICK:g} c thickness: dx = {dx:g} c, " f"{m + 1} x {n + 1} = {(m + 1) * (n + 1) / 1e6:.2f} M cells, {nt} steps of dt = {dt:.2e}")
    else:
        print(json.dumps(case, indent=4))
