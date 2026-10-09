#!/usr/bin/env python3
"""
NACA 0012 at M = 0.3, alpha = 2 deg, 200 cells per chord: the 2D_ibm_airfoil_surface_pressure example
switched to a viscous no-slip wall (Re = 1000 per chord), with frequent snapshots for time averaging.
Env: RE, T_STOP, T_SAVE, SECOND_ORDER=1 (PR branch only: second-order IB velocities).
"""

import json
import math
import os

Ma = 0.3
alpha_deg = 2.0
gamma = 1.4
rho_inf, U_inf, chord = 1.0, 1.0, 1.0
Re = float(os.environ.get("RE", 1000.0))
P_inf = rho_inf * U_inf**2 / (gamma * Ma**2)

case = {
  # --- Output ---
  "run_time_info": "T",
  "format": 2,
  "precision": 2,
  "parallel_io": "T",
  "prim_vars_wrt": "T",
    # --- Domain: pre-stretch extents; the stretching maps them outward ---
  "x_domain%beg": -3.0,
  "x_domain%end": 4.0,
  "y_domain%beg": -3.0,
  "y_domain%end": 3.0,
  "m": 1399,
  "n": 1199,
  "p": 0,
  "cyl_coord": "F",
  "stretch_x": "T",
  "a_x": 15.0,
  "x_a": -0.8,
  "x_b": 1.8,
  "loops_x": 2,
  "stretch_y": "T",
  "a_y": 15.0,
  "y_a": -0.7,
  "y_b": 0.7,
  "loops_y": 2,
  # --- Time stepping ---
  "cfl_adap_dt": "T",
  "cfl_target": 0.5,
  "n_start": 0,
  "t_save": float(os.environ.get("T_SAVE", 0.2)),
  "t_stop": float(os.environ.get("T_STOP", 6.0)),
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
  "weno_Re_flux": "F",
  "weno_avg": "T",
  "avg_state": 2,
  "mapped_weno": "T",
  "null_weights": "F",
  "mp_weno": "F",
  "riemann_solver": 2,
  "low_Mach": 2,
  "wave_speeds": 1,
  "viscous": "T",
  "fd_order": 4,
  # --- Uniform freestream ---
  "patch_icpp(1)%geometry": 3,
  "patch_icpp(1)%x_centroid": 0.0,
  "patch_icpp(1)%y_centroid": 0.0,
  "patch_icpp(1)%length_x": 1.0e3,
  "patch_icpp(1)%length_y": 1.0e3,
  "patch_icpp(1)%vel(1)": U_inf,
  "patch_icpp(1)%vel(2)": 0.0,
  "patch_icpp(1)%pres": P_inf,
  "patch_icpp(1)%alpha_rho(1)": rho_inf,
  "patch_icpp(1)%alpha(1)": 1.0,
  "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
  "fluid_pp(1)%eos": "ideal_gas",
  "fluid_pp(1)%Re(1)": Re,
  # --- Characteristic far-field boundaries ---
  "bc_x%beg": -7,
  "bc_x%grcbc_in": "T",
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
  # --- IB: NACA 0012 (m must be > 0, so a negligible camber stands in for 0) ---
  "ib": "T",
  "num_ibs": 1,
  "patch_ib(1)%geometry": 4,
  "patch_ib(1)%x_centroid": 0.0,
  "patch_ib(1)%y_centroid": 0.0,
  "patch_ib(1)%airfoil_id": 1,
  "patch_ib(1)%angles(3)": -math.radians(alpha_deg),
  "patch_ib(1)%slip": "F",
  "patch_ib(1)%moving_ibm": 0,
  "ib_airfoil(1)%c": chord,
  "ib_airfoil(1)%t": 0.12,
  "ib_airfoil(1)%p": 0.4,
  "ib_airfoil(1)%m": 1.0e-9,
}

# Second-order IB velocity correction: only on origin/mittal-second-order-ibm-velocities (master rejects these keys)
if os.environ.get("SECOND_ORDER", "0") == "1":
  case.update({"ib_second_order_vel": "T", "ib_ip_min_dist": 1.5})

case.update(json.loads(os.environ.get("EXTRA", "{}")))

print(json.dumps(case, indent=4))