#!/usr/bin/env python3
"""
2D periodic shock + transverse shear (single fluid, moderate Mach) — TRANSVERSE-MOMENTUM
validation for the split-explicit low-Mach two-tier acoustic substep (acoustic_substepping).

Why this case exists
--------------------
The 1D shock-stability gate (examples/1D_shock_periodic_lowmach) only exercises the
FACE-NORMAL convective momentum of the robust (full-HLLC) tier. The robust tier also maps
the live HLLC rotated-flux TRANSVERSE slot (d != i) back into the stored global-component
frozen flux; that path is correct by inspection but was empirically unexercised. This case
forces it: an x-normal Sod-like driver launches a right-running shock through a field that
carries a TRANSVERSE velocity component (vel(2)) which JUMPS across the driver interface.
At every flagged x-face the transverse momentum rho*u*v is then nonzero and varies, so the
d=2 (transverse) flux slot of the live HLLC flux is actively transported.

Stable-BC by construction: PERIODIC in both directions (no outflow boundary, so the
separate centered-smooth-tier outflow instability cannot occur). A uniform background mean
flow u_bg bounds the advective outer step (~dx/|u|). The field is uniform in y, so the
exact transverse-momentum behaviour is a clean diagnostic: with no y-gradients and no
transverse pressure gradient, rho*v is advected as a contact quantity and must stay uniform
in y. Any frame/sign error in the d != i mapping would break y-uniformity or disagree with
the full-HLLC baseline.

Toggle via MFC_SPLIT env var:  MFC_SPLIT=T -> split, otherwise -> baseline full HLLC.
"""

import json
import math
import os

split = os.environ.get("MFC_SPLIT", "F").upper()
acoustic_substepping = "T" if split == "T" else "F"

gamma = 1.4

# Sod-like driver vs ambient
p_drive, rho_drive = 1.0, 1.0
p_amb, rho_amb = 0.1, 0.125

c_amb = math.sqrt(gamma * p_amb / rho_amb)
u_bg = 0.4 * c_amb  # background mean flow in x (ambient Mach ~0.4)
# Transverse velocity JUMPS across the driver interface (shear): this is what makes the
# transverse momentum nonzero and varying at the flagged x-faces.
v_amb = 0.15 * c_amb
v_drive = -0.15 * c_amb

L = 1.0
Nx = 199  # 200 cells in x
Ny = 31  # 32 cells in y (uniform-in-y field; >= num_stcls_min*weno_order cells per dim)
dx = L / (Nx + 1)

CFL = 0.4
umag = math.sqrt(u_bg**2 + max(abs(v_amb), abs(v_drive)) ** 2)
dt_init = CFL * dx / (umag + math.sqrt(gamma * p_drive / rho_drive))

T_stop = 0.06
t_save = T_stop / 6.0

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "y_domain%beg": 0.0,
            "y_domain%end": L * (Ny + 1) / (Nx + 1),
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": dt_init,
            "cfl_adap_dt": "T",
            "cfl_target": CFL,
            "t_step_start": 0,
            "n_start": 0,
            "t_stop": T_stop,
            "t_save": t_save,
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 1,
            "alt_soundspeed": "F",
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": "rk3",
            "weno_order": 5,
            "weno_eps": 1.0e-6,
            "mp_weno": "F",
            "weno_avg": "F",
            "mapped_weno": "F",
            "null_weights": "F",
            "riemann_solver": "hllc",
            "wave_speeds": "direct",
            "avg_state": "arithmetic",
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "acoustic_substepping": acoustic_substepping,
            "n_acoustic_substeps": 0,
            "acoustic_div_damp": 0.1,
            "format": "binary",
            "precision": "double",
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: ambient filling the domain (transverse velocity v_amb)
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%y_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": 2.0 * L,
            "patch_icpp(1)%vel(1)": u_bg,
            "patch_icpp(1)%vel(2)": v_amb,
            "patch_icpp(1)%pres": p_amb,
            "patch_icpp(1)%alpha_rho(1)": rho_amb,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: high-pressure driver band (0.15 < x < 0.35), transverse velocity v_drive
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.25 * L,
            "patch_icpp(2)%y_centroid": 0.5 * L,
            "patch_icpp(2)%length_x": 0.2 * L,
            "patch_icpp(2)%length_y": 2.0 * L,
            "patch_icpp(2)%vel(1)": u_bg,
            "patch_icpp(2)%vel(2)": v_drive,
            "patch_icpp(2)%pres": p_drive,
            "patch_icpp(2)%alpha_rho(1)": rho_drive,
            "patch_icpp(2)%alpha(1)": 1.0,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
