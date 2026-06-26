#!/usr/bin/env python3
"""
1D periodic shock tube (single fluid, moderate Mach) — interior shock-stability gate
for the split-explicit low-Mach two-tier acoustic substep (acoustic_substepping).

Stable-BC by construction: PERIODIC boundaries (bc_x = -1), so there is no outflow
boundary and the separate centered-smooth-tier outflow instability cannot occur. A
uniform background mean flow u_bg bounds the advective CFL (the outer step is ~dx/|u|,
unbounded for a quiescent tube). A Sod-like high-pressure region launches a right-running
shock into the interior; over the short test window the shock stays well inside the domain.

This isolates the convective-momentum fix: split mode (acoustic_substepping='T') must no
longer NaN / collapse dt in the INTERIOR, and the shock must be captured. Compare against
the full-HLLC baseline (acoustic_substepping='F') at the same resolution.

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
u_bg = 0.4 * c_amb  # background mean flow, ambient Mach ~0.4 (bounds the advective CFL)

L = 1.0
Nx = 199  # 200 cells, tiny
dx = L / (Nx + 1)

CFL = 0.4
dt_init = CFL * dx / (u_bg + math.sqrt(gamma * p_drive / rho_drive))

# Short window: shock has formed and propagated a few tens of cells but stays interior.
T_stop = 0.06
t_save = T_stop / 6.0

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "m": Nx,
            "n": 0,
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
            "acoustic_substepping": acoustic_substepping,
            "n_acoustic_substeps": 0,
            "acoustic_div_damp": 0.1,
            "format": "binary",
            "precision": "double",
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: ambient filling the domain
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%vel(1)": u_bg,
            "patch_icpp(1)%pres": p_amb,
            "patch_icpp(1)%alpha_rho(1)": rho_amb,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: high-pressure driver (0.15 < x < 0.35)
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.25 * L,
            "patch_icpp(2)%length_x": 0.2 * L,
            "patch_icpp(2)%vel(1)": u_bg,
            "patch_icpp(2)%pres": p_drive,
            "patch_icpp(2)%alpha_rho(1)": rho_drive,
            "patch_icpp(2)%alpha(1)": 1.0,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
