#!/usr/bin/env python3
"""
3D smooth low-Mach acoustic pulse (single fluid, periodic) — the SMOOTH-FLOW companion to
examples/3D_shock_lowmach_gpu for the split-explicit low-Mach acoustic substep
(acoustic_substepping) on V100/OpenACC.

Why this case exists
--------------------
A localized strong shock raises the local Mach number and limits the advective-dt speedup
of the split scheme (it can be near break-even in shock-dominated cases). The big speedup
lives in the SMOOTH low-Mach bulk, where the robust (full-HLLC) shock-capture tier never
fires and the whole domain advances on the cheap centered/acoustic-substep path with an
advective-CFL outer step (~c/u ~ 1/M larger than the acoustic-CFL step of full HLLC).

Setup: a single uniform low-Mach background (rho0=1, p0=1, mean flow u_bg = 0.1*c0 in x) on
a periodic cube, with a tiny smooth Gaussian pressure pulse (dp = 1e-3 * p0, sigma = 0.05)
added by hardcoded IC hcid=304. The perturbation is smooth and infinitesimal, so the shock
sensor never trips: this is the large-speedup regime.

Toggle via MFC_SPLIT env var:  MFC_SPLIT=T -> split, otherwise -> baseline full HLLC.
"""

import json
import math
import os

split = os.environ.get("MFC_SPLIT", "F").upper()
acoustic_substepping = "T" if split == "T" else "F"

gamma = 1.4
p0, rho0 = 1.0, 1.0
c0 = math.sqrt(gamma * p0 / rho0)  # ~1.1832
u_bg = 0.1 * c0  # background Mach 0.1

# Domain: small cube so the fixed sigma=0.05 pulse (set inside hcid=304) is well resolved.
L = 0.5
Nx = 47  # 48^3 = 110592 cells; small for fast GPU turnaround
Ny = 47
Nz = 47
dx = L / (Nx + 1)

CFL = 0.4
# Acoustic-CFL outer dt (baseline / initial guess); split mode grows dt to the advective CFL.
dt_init = CFL * dx / (u_bg + c0)

T_stop = 0.6
t_save = T_stop / 4.0

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "y_domain%beg": 0.0,
            "y_domain%end": L,
            "z_domain%beg": 0.0,
            "z_domain%end": L,
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "dt": dt_init,
            "cfl_adap_dt": "T",
            "cfl_target": CFL,
            "t_step_start": 0,
            "n_start": 0,
            "t_stop": T_stop,
            "t_save": t_save,
            "num_patches": 1,
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
            "bc_z%beg": -1,
            "bc_z%end": -1,
            "acoustic_substepping": acoustic_substepping,
            "n_acoustic_substeps": 0,
            "acoustic_div_damp": 0.1,
            "format": "binary",
            "precision": "double",
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: uniform low-Mach background filling the periodic cube.
            # hcid=304 overlays a smooth Gaussian pressure pulse (dp=1e-3*pres, sigma=0.05)
            # centered at the patch centroid.
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%y_centroid": 0.5 * L,
            "patch_icpp(1)%z_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": L,
            "patch_icpp(1)%length_z": L,
            "patch_icpp(1)%vel(1)": u_bg,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%hcid": 304,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
