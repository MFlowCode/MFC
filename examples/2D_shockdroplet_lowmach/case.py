#!/usr/bin/env python3
"""
2D periodic shock-droplet (two fluids, moderate Mach) — MULTIFLUID shock-interface
validation for the split-explicit low-Mach two-tier acoustic substep (acoustic_substepping).

Why this case exists
--------------------
The robust (full-HLLC) tier rebuilds the stiffened-gas mixture EOS at a flagged face from
the CELL volume fractions (exact only for num_fluids == 1). This case stresses that
approximation at a genuine num_fluids > 1 shock-interface and also provides a genuinely
2D transverse-momentum exercise: a planar shock launched in x diffracts around a circular
heavy-fluid droplet, producing nonzero transverse (y) momentum at the curved, flagged
shock + interface faces.

Stable-BC by construction: PERIODIC in both directions (no outflow boundary, avoiding the
separate centered-smooth-tier outflow instability). A uniform background mean flow u_bg in
x bounds the advective outer step. The droplet is sized/placed so the shock interacts with
it well inside the short test window.

Diagnostics: stability (no NaN, no dt collapse), non-oscillatory capture of the shock and
the droplet interface, per-species mass conservation, and L2/Linf of all conserved fields
vs the full-HLLC baseline at the same resolution.

Toggle via MFC_SPLIT env var:  MFC_SPLIT=T -> split, otherwise -> baseline full HLLC.
"""

import json
import math
import os

split = os.environ.get("MFC_SPLIT", "F").upper()
acoustic_substepping = "T" if split == "T" else "F"

gamma1 = 1.4  # ambient / driver (light gas)
gamma2 = 1.6  # droplet fluid (heavier gas)
eps = 1.0e-6

# Sod-like driver vs ambient (drives a moderate, right-running shock)
p_drive, rho_drive = 1.0, 1.0
p_amb, rho_amb = 0.1, 0.125
# Heavy droplet fluid at ambient pressure (density contrast ~16 vs ambient)
p_int, rho_int = 0.1, 2.0

c_amb = math.sqrt(gamma1 * p_amb / rho_amb)
u_bg = 0.3 * c_amb  # background mean flow in x (ambient Mach ~0.3)

Lx = 2.0
Ly = 1.0
Nx = 159  # 160 cells in x
Ny = 79  # 80 cells in y
dx = Lx / (Nx + 1)

R = 0.15  # droplet radius
x_drop = 1.2  # droplet centre (right half; shock launched from left reaches it)
y_drop = 0.5 * Ly

CFL = 0.4
dt_init = CFL * dx / (u_bg + math.sqrt(gamma1 * p_drive / rho_drive))

T_stop = 0.5
t_save = T_stop / 5.0

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": Lx,
            "y_domain%beg": 0.0,
            "y_domain%end": Ly,
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
            "num_patches": 3,
            "model_eqns": 2,
            "num_fluids": 2,
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
            # Patch 1: light ambient fluid 1 filling the domain
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * Lx,
            "patch_icpp(1)%y_centroid": 0.5 * Ly,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%length_y": Ly,
            "patch_icpp(1)%vel(1)": u_bg,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p_amb,
            "patch_icpp(1)%alpha_rho(1)": (1.0 - eps) * rho_amb,
            "patch_icpp(1)%alpha_rho(2)": eps * rho_int,
            "patch_icpp(1)%alpha(1)": 1.0 - eps,
            "patch_icpp(1)%alpha(2)": eps,
            # Patch 2: high-pressure driver band (x < 0.4), fluid 1
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.2,
            "patch_icpp(2)%y_centroid": 0.5 * Ly,
            "patch_icpp(2)%length_x": 0.4,
            "patch_icpp(2)%length_y": Ly,
            "patch_icpp(2)%vel(1)": u_bg,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": p_drive,
            "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_drive,
            "patch_icpp(2)%alpha_rho(2)": eps * rho_int,
            "patch_icpp(2)%alpha(1)": 1.0 - eps,
            "patch_icpp(2)%alpha(2)": eps,
            # Patch 3: heavy circular droplet of fluid 2
            "patch_icpp(3)%geometry": 2,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%x_centroid": x_drop,
            "patch_icpp(3)%y_centroid": y_drop,
            "patch_icpp(3)%radius": R,
            "patch_icpp(3)%vel(1)": u_bg,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%pres": p_int,
            "patch_icpp(3)%alpha_rho(1)": eps * rho_amb,
            "patch_icpp(3)%alpha_rho(2)": (1.0 - eps) * rho_int,
            "patch_icpp(3)%alpha(1)": eps,
            "patch_icpp(3)%alpha(2)": 1.0 - eps,
            "fluid_pp(1)%gamma": 1.0 / (gamma1 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(2)%gamma": 1.0 / (gamma2 - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
        }
    )
)
