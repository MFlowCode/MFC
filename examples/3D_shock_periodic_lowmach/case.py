#!/usr/bin/env python3
"""
3D periodic shock + dual transverse shear (single fluid, moderate Mach) — the 3D
TRANSVERSE-MOMENTUM validation for the split-explicit low-Mach two-tier acoustic
substep (acoustic_substepping).

Why this case exists
--------------------
The 2D shock-transverse gate (examples/2D_shock_transverse_lowmach) exercises only ONE
off-normal momentum slot per face. A 3D x-normal face has TWO off-normal components
(rho*u*v and rho*u*w), and the robust (full-HLLC) tier must map BOTH live HLLC
rotated-flux transverse slots (d != i) back into the stored global-component frozen
flux. A direction-indexing error in that 3D transverse mapping, in the per-direction
reconstruction, or in the 3-direction delta apply would only surface here.

This case forces both transverse slots: an x-normal Sod-like driver launches a
right-running shock through a field that carries transverse velocity components in BOTH
y (vel(2)) and z (vel(3)) which JUMP across the driver interface. At every flagged
x-face the transverse momenta rho*u*v and rho*u*w are then nonzero and varying, so the
d=2 AND d=3 transverse slots of the live HLLC flux are actively transported.

Stable-BC by construction: PERIODIC in all three directions (no outflow boundary, so the
separate centered-smooth-tier outflow instability cannot occur). A uniform background
mean flow u_bg in x bounds the advective outer step (~dx/|u|). The field is UNIFORM in y
and z, so the exact transverse behaviour is a clean diagnostic: with no y/z gradients and
no transverse pressure gradient, rho*v and rho*w are advected as contact quantities and
must stay UNIFORM in y and z. Any frame/sign/index error in the 3D d != i mapping would
break y/z-uniformity or disagree with the full-HLLC baseline.

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
# Transverse velocities JUMP across the driver interface (shear in BOTH y and z): this is
# what makes BOTH off-normal momentum slots nonzero and varying at the flagged x-faces.
v_amb = 0.15 * c_amb
v_drive = -0.15 * c_amb
w_amb = -0.10 * c_amb
w_drive = 0.10 * c_amb

L = 1.0
Nx = 47  # 48 cells in x
Ny = 31  # 32 cells in y (uniform-in-y,z field)
Nz = 31  # 32 cells in z
dx = L / (Nx + 1)
Ly = L * (Ny + 1) / (Nx + 1)  # cubic cells (dy = dz = dx)
Lz = L * (Nz + 1) / (Nx + 1)

CFL = 0.4
umag = math.sqrt(u_bg**2 + max(abs(v_amb), abs(v_drive)) ** 2 + max(abs(w_amb), abs(w_drive)) ** 2)
dt_init = CFL * dx / (umag + math.sqrt(gamma * p_drive / rho_drive))

T_stop = 0.04
t_save = T_stop / 4.0

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "y_domain%beg": 0.0,
            "y_domain%end": Ly,
            "z_domain%beg": 0.0,
            "z_domain%end": Lz,
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
            # Patch 1: ambient filling the domain (transverse velocities v_amb, w_amb)
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%y_centroid": 0.5 * Ly,
            "patch_icpp(1)%z_centroid": 0.5 * Lz,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": 2.0 * Ly,
            "patch_icpp(1)%length_z": 2.0 * Lz,
            "patch_icpp(1)%vel(1)": u_bg,
            "patch_icpp(1)%vel(2)": v_amb,
            "patch_icpp(1)%vel(3)": w_amb,
            "patch_icpp(1)%pres": p_amb,
            "patch_icpp(1)%alpha_rho(1)": rho_amb,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: high-pressure driver slab (0.15 < x < 0.35), transverse v_drive, w_drive
            "patch_icpp(2)%geometry": 9,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.25 * L,
            "patch_icpp(2)%y_centroid": 0.5 * Ly,
            "patch_icpp(2)%z_centroid": 0.5 * Lz,
            "patch_icpp(2)%length_x": 0.2 * L,
            "patch_icpp(2)%length_y": 2.0 * Ly,
            "patch_icpp(2)%length_z": 2.0 * Lz,
            "patch_icpp(2)%vel(1)": u_bg,
            "patch_icpp(2)%vel(2)": v_drive,
            "patch_icpp(2)%vel(3)": w_drive,
            "patch_icpp(2)%pres": p_drive,
            "patch_icpp(2)%alpha_rho(1)": rho_drive,
            "patch_icpp(2)%alpha(1)": 1.0,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
