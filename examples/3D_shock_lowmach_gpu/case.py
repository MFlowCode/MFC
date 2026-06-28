#!/usr/bin/env python3
"""
3D low-Mach background + localized shock driver (single fluid) — the V100/OpenACC GPU
demonstrator for the split-explicit low-Mach two-tier acoustic substep
(acoustic_substepping).

This is a GPU-turnaround variant of examples/3D_shock_periodic_lowmach: same periodic,
single-fluid, transverse-shear, Sod-like-driver structure, but the BACKGROUND is made
genuinely low-Mach (u_bg ~ 0.1*c_amb, total background Mach ~0.12) so the smooth bulk sits
in the low-Mach regime where the split scheme's advective outer step is large. A localized
high-pressure slab still launches a right-running shock, so the robust (full-HLLC) tier
flags at the shock faces and the two-tier shock capture is exercised in 3D.

Demonstrates (compare MFC_SPLIT=T vs =F at identical resolution):
  - Stability: runs to completion with acoustic_substepping='T', no NaN, no dt collapse.
  - Correctness: split solution matches full-HLLC (L2/Linf of conserved fields); shock
    captured non-oscillatorily (no spurious overshoot).
  - Conservation: total mass, energy, and all three momenta to ~machine precision.
  - GPU performance: per-step "Time Avg" and total time split vs full-HLLC on the V100.

Toggle via MFC_SPLIT env var:  MFC_SPLIT=T -> split, otherwise -> baseline full HLLC.
"""

import json
import math
import os

split = os.environ.get("MFC_SPLIT", "F").upper()
acoustic_substepping = "T" if split == "T" else "F"

gamma = 1.4

# Moderate (2:1) pressure driver vs ambient. A full 10:1 Sod jump pushes the local Mach to
# order 1 (no longer "low-Mach") and is smeared markedly more by the substepped scheme than
# by tiny-dt full HLLC; a moderate jump keeps the local Mach modest, still launches a genuine
# shock that flags the robust (full-HLLC) tier, and keeps split and full in close agreement.
p_drive, rho_drive = 1.0, 1.0
p_amb, rho_amb = 0.5, 0.5

c_amb = math.sqrt(gamma * p_amb / rho_amb)
# Genuinely low-Mach background: mean flow Mach ~0.1, modest transverse shear ~0.05.
u_bg = 0.1 * c_amb
# Transverse velocities JUMP across the driver interface (shear in BOTH y and z): keeps both
# off-normal momentum slots nonzero and varying at the flagged x-faces, while staying low Mach.
v_amb = 0.05 * c_amb
v_drive = -0.05 * c_amb
w_amb = -0.05 * c_amb
w_drive = 0.05 * c_amb

L = 1.0
Nx = 47  # 48 cells in x
Ny = 47  # 48 cells in y
Nz = 47  # 48 cells in z (48^3 = 110592 cells; small for fast GPU turnaround)
dx = L / (Nx + 1)
Ly = L * (Ny + 1) / (Nx + 1)  # cubic cells (dy = dz = dx)
Lz = L * (Nz + 1) / (Nx + 1)

# Outer-step CFL (env-overridable). With a strong embedded shock the split scheme's
# advective-CFL outer step must stay modest for stability, so the demonstrator uses a
# conservative CFL shared by both split and full-HLLC for an apples-to-apples comparison.
CFL = float(os.environ.get("MFC_CFL", "0.4"))
umag = math.sqrt(u_bg**2 + max(abs(v_amb), abs(v_drive)) ** 2 + max(abs(w_amb), abs(w_drive)) ** 2)
dt_init = CFL * dx / (umag + math.sqrt(gamma * p_drive / rho_drive))

# Keep the run short enough that the shock stays localized (does not traverse the periodic
# domain), so the split-vs-full pointwise comparison is clean (no wrap-around mixing).
T_stop = 0.3
t_save = T_stop / 6.0

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
