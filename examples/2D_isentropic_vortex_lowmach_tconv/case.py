#!/usr/bin/env python3
"""
Temporal-convergence study for the split-explicit low-Mach integrator
(acoustic_substepping). Smooth translating isentropic vortex at M_inf=0.1 on a
periodic [-3,3]^2 domain — the same physics as 2D_isentropic_vortex_lowmach, but
configured as a fixed-grid Delta t-refinement harness.

Method (successive-difference Richardson, fixed 32x32 grid):
  The spatial discretization is identical across runs, so it cancels in the
  difference of two solutions and ||u_{Δt} - u_{Δt/2}|| isolates the TEMPORAL
  error. cfl_const_dt gives uniform macro steps; the final step is clipped to
  land exactly on t_stop, so every run reaches the same final time.

Two regimes (select via env vars; both refine the advective/macro Δt by 2x):

  1. ACOUSTIC CFL HELD FIXED (the design measurement).
     Hold the acoustic micro-step Δτ fixed by scaling n_acoustic_substeps with
     Δt, so only the macro (Wicker–Skamarock RK3) coupling is refined:
        MFC_CFL=0.40 MFC_NSUB=96   MFC_DAMP=0.0
        MFC_CFL=0.20 MFC_NSUB=48   MFC_DAMP=0.0
        MFC_CFL=0.10 MFC_NSUB=24   MFC_DAMP=0.0
        MFC_CFL=0.05 MFC_NSUB=12   MFC_DAMP=0.0
     Measured L2 rate of the conserved vector: ~2.00, 1.98 -> SECOND ORDER.
     The split-explicit WS-RK3 macro scheme is 2nd-order in the advective step.

  2. AUTO SUBSTEPS (Δτ ∝ Δt, n_substeps≈(|u|+c)/|u| set by the Mach number).
        MFC_CFL=0.4,0.2,0.1,0.05  MFC_NSUB=0  MFC_DAMP=0.0
     Measured rate trends to ~1.0: the forward-backward (symplectic-Euler)
     acoustic micro-step is 1st-order in Δτ, so refining Δτ together with Δt
     limits the practical rate to first order. This is the expected behaviour of
     a split-explicit / forward-backward scheme and is documented in
     docs/documentation/acousticSubstepping.md.

  acoustic_div_damp must be 0 for the study: it is a per-micro-step fixed
  fraction (not Δτ-scaled), so it does not vanish under refinement and would
  otherwise mask the temporal error.

Compare cons.<var>.00.000001.dat (final dump, ASCII "x y value") across runs.
"""

import json
import math
import os

gamma = 1.4
pi_inf = 0.0
c0 = math.sqrt(gamma)
M_inf = 0.1
U0 = M_inf * c0
V0 = 0.0
epsilon = 0.5

Ldom = 6.0
Nx = 32
Ny = 32
dx = Ldom / Nx

CFL = float(os.environ.get("MFC_CFL", "0.4"))
NSUB = int(os.environ.get("MFC_NSUB", "96"))
DAMP = float(os.environ.get("MFC_DAMP", "0.0"))
T_stop = float(os.environ.get("MFC_TSTOP", "5.0"))
dt_init = CFL * dx / c0

print(
    json.dumps(
        {
            "run_time_info": "F",
            "x_domain%beg": -3.0,
            "x_domain%end": 3.0,
            "y_domain%beg": -3.0,
            "y_domain%end": 3.0,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": dt_init,
            "cfl_const_dt": "T",
            "cfl_target": CFL,
            "t_step_start": 0,
            "n_start": 0,
            "t_stop": T_stop,
            "t_save": T_stop,
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "alt_soundspeed": "F",
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": "rk3",
            "weno_order": 5,
            "weno_eps": 1.0e-16,
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
            "acoustic_substepping": "T",
            "n_acoustic_substeps": NSUB,
            "acoustic_div_damp": DAMP,
            "format": "binary",
            "precision": "double",
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "parallel_io": "F",
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": Ldom,
            "patch_icpp(1)%length_y": Ldom,
            "patch_icpp(1)%hcid": 283,
            "patch_icpp(1)%epsilon": epsilon,
            "patch_icpp(1)%vel(1)": U0,
            "patch_icpp(1)%vel(2)": V0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": pi_inf,
        }
    )
)
