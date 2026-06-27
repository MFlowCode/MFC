#!/usr/bin/env python3
"""
1D Euler-Euler bubble, acoustically driven, at low Mach -- acoustic-substep validation.

Physics:
  - Single carrier fluid (model_eqns=2, num_fluids=1) seeded with a dilute population of
    Euler-Euler bubbles (bubbles_euler='T', adv_n='T', polytropic Rayleigh-Plesset).
  - Uniform bubbly background (void fraction vf0) on a periodic [0,1] domain with a small
    Gaussian pressure pulse centered at x=0.5. The pulse drives a gentle bubble response
    (compression/expansion) -- a MODERATE forcing, the regime the co-subcycle targets.
  - Low Mach: u0 = M*c with M=0.1 (c=sqrt(gamma)~1.18, normalised p0=rho0=1), so the split
    integrator takes the advective-CFL outer step and resolves acoustics + co-subcycles the
    bubble dynamics on the micro-steps.

Primary validation (R(t)/alpha(t) match):
  Run split (MFC_SPLIT=T -> acoustic_substepping='T') vs standard (MFC_SPLIT=F, the full
  RK3 + adaptive bubble Strang split) at the SAME resolution and physical t_stop, and
  compare the acoustically driven bubble response at the pulse-center probe. The
  co-subcycled equivalent radius R(t) tracks the standard solver in sign, magnitude, and
  time-trend; the void fraction alpha(t) under-responds because the bubble number density
  rides the slow advective transport (see docs/documentation/acousticSubstepping.md,
  "Euler-Euler bubbles"). Conservation (carrier mass + total energy) is exact and the void
  fraction stays bounded and positive.

  Diagnostics (read from <case>/D/ 1D ASCII): void fraction (prim slot eqn_idx%alf) and the
  bubble radius moment at the center cell over the save sequence.

Toggle baseline vs split via the MFC_SPLIT env var:
  MFC_SPLIT=F  -> acoustic_substepping='F' (standard RK3 + adaptive bubble integrator) [default]
  MFC_SPLIT=T  -> acoustic_substepping='T' (split-explicit, bubble dynamics co-subcycled)

NOTE on scope: this is MODERATE forcing (0.5% pulse). A violent single-bubble collapse
(strong forcing) drives the adaptive bubble sub-integrator to a stiffness convergence
failure regardless of the acoustic substep count -- documented out of scope (see
docs/documentation/acousticSubstepping.md, "Euler-Euler bubbles").
"""

import json
import math
import os

split = os.environ.get("MFC_SPLIT", "F").upper()
acoustic_substepping = "T" if split == "T" else "F"

gamma = 1.4
c0 = math.sqrt(gamma)  # carrier sound speed ~1.183 (normalised p0=rho0=1)
Mach = 0.1
u0 = Mach * c0  # uniform mean flow

L = 1.0
Nx = 99  # 100 cells, dx = 1/100
dx = L / (Nx + 1)

CFL = 0.5
dt_init = CFL * dx / c0  # acoustic-CFL initial guess (cfl_adap_dt refines it)

vf0 = 1.0e-3  # background void fraction (dilute)
pulse = 0.005  # 0.5% Gaussian pressure pulse -> MODERATE bubble forcing

print(
    json.dumps(
        {
            "run_time_info": "T",
            "m": Nx,
            "n": 0,
            "p": 0,
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "dt": dt_init,
            "cfl_adap_dt": "T",
            "cfl_target": CFL,
            "n_start": 0,
            "t_stop": 0.15,
            "t_save": 0.01,
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "alt_soundspeed": "F",
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "F",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "acoustic_substepping": acoustic_substepping,
            "n_acoustic_substeps": 20,
            "acoustic_div_damp": 0.4,
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "fd_order": 1,
            "probe_wrt": "T",
            "num_probes": 1,
            "probe(1)%x": 0.5,
            # Patch 1: uniform bubbly background + Gaussian pressure pulse
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%vel(1)": u0,
            "patch_icpp(1)%pres": f"1.0 + {pulse} * exp(-((x - 0.5) / 0.15)**2)",
            "patch_icpp(1)%alpha_rho(1)": (1.0 - vf0) * 1.0,
            "patch_icpp(1)%alpha(1)": vf0,
            "patch_icpp(1)%r0": 1.0,
            "patch_icpp(1)%v0": 0.0,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            # Euler-Euler bubbles (co-subcycled with the acoustic substep when split)
            "bubbles_euler": "T",
            "bubble_model": 2,
            "polytropic": "T",
            "adv_n": "T",
            "polydisperse": "F",
            "thermal": 3,
            "nb": 1,
            "bub_pp%R0ref": 1.0,
            "bub_pp%p0ref": 1.0,
            "bub_pp%rho0ref": 1.0,
            "bub_pp%ss": 0.0,
            "bub_pp%pv": 0.0,
            "bub_pp%mu_l": 0.1,
            "bub_pp%gam_g": 1.4,
        }
    )
)
