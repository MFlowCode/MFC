#!/usr/bin/env python3
"""
2D Gresho vortex at M ≈ 0.01 — low-Mach validation for acoustic_substepping.

Physics:
  - Single-fluid ideal gas (gamma=1.4, pi_inf=0) on a [0,1]^2 periodic domain.
  - Gresho vortex IC (hcid=203): rotating vortex with pressure balanced by
    centrifugal acceleration.  The exact incompressible solution is steady, so
    any temporal drift measures numerical dissipation/dispersion.
  - Mach ≈ 0.001 set via vel(2)=Ma0 (the hcid=203 Mach-number parameter).

Usage:
  Baseline (acoustic-CFL dt, acoustic_substepping='F'):
    ./mfc.sh run examples/2D_gresho_lowmach/case.py -n 1 -t pre_process simulation

  Split-explicit (acoustic_substepping='T'):
    Edit ACOUSTIC_SUBSTEPPING to 'T' and re-run.

Concerns probed:
  1. Multi-D substep count: auto n_substeps uses sum(|vel|) in the denominator
     but the outer dt is set by the per-dimension max wave speed. In a 2D
     rotating vortex (|u| ~ |v|) this may under-count substeps => acoustic-CFL
     violation. Test with N_ACOUSTIC_SUBSTEPS=0 (auto), then 2x and 4x manual.
  2. Odd-even / checkerboard: acoustic_div_damp (default 0.1) is the mitigation.
     Increase to 0.2 if checkerboard arises in the split-mode pressure field.
"""

import json
import math

# ── Physical parameters ───────────────────────────────────────────────────────
gamma = 1.4
pi_inf = 0.0

# hcid=203 uses vel(2) as the vortex Mach number Ma0:
#   umax = 2*pi*rmax * (vel(2) * c0), where c0=sqrt(gamma) in normalised units.
# Setting Ma0=1e-3 gives M ~ 0.001 (low Mach, acoustic substepping relevant).
Ma0 = 1.0e-3  # vortex Mach number

# Normalised sound speed: with p0=1, rho0=1 the EOS gives c0=sqrt(gamma*p/rho)
c0 = math.sqrt(gamma)  # ≈ 1.183 in normalised units

# ── Grid ─────────────────────────────────────────────────────────────────────
# 32x32 keeps baseline feasible; increase to 64+ for production accuracy.
Nx = 32
Ny = 32
dx = 1.0 / Nx

# ── CFL-based time stepping ───────────────────────────────────────────────────
# BASELINE mode: outer dt ~ CFL*dx/(c0+umax) ≈ CFL*dx/c0 (umax << c0).
# SPLIT mode: outer dt ~ CFL*dx/umax = CFL*dx/(Ma0*c0), ~1/Ma0 = 1000x larger.
# n_substeps (auto) ~ 1/Ma0 = 1000 per outer step in split mode.
CFL = 0.5
dt_acoustic = CFL * dx / c0  # ~0.013 (acoustic CFL, baseline outer dt)
dt_advect = CFL * dx / (Ma0 * c0)  # ~13.2 (advective CFL, split-mode outer dt)

# Physical time: run for 5 split-mode outer steps (= 5 * dt_advect ≈ 66).
# Baseline will do ~5000 acoustic-CFL steps over the same physical time.
T_stop = 5.0 * dt_advect  # ~66 normalised time units
t_save = T_stop / 5.0  # 5 save points

# Initial dt guess (recomputed each step by the solver from cfl_target)
dt_init = dt_acoustic

# ── Acoustic substepping toggle ───────────────────────────────────────────────
# ACOUSTIC_SUBSTEPPING = 'F'  -> baseline (acoustic-CFL outer step, standard RK3)
# ACOUSTIC_SUBSTEPPING = 'T'  -> split-explicit (advective-CFL outer step, substepped)
ACOUSTIC_SUBSTEPPING = "F"  # change to 'T' for split mode (requires cfl_adap_dt)

# N_ACOUSTIC_SUBSTEPS: 0 = auto (from local u/c ratio); >0 = manual override.
# Concern #1: in 2D rotating vortex auto may under-count. Try 0, then 2x, 4x.
N_ACOUSTIC_SUBSTEPS = 0  # 0 = auto (adapt to outer dt each step)

# Divergence damping coefficient (default 0.1).
# Concern #2: increase if checkerboard noise appears in pressure/velocity.
ACOUSTIC_DIV_DAMP = 0.1

# ── Case dictionary ───────────────────────────────────────────────────────────
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "m": Nx,
            "n": Ny,
            "p": 0,
            # Time stepping (CFL-based — required for acoustic_substepping)
            "dt": dt_init,
            "cfl_adap_dt": "T",
            "cfl_target": CFL,
            "t_step_start": 0,
            "n_start": 0,
            "t_stop": T_stop,
            "t_save": t_save,
            # Simulation algorithm
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
            # Periodic boundary conditions
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            # Acoustic substepping (split-explicit low-Mach)
            "acoustic_substepping": ACOUSTIC_SUBSTEPPING,
            "n_acoustic_substeps": N_ACOUSTIC_SUBSTEPS,
            "acoustic_div_damp": ACOUSTIC_DIV_DAMP,
            # Output
            "format": "binary",
            "precision": "double",
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: Gresho vortex with density correction (hcid=203)
            # vel(2)=Ma0 is the Mach-number parameter used by the IC formula.
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%hcid": 203,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": Ma0,
            "patch_icpp(1)%pres": 1.0,
            # Fluid physical parameters (stiffened-gas, ideal limit)
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": pi_inf,
        }
    )
)
