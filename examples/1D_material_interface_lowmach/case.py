#!/usr/bin/env python3
"""
1D material-interface advection at low Mach — two-tier acoustic-substep validation
(split-explicit low-Mach integrator).

Physics:
  - Two ideal gases (model_eqns=2, num_fluids=2) on a periodic [0,1] domain.
  - A heavy-fluid slab (fluid 2, rho=3) in a light ambient (fluid 1, rho=1) forms a
    crisp material interface (sharp alpha / density contrast). The interface is RESOLVED
    over a few cells (tanh, half-width w ~ 4 dx): a centered transport scheme advects a
    resolved profile without Gibbs ringing, whereas a 1-cell-sharp contact would ring
    (the substep's convective mass/energy transport is centered, not upwinded — see the
    task report; the robust tier protects the acoustic pressure flux, not mass transport).
  - Pressure and velocity are UNIFORM (p0, u0), so the exact solution is the slab advecting
    at u0 with p kept flat; any pressure deviation at the interface is numerical ringing.
  - Low Mach: u0 = M*c with M=0.05 (c ~ O(1)), so the split integrator takes the
    advective-CFL outer step and resolves acoustics by the substep loop.

Two-tier check: the alpha-jump criterion flags the interface band, so the robust
WENO+HLLC-acoustic tier acts there while the uniform far field uses the cheap centered tier.

Diagnostics (read from <case>/D/ ASCII files, 1D):
  1. Non-oscillatory: pressure (prim slot E) stays ~p0 across the interface (no ringing).
  2. Per-species mass conservation: sum_x of each partial density (cons slots 1,2) is
     invariant in time to ~machine precision (periodic BCs, pure advection).

Toggle baseline (full HLLC) vs split via the MFC_SPLIT env var:
  MFC_SPLIT=F  -> acoustic_substepping='F' (standard RK3 + full HLLC)  [default]
  MFC_SPLIT=T  -> acoustic_substepping='T' (split-explicit two-tier)
Run both at the SAME resolution for the baseline-vs-split comparison.
"""

import json
import math
import os

split = os.environ.get("MFC_SPLIT", "F").upper()
acoustic_substepping = "T" if split == "T" else "F"

# Two distinct ideal gases (exercises the mixture gamma/pi_inf sums)
gamma1 = 1.4
gamma2 = 2.4
eps = 1.0e-6  # minority volume fraction (keeps mixture EOS well-defined)

p0 = 1.0
rho1 = 1.0  # light ambient fluid
rho2 = 3.0  # heavy slab fluid (material density contrast)

c1 = math.sqrt(gamma1 * p0 / rho1)  # light-fluid sound speed ~1.183
Mach = 0.05
u0 = Mach * c1  # uniform advection velocity

L = 1.0
Nx = 199  # 200 cells, dx = 1/200
dx = L / (Nx + 1)
w = 4.0 * dx  # interface half-width (resolved over a few cells)

CFL = 0.5
dt_init = CFL * dx / c1

# Advect the slab ~0.25 of the domain so the interface is well resolved; short run.
T_stop = 0.25 / u0
t_save = T_stop / 5.0

# Resolved slab indicator phi(x) ~ 1 inside [0.4,0.6], ~0 outside (tanh shoulders).
phi = f"(0.5*(tanh((x - 0.4)/{w}) - tanh((x - 0.6)/{w})))"
a2 = f"({eps} + {1.0 - 2.0 * eps}*{phi})"  # fluid-2 volume fraction
a1 = f"(1.0 - {a2})"

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
            "num_patches": 1,
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
            "acoustic_substepping": acoustic_substepping,
            "n_acoustic_substeps": 0,
            "acoustic_div_damp": 0.1,
            "format": "binary",
            "precision": "double",
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: resolved heavy slab in light ambient (analytic tanh profiles)
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%vel(1)": u0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": f"{rho1}*{a1}",
            "patch_icpp(1)%alpha_rho(2)": f"{rho2}*{a2}",
            "patch_icpp(1)%alpha(1)": a1,
            "patch_icpp(1)%alpha(2)": a2,
            "fluid_pp(1)%gamma": 1.0 / (gamma1 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(2)%gamma": 1.0 / (gamma2 - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
        }
    )
)
