#!/usr/bin/env python3
"""
1D shock-interface interaction (moderate Mach) — two-tier acoustic-substep robust-tier
validation (split-explicit low-Mach integrator).

Physics:
  - Two ideal gases (model_eqns=2, num_fluids=2) on [0,1].
  - A Sod-type pressure discontinuity at x=0.4 launches a right-running shock
    (shock Mach ~1.6, moderate) into a light ambient (fluid 1). The shock then impinges
    on a heavy material interface (fluid 2, rho=1.0 vs ambient 0.125) at x=0.65, producing
    reflected + transmitted waves through the interface.
  - The whole field carries a moderate background mean flow u_bg (ambient Mach ~0.3). This
    is REQUIRED for split mode: the outer step is the advective CFL ~dx/|u|, which is
    unbounded for a quiescent (u=0) shock tube; a nonzero mean flow bounds it. It also makes
    this a genuinely discontinuity-rich, moderate-Mach blended-cost case.
  - This is the discontinuity-rich stress case for the two-tier substep.

Two-tier robust flux: both the BASELINE full-HLLC solver (acoustic_substepping='F', the
default here) and SPLIT mode (acoustic_substepping='T') capture this shock-interface
interaction. The robust tier makes the mass, energy, AND convective-momentum transport full
HLLC at flagged faces; the convective momentum is made live through a live-minus-frozen delta
(the slow path's stored stage-entry convective flux is subtracted and the live HLLC flux
added), which keeps the slow/frozen cancellation exact on every RK stage and stays stable
when momentum advection is leading order — i.e. at a propagating shock. See
docs/documentation/acousticSubstepping.md.

Note on boundaries: the centered SMOOTH tier is unstable at extrapolation/outflow boundaries
carrying a mean flow (independent of shock capture; see the docs). This case uses
non-reflecting BCs (bc_x=-3), so for a clean split-mode run keep the structure interior over
the test window, or prefer a stable-BC case such as examples/1D_shock_periodic_lowmach.

Diagnostics:
  - Read from <case>/D/ ASCII (1D): prim.4 = pressure, prim.3 = velocity,
    prim.1/prim.2 = partial densities, prim.5/prim.6 = volume fractions.

Toggle baseline (full HLLC) vs split via the MFC_SPLIT env var:
  MFC_SPLIT=F  -> acoustic_substepping='F' (standard RK3 + full HLLC)  [default]
  MFC_SPLIT=T  -> acoustic_substepping='T' (split-explicit two-tier)
"""

import json
import math
import os

split = os.environ.get("MFC_SPLIT", "F").upper()
acoustic_substepping = "T" if split == "T" else "F"

gamma1 = 1.4  # ambient / driver (air-like)
gamma2 = 1.6  # interface fluid
eps = 1.0e-6

# Sod-like driver vs ambient (drives a moderate, right-running shock)
p_drive, rho_drive = 1.0, 1.0
p_amb, rho_amb = 0.1, 0.125
# Heavy material interface fluid at the ambient pressure
p_int, rho_int = 0.1, 1.0

c_amb = math.sqrt(gamma1 * p_amb / rho_amb)  # ambient sound speed ~1.06
u_bg = 0.3 * c_amb  # moderate background mean flow (ambient Mach ~0.3)

L = 1.0
Nx = 399  # 400 cells
dx = L / (Nx + 1)

CFL = 0.4
dt_init = CFL * dx / (u_bg + math.sqrt(gamma1 * p_drive / rho_drive))

# Run until the shock has crossed the interface and the reflected/transmitted
# structure is established (lab-frame shock reaches x=0.65 around t~0.12).
T_stop = 0.15
t_save = T_stop / 5.0

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
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "acoustic_substepping": acoustic_substepping,
            "n_acoustic_substeps": 0,
            "acoustic_div_damp": 0.1,
            "format": "binary",
            "precision": "double",
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: light ambient fluid 1 filling the domain
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%vel(1)": u_bg,
            "patch_icpp(1)%pres": p_amb,
            "patch_icpp(1)%alpha_rho(1)": (1.0 - eps) * rho_amb,
            "patch_icpp(1)%alpha_rho(2)": eps * rho_int,
            "patch_icpp(1)%alpha(1)": 1.0 - eps,
            "patch_icpp(1)%alpha(2)": eps,
            # Patch 2: high-pressure driver (x < 0.4), fluid 1
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.2 * L,
            "patch_icpp(2)%length_x": 0.4 * L,
            "patch_icpp(2)%vel(1)": u_bg,
            "patch_icpp(2)%pres": p_drive,
            "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_drive,
            "patch_icpp(2)%alpha_rho(2)": eps * rho_int,
            "patch_icpp(2)%alpha(1)": 1.0 - eps,
            "patch_icpp(2)%alpha(2)": eps,
            # Patch 3: heavy material interface fluid 2 (x > 0.65)
            "patch_icpp(3)%geometry": 1,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%x_centroid": 0.825 * L,
            "patch_icpp(3)%length_x": 0.35 * L,
            "patch_icpp(3)%vel(1)": u_bg,
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
