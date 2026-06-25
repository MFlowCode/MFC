#!/usr/bin/env python3
"""
2D translating isentropic vortex at M ≈ 0.1 — low-Mach validation for
acoustic_substepping (split-explicit low-Mach time integrator).

Physics:
  - Single-fluid ideal gas (gamma=1.4, pi_inf=0) on a [-3,3]^2 periodic domain.
  - Isentropic vortex IC (hcid=283): Gauss-averaged cell values eliminating
    O(h^2) quadrature error in the nonlinear prim->cons conversion.
  - Mean flow at M_inf=0.1 (vel(1)=U0) translates the vortex.
  - Vortex strength eps=0.5 gives max perturbation velocity ~M_pert≈0.08,
    keeping the total Mach number below 0.2 everywhere.

Baseline mode (acoustic_substepping='F'):
  Outer dt is set by the acoustic CFL: dt ≈ CFL*dx/c0 ~ 0.01 (M-independent).

Split-explicit mode (acoustic_substepping='T'):
  Outer dt is set by the advective CFL: dt ≈ CFL*dx/U0 ~ 0.10 (~10x larger).
  Acoustic waves are resolved by ~10 forward-backward substeps per outer step.
  Expected wall-clock speedup ≈ M_inf^{-1} = 10 for this case.

Divergence damping (acoustic_div_damp=0.1) is active: the translating vortex
has nonzero divergence in the compressible regime, which the divergence-damping
term in the acoustic substep suppresses.

Usage:
  Baseline (acoustic-CFL dt, acoustic_substepping='F'):
    ./mfc.sh run examples/2D_isentropic_vortex_lowmach/case.py -n 1 -t pre_process simulation

  Split-explicit (acoustic_substepping='T'):
    Set ACOUSTIC_SUBSTEPPING = 'T' below and re-run.

Diagnostics probed:
  1. Stability: no NaN in any field.
  2. Mass conservation: |rho_total(t) - rho_total(0)| / rho_total(0) < 1e-12.
  3. Vortex accuracy: vortex center position vs. U0*t; kinetic energy drift.
  4. Divergence damping active: nonzero div(u) in baseline translating flow.
"""

import json
import math

# ── Physical parameters ───────────────────────────────────────────────────────
gamma = 1.4
pi_inf = 0.0

# Normalised sound speed: p0=1, rho0=1 → c0 = sqrt(gamma)
c0 = math.sqrt(gamma)  # ≈ 1.1832

# Mean flow: M_inf = 0.1
M_inf = 0.1
U0 = M_inf * c0  # ≈ 0.1183 (x-velocity)
V0 = 0.0  # no mean y-velocity

# Vortex strength (epsilon): max velocity perturbation ≈ eps*exp(0.5)/(2*pi*sqrt(2)) ≈ 0.186*eps
# eps=0.5 → M_pert ≈ 0.08, total M_max ≈ 0.18
epsilon = 0.5

# ── Grid ─────────────────────────────────────────────────────────────────────
# Domain: [-3, 3]^2 — vortex decays to exp(1-9)≈3e-4 at midpoint edges (r=3)
# so periodic BCs are satisfied to high accuracy.
Ldom = 6.0  # domain side length
Nx = 32
Ny = 32
dx = Ldom / Nx  # 0.1875 — vortex core (r<2) resolved by ~10 cells

# ── CFL-based time stepping ───────────────────────────────────────────────────
CFL = 0.5
dt_acoustic = CFL * dx / c0  # ≈ 0.0792 (baseline outer dt)
dt_advect = CFL * dx / U0  # ≈ 0.792  (split-mode outer dt, ~10x acoustic)

# Physical time: 0.5 (short, ~6 acoustic-CFL steps baseline or ~1 advective step split)
T_stop = 0.5
t_save = 0.1  # 5 save points

# Initial dt guess (adapted each step from cfl_target)
dt_init = dt_acoustic

# ── Acoustic substepping toggle ───────────────────────────────────────────────
# 'F' = baseline (acoustic-CFL outer step, standard SSP-RK3)
# 'T' = split-explicit (advective-CFL outer step, acoustic substep loop)
ACOUSTIC_SUBSTEPPING = "F"

# 0 = auto-compute n_substeps from outer dt / acoustic-CFL dt each step
N_ACOUSTIC_SUBSTEPS = 0

# Divergence damping coefficient (0.1 default); active for translating vortex
# where div(u) != 0 away from the vortex center.
ACOUSTIC_DIV_DAMP = 0.1

# ── Case dictionary ───────────────────────────────────────────────────────────
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain
            "x_domain%beg": -3.0,
            "x_domain%end": 3.0,
            "y_domain%beg": -3.0,
            "y_domain%end": 3.0,
            "m": Nx,
            "n": Ny,
            "p": 0,
            # CFL-based time stepping (required for acoustic_substepping)
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
            # Patch 1: Gauss-averaged isentropic vortex IC (hcid=283)
            #   vel(1) = U0 added to the vortex velocity by hcid=283 code
            #   vel(2) = V0 = 0 (no mean y-flow)
            #   epsilon = vortex strength (read by hcid=283 from patch_icpp%epsilon)
            #   Background state: rho0=1, p0=1 (overridden by hcid for r << 3)
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
            # Fluid: ideal gas (stiffened-gas limit with pi_inf=0)
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": pi_inf,
        }
    )
)
