#!/usr/bin/env python3
"""
3D spherical blast wave — Sod-like pressure jump on a sphere.
High-pressure sphere (p=1, rho=1) explodes into low-pressure gas (p=0.1, rho=0.125).
Produces a spherical shock front + contact discontinuity + rarefaction fan.

Grid: 63×63×63 (64 cells/dim).
MPI: 8 ranks (2×2×2) → 32 cells/rank/dim  ≥ 25 (WENO5 minimum).
"""
import json
import math

# Domain: unit cube centred at origin
L   = 1.0
N   = 63          # m=n=p → m+1=64

# Ideal gas γ=1.4
gamma = 1.4

# Inside sphere (r < R0)
rho_hi, p_hi = 1.0, 1.0
# Outside sphere
rho_lo, p_lo = 0.125, 0.1
R0 = 0.15         # sphere radius

# Shock speed estimate (Rankine-Hugoniot) for CFL
c_lo    = math.sqrt(gamma * p_lo / rho_lo)   # ~1.058
Ms      = math.sqrt((gamma + 1) / (2 * gamma) * (p_hi / p_lo)
                    + (gamma - 1) / (2 * gamma))  # ~2.0
v_shock = Ms * c_lo                               # ~2.12
c_hi    = math.sqrt(gamma * p_hi / rho_hi)        # ~1.183
max_wave = max(v_shock, c_hi)                      # ~2.12

dx = L / (N + 1)                       # ~0.01563
dt = 0.4 * dx / max_wave               # ~0.00295  (CFL=0.4)

Nt = 100   # simulate until t ≈ 0.295 s — shock reaches r ≈ 0.63
Ns = 10    # save every 10 steps → 11 frames

print(json.dumps({
    # Logistics
    "run_time_info": "T",

    # Domain
    "x_domain%beg": -L/2, "x_domain%end": L/2,
    "y_domain%beg": -L/2, "y_domain%end": L/2,
    "z_domain%beg": -L/2, "z_domain%end": L/2,
    "m": N, "n": N, "p": N,

    # Time
    "dt": dt, "t_step_start": 0, "t_step_stop": Nt, "t_step_save": Ns,

    # Numerics
    "model_eqns": 2, "num_fluids": 1,
    "alt_soundspeed": "F", "mpp_lim": "F", "mixture_err": "F",
    "time_stepper": 3, "weno_order": 5,
    "weno_eps": 1.0e-16, "teno": "T", "teno_CT": 1e-8,
    "null_weights": "F", "mp_weno": "F",
    "riemann_solver": 2, "wave_speeds": 1, "avg_state": 2,

    # Outflow BCs on all faces
    "bc_x%beg": -6, "bc_x%end": -6,
    "bc_y%beg": -6, "bc_y%end": -6,
    "bc_z%beg": -6, "bc_z%end": -6,

    # Output: binary format for the viz tool
    "format": 2, "precision": 2,
    "prim_vars_wrt": "T", "parallel_io": "T",

    # Patch 1 — low-pressure background (entire domain)
    "num_patches": 2,
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%length_y": L,
    "patch_icpp(1)%length_z": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p_lo,
    "patch_icpp(1)%alpha_rho(1)": rho_lo,
    "patch_icpp(1)%alpha(1)": 1.0,

    # Patch 2 — high-pressure sphere at origin
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": R0,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p_hi,
    "patch_icpp(2)%alpha_rho(1)": rho_hi,
    "patch_icpp(2)%alpha(1)": 1.0,

    # Ideal gas fluid properties
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}))
