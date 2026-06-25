#!/usr/bin/env python3
"""
1D uniform advection — slow-flux HLLC regression for acoustic_substepping
(split-explicit low-Mach time integrator), Task 3.

Physics / purpose:
  - Single-fluid stiffened-gas air (gamma=1.4, pi_inf=0) on a periodic domain.
  - Perfectly uniform initial state: constant density rho0, constant pressure p0,
    constant velocity u0 everywhere (no perturbation of any kind).
  - With acoustic_substepping = 'T', the HLLC flux returns ONLY the slow
    (advective) part: mass and total-energy face fluxes are zero (deferred to a
    later acoustic substep), and the momentum face flux is the pure contact-
    upwinded convective flux rho*u0**2 (no pressure term, no acoustic u+-c
    dissipation, no low-Mach pcorr).
  - For a uniform field every face sees identical left/right states, so the
    momentum flux is the same constant rho0*u0**2 at every face; its divergence
    is zero. Mass and energy fluxes are identically zero. Hence the uniform state
    must remain uniform to round-off after a few steps (pure, distortion-free
    translation). NOTE: at this task there is no acoustic substep yet, so only the
    slow momentum flux is exercised; full split-explicit dynamics arrive later.
"""

import json
import math

# Physical parameters
gamma = 1.4
p0 = 101325.0  # background pressure [Pa]
rho0 = 1.204  # background density [kg/m^3]
c0 = math.sqrt(gamma * p0 / rho0)  # ~343 m/s
Mach = 0.01
u0 = Mach * c0  # background velocity ~3.43 m/s (low-Mach)

# Domain
L = 1.0  # domain length [m]
Nx = 99  # cells (m parameter)
dx = L / (Nx + 1)

# CFL-based time stepping. Split-explicit (acoustic_substepping) mode requires
# CFL-based dt: the advective CFL sets dt and n_substeps is derived from the
# acoustic-to-advective wave-speed ratio. dt below is just the initial guess.
CFL = 0.5
dt = CFL * dx / c0
# Run a few advective steps (t_stop ~ a couple of advective dt).
t_stop = 3.0 * (CFL * dx / u0)

# Stiffened-gas EOS: Gamma = 1/(gamma-1), pi_inf = 0 for ideal gas
Gamma = 1.0 / (gamma - 1.0)
pi_inf = 0.0

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": dt,
            "cfl_adap_dt": "T",
            "cfl_target": CFL,
            "t_step_start": 0,
            "n_start": 0,
            "t_stop": t_stop,
            "t_save": t_stop / 2.0,
            # Simulation Algorithm Parameters
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
            # Acoustic substepping (split-explicit low-Mach): ON — exercises slow flux.
            "acoustic_substepping": "T",
            # Formatted Database Files Structure Parameters
            "format": "binary",
            "precision": "double",
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: perfectly uniform background (constant rho, p, u)
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%vel(1)": u0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluid physical parameters (stiffened-gas air, ideal limit)
            "fluid_pp(1)%gamma": Gamma,
            "fluid_pp(1)%pi_inf": pi_inf,
        }
    )
)
