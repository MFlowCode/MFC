#!/usr/bin/env python3
"""
1D weakly-compressible acoustic pulse (M ≈ 0.01) — baseline harness for
acoustic_substepping (split-explicit low-Mach time integrator).

Physics:
  - Single-fluid stiffened-gas air (gamma=1.4, pi_inf=0) on a periodic domain.
  - Background state: p0=101325 Pa, rho0=1.204 kg/m^3.
  - Sound speed: c0 = sqrt(gamma*p0/rho0) ≈ 343 m/s.
  - Background velocity: u0 = M*c0 with M=0.01 ≈ 3.43 m/s (low-Mach).
  - Initial condition: small Gaussian pressure perturbation dp = 1e-3 * p0.

Set acoustic_substepping = 'T' (and build the numerics in later tasks) to
exercise the split-explicit pathway.  Default here is 'F' (baseline RK3).
"""

import json
import math

# Physical parameters
gamma = 1.4
p0 = 101325.0  # background pressure [Pa]
rho0 = 1.204  # background density [kg/m^3]
c0 = math.sqrt(gamma * p0 / rho0)  # ~343 m/s
Mach = 0.01
u0 = Mach * c0  # background velocity ~3.43 m/s

# Domain
L = 1.0  # domain length [m]
Nx = 199  # cells (m parameter)
dx = L / (Nx + 1)

# Time step: CFL ~ 0.5 based on acoustic wave speed
CFL = 0.5
dt = CFL * dx / c0
t_end = L / c0  # one acoustic transit
Nt = int(math.ceil(t_end / dt))
t_save = max(1, Nt // 10)

# Gaussian pressure pulse: dp = eps_p * p0 * exp(-((x-0.5)/sigma)^2)
# sigma = 0.05 * L (narrow pulse, well-resolved on the grid)
# Written as an analytic IC expression evaluated at cell centers.
eps_p = 1.0e-3
sigma = 0.05 * L

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
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": t_save,
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
            # Acoustic substepping (split-explicit low-Mach): off by default (baseline)
            "acoustic_substepping": "F",
            # Formatted Database Files Structure Parameters
            "format": "binary",
            "precision": "double",
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: uniform background + small Gaussian pressure pulse
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%vel(1)": u0,
            "patch_icpp(1)%pres": f"{p0} + {eps_p * p0} * exp(-((x - {0.5 * L}) / {sigma})**2)",
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluid physical parameters (stiffened-gas air, ideal limit)
            "fluid_pp(1)%gamma": Gamma,
            "fluid_pp(1)%pi_inf": pi_inf,
        }
    )
)
