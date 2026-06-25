#!/usr/bin/env python3
"""
2D two-fluid low-Mach acoustic pulse — correctness check for acoustic_substepping
with num_fluids > 1 (split-explicit low-Mach time integrator).

Physics:
  - Two ideal-gas fluids (model_eqns=2, num_fluids=2) on a periodic [-1,1]^2 domain.
  - A spatially UNIFORM two-fluid mixture (alpha_1 = alpha_2 = 0.5 everywhere),
    advecting at a uniform low-Mach mean flow, with a small Gaussian acoustic
    pressure pulse (hcid=264, dp/p = 1e-3) superposed.
  - The uniform mixture exercises the generalized mixture stiffened-gas EOS
    (gamma = sum_k alpha_k*gammas(k), with gamma_1 != gamma_2) and the new
    num_fluids>1 volume-fraction advance in the acoustic substep, without a sharp
    material interface — so per-fluid mass conservation is a clean diagnostic.

Diagnostics:
  1. Stability: no NaN in any field over a few hundred outer steps.
  2. Per-fluid mass conservation: sum_x alpha_k*rho_k constant (periodic BCs).
  3. Volume fraction stays uniform in [0,1] (alpha_k ~ 0.5, advected sanely).
"""

import json
import math

# Two ideal gases with distinct gamma (exercises the mixture gamma/pi_inf/qv sums)
gamma1 = 1.4
gamma2 = 2.4
piinf1 = 0.0
piinf2 = 0.0

p0 = 1.0
rho1 = 1.0
rho2 = 1.0
a1 = 0.5  # uniform volume fraction of fluid 1
a2 = 0.5  # uniform volume fraction of fluid 2

# Mixture EOS gamma_mix = a1*1/(g1-1) + a2*1/(g2-1); effective sound speed ~ O(1)
c0 = math.sqrt(gamma1 * (p0 + piinf1) / rho1)
M_inf = 0.1
U0 = M_inf * c0
V0 = 0.5 * U0

Ldom = 2.0
Nx = 64
Ny = 64
dx = Ldom / Nx

CFL = 0.5
dt_advect = CFL * dx / math.sqrt(U0**2 + V0**2)
dt_init = dt_advect

T_stop = 300.0 * dt_advect  # ~ a few hundred outer steps (each runs the acoustic substep loop)
t_save = T_stop / 5.0

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": -1.0,
            "x_domain%end": 1.0,
            "y_domain%beg": -1.0,
            "y_domain%end": 1.0,
            "m": Nx,
            "n": Ny,
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
            # Acoustic substepping (split-explicit low-Mach, multi-fluid)
            "acoustic_substepping": "T",
            "n_acoustic_substeps": 0,
            "acoustic_div_damp": 0.1,
            "format": "binary",
            "precision": "double",
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: uniform two-fluid mixture filling the whole domain,
            #   with a small Gaussian acoustic pressure pulse via hcid=264.
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": Ldom,
            "patch_icpp(1)%length_y": Ldom,
            "patch_icpp(1)%hcid": 264,
            "patch_icpp(1)%vel(1)": U0,
            "patch_icpp(1)%vel(2)": V0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": a1 * rho1,
            "patch_icpp(1)%alpha_rho(2)": a2 * rho2,
            "patch_icpp(1)%alpha(1)": a1,
            "patch_icpp(1)%alpha(2)": a2,
            "fluid_pp(1)%gamma": 1.0 / (gamma1 - 1.0),
            "fluid_pp(1)%pi_inf": gamma1 * piinf1 / (gamma1 - 1.0),
            "fluid_pp(2)%gamma": 1.0 / (gamma2 - 1.0),
            "fluid_pp(2)%pi_inf": gamma2 * piinf2 / (gamma2 - 1.0),
        }
    )
)
