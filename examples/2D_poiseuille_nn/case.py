#!/usr/bin/env python3
"""
2D power-law (Herschel-Bulkley with tau0=0) Poiseuille channel.

Validation case for the non-Newtonian (shear-dependent) viscosity in the HLLC
viscous flux. A constant body acceleration g_x drives a fully-developed channel
flow between two stationary no-slip walls; the steady velocity profile is
compared against the closed-form power-law analytic solution:

    u(y) = (n/(n+1)) * (rho*g/K)^(1/n) * ( H^((n+1)/n) - |y-H|^((n+1)/n) )

with channel height L_y, half-height H = L_y/2, centerline at y = H, walls at
y = 0 and y = L_y. The effective dynamic viscosity in MFC is mu = K*gamma_dot^(n-1)
(tau0 = 0), clamped to [mu_min, mu_max]; for the committed nn = 0.7 the clamp is
never active (mu stays well below mu_max away from the centerline).

Parameters (nondimensional MFC units):
    K   = 2e-2   (consistency index)
    nn  = 0.7    (flow index; shear-thinning, profile blunter than parabolic)
    g_x = 8e-2   (driving body acceleration)
    rho = 1.0, pres = 10  -> sound speed ~3.74, u_max ~0.011 => Mach ~3e-3
    grid: m = 24 (x, periodic), n = 127 (y), L_x = L_y = 0.2, H = 0.1
          (m >= 24 so a 2-rank y-split satisfies the WENO5 decomposition minimum;
           the compact channel keeps the diffusive time small, ~0.12, so a short
           acoustic-CFL-limited run reaches steady state in ~1e4 steps)
    cfl_adap_dt with cfl_target = 0.3, integrated to t_stop = 1.2 (~10 diffusive times)

Result (compare_analytic.py, 2-rank run, last of 11 saves; flow is steady, u_max
converged monotonically with <0.1% change over the final steps):
    relative L2 error vs. analytic power-law profile : 6.8e-3 (0.68%)
    u_max numeric / analytic                         : 1.1149e-2 / 1.1121e-2
    mean/peak bluntness (parabola 0.667; n=0.7: 0.708): 0.706
    local balance K|du/dy|^n vs rho*g*(H-y), ratio    : 1.000 across the channel
The pointwise momentum balance holding to ~0.1% confirms the shear-dependent
viscosity in the HLLC viscous flux is correct.

A shear-thinning index nn = 0.7 is committed here (rather than a more singular
nn = 0.5): for n < 1 the effective viscosity mu = K*gamma_dot^(n-1) diverges at the
shear-free centerline, where any regularized solver must cap it, so the analytic
power-law profile is only an exact reference for n not too far below 1. At nn = 0.5
the near-wall momentum balance is still exact but u_max overshoots the (singular)
analytic by ~25% purely from the centerline cap; nn = 0.7 keeps a clearly
shear-thinning, blunter-than-parabolic profile while matching the analytic to <1%.
A Newtonian control (non_newtonian = F, constant mu) reproduces the parabola to
better than 0.01%.
"""

import json

# Channel / fluid parameters
L_x = 0.2
L_y = 0.2
rho = 1.0
pres = 10.0
K = 2.0e-2
nn = 0.7
g_x = 8.0e-2

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": L_x,
            "y_domain%beg": 0.0,
            "y_domain%end": L_y,
            "m": 24,
            "n": 127,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.3,
            "n_start": 0,
            "t_save": 0.12,
            "t_stop": 1.2,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1e-16,
            "mapped_weno": "T",
            "weno_Re_flux": "T",
            "mp_weno": "T",
            "weno_avg": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            # x: periodic channel; y: stationary no-slip walls
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -16,
            "bc_y%end": -16,
            "viscous": "T",
            # Constant body acceleration in +x: accel = g_x + k_x*sin(w_x*t - p_x)
            "bf_x": "T",
            "g_x": g_x,
            "k_x": 0.0,
            "w_x": 0.0,
            "p_x": 0.0,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "fd_order": 4,
            "parallel_io": "T",
            # Patch 1: full domain, quiescent
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * L_x,
            "patch_icpp(1)%y_centroid": 0.5 * L_y,
            "patch_icpp(1)%length_x": L_x,
            "patch_icpp(1)%length_y": L_y,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": pres,
            "patch_icpp(1)%alpha_rho(1)": rho,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluids Physical Parameters: single power-law (HB with tau0 = 0) fluid
            "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": 1.0 / K,
            "fluid_pp(1)%non_newtonian": "T",
            "fluid_pp(1)%K": K,
            "fluid_pp(1)%nn": nn,
            "fluid_pp(1)%tau0": 0.0,
            "fluid_pp(1)%hb_m": 1000.0,
            "fluid_pp(1)%mu_min": 1e-6,
            "fluid_pp(1)%mu_max": 10.0,
        }
    )
)
