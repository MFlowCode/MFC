#!/usr/bin/env python3
"""
2D power-law (Herschel-Bulkley with tau0=0) shear-THICKENING Poiseuille channel.

Companion to examples/2D_poiseuille_nn (shear-thinning, nn = 0.7); this case sets
nn = 1.5 > 1 so the effective viscosity mu = K*gamma_dot^(n-1) INCREASES with shear
rate. The signature of a correct shear-thickening power-law term is a velocity
profile sharper (more pointed) at the centerline than a parabola.

A constant body acceleration g_x drives a fully-developed channel flow between two
stationary no-slip walls; the steady velocity profile is compared against the
closed-form power-law analytic solution:

    u(y) = (n/(n+1)) * (rho*g/K)^(1/n) * ( H^((n+1)/n) - |y-H|^((n+1)/n) )

with channel height L_y, half-height H = L_y/2, centerline at y = H, walls at
y = 0 and y = L_y. For n > 1 the effective viscosity mu = K*gamma_dot^(n-1) -> 0 at
the shear-free centerline (rather than diverging as for n < 1), so the clamp to
[mu_min, mu_max] is never active and the analytic power-law profile is an exact
reference everywhere.

Parameters (nondimensional MFC units):
    K   = 5e-2   (consistency index)
    nn  = 1.5    (flow index; shear-thickening, profile more pointed than parabolic)
    g_x = 5e-2   (driving body acceleration)
    rho = 1.0, pres = 10  -> sound speed ~3.74, u_max ~0.013 => Mach ~3e-3
    grid: m = 24 (x, periodic), n = 63 (y), L_x = L_y = 0.2, H = 0.1
          (m >= 24 so a 2-rank y-split satisfies the WENO5 decomposition minimum;
           n = 63 keeps the explicit viscous CFL dt ~ dy^2 rho/mu_max tractable)
    For n > 1 the maximum physical viscosity is at the WALL (highest shear rate):
          mu_wall = K^(1/n) * (rho*g*H)^((n-1)/n) ~ 0.023.
    mu_max = 0.035 ~ 1.5*mu_wall: just above the physical max so the clamp stays
    inactive (analytic profile exact) while keeping the viscous timestep large --
    the explicit dt scales as 1/mu_max, so the previous mu_max = 10 was the source
    of the catastrophically small dt.
    cfl_adap_dt with cfl_target = 0.3, integrated to t_stop = 0.9 (~2 wall-viscous
    diffusion times tau = H^2 rho/mu_wall ~ 0.43)

Compare with compare_analytic.py.
"""

import json

# Channel / fluid parameters
L_x = 0.2
L_y = 0.2
rho = 1.0
pres = 10.0
K = 5.0e-2
nn = 1.5
g_x = 5.0e-2

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
            "n": 63,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.3,
            "n_start": 0,
            "t_save": 0.09,
            "t_stop": 0.9,
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
            "fluid_pp(1)%mu_max": 0.035,
        }
    )
)
