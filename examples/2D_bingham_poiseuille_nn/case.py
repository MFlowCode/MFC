#!/usr/bin/env python3
"""
2D Bingham (Herschel-Bulkley with flow index n = 1, yield stress tau0 > 0)
Poiseuille channel.

Validation case for the YIELD-STRESS term of the non-Newtonian (Herschel-Bulkley)
viscosity. The companion examples/2D_poiseuille_nn validates the power-law term
(tau0 = 0); this case isolates the yield term by fixing n = 1 (so K = mu is a plain
Newtonian consistency) and turning on tau0 > 0. The signature of a correct yield
term is a rigid PLUG of uniform velocity near the centerline, where the shear
stress |tau| = rho*g*(H - y) drops below tau0 and the fluid stops yielding.

A constant body acceleration g_x drives a fully-developed channel flow between two
stationary no-slip walls (y = 0, L_y; half-height H = L_y/2, centerline y = H).
The steady Bingham profile (n = 1, K = mu, tau_w = rho*g*H):

  plug half-width from centerline   : y0   = tau0/(rho*g)
  sheared region (0 <= y <= H - y0) : u(y) = (1/(2*mu*rho*g)) *
                                              [ (tau_w - tau0)^2 - (rho*g*(H-y) - tau0)^2 ]
  plug (H - y0 <= y <= H)           : u_plug = (1/(2*mu*rho*g)) * (tau_w - tau0)^2
  upper half mirrors about y = H.

Requires tau_w = rho*g*H > tau0 for any flow.

Parameters (nondimensional MFC units):
    rho  = 1.0, H = 0.1 (L_y = 0.2)
    g_x  = 0.1          -> tau_w = rho*g*H = 1.0e-2
    tau0 = 4.0e-3       -> y0 = tau0/(rho*g) = 4.0e-2 = 0.4 H  (clear plug + shear)
    K    = mu = 5.0e-2  (n = 1 Bingham consistency; short diffusive time t_d = H^2 rho/mu = 0.2)
    hb_m = 1.0e4        (sharp Papanastasiou yield; uncapped plug mu = tau0*hb_m = 40)
    mu_min = 1e-6, mu_max = 1.0  (plug viscosity = 20*K, effectively rigid plug;
          kept modest because the explicit viscous timestep scales as 1/mu_max)
    pres = 10 -> sound speed ~3.74; u_plug ~ 3.6e-3 => Mach ~1e-3 (low Mach)
    grid: m = 24 (x, periodic; >= 24 for a 2-rank y-split WENO5 decomposition),
          n = 63 (y resolution of the plug; the viscous CFL dt ~ dy^2 rho/mu_max
          makes a finer grid prohibitively slow for an explicit solver), L_x = L_y = 0.2
    cfl_adap_dt, cfl_target = 0.3, t_stop = 1.5 (~7.5 diffusive times t_d = 0.2)

Compare with compare_analytic.py.
"""

import json

# Channel / fluid parameters
L_x = 0.2
L_y = 0.2
rho = 1.0
pres = 10.0
K = 5.0e-2  # n = 1 -> K is the plain dynamic viscosity mu
nn = 1.0
tau0 = 4.0e-3
g_x = 0.1

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
            "t_save": 0.15,
            "t_stop": 1.5,
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
            # Fluids Physical Parameters: single Bingham (HB with n = 1, tau0 > 0) fluid
            "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": 1.0 / K,
            "fluid_pp(1)%non_newtonian": "T",
            "fluid_pp(1)%K": K,
            "fluid_pp(1)%nn": nn,
            "fluid_pp(1)%tau0": tau0,
            "fluid_pp(1)%hb_m": 1.0e4,
            "fluid_pp(1)%mu_min": 1e-6,
            "fluid_pp(1)%mu_max": 1.0,
        }
    )
)
