#!/usr/bin/env python3
"""
2D general Herschel-Bulkley (flow index n != 1 AND yield stress tau0 > 0)
Poiseuille channel.

This case combines both non-Newtonian effects at once: a shear-thinning power-law
term (nn = 0.5) AND a yield stress (tau0 > 0). The companion examples isolate each
effect -- 2D_poiseuille_nn / 2D_poiseuille_thickening_nn (power-law only, tau0 = 0)
and 2D_bingham_poiseuille_nn (yield only, n = 1). The signature of a correct general
Herschel-Bulkley model is a rigid PLUG near the centerline (where |tau| < tau0)
joined to a shear-thinning sheared profile near the walls.

A constant body acceleration g_x drives a fully-developed channel flow between two
stationary no-slip walls (y = 0, L_y; half-height H = L_y/2, centerline y = H).
With tau_w = rho*g*H > tau0, the steady Herschel-Bulkley profile is

  plug half-width from centerline   : y0   = tau0/(rho*g)
  sheared region (0 <= y <= H - y0) : u(y) = (n/((n+1)*rho*g)) * K^(-1/n) *
                                        [ (tau_w - tau0)^((n+1)/n)
                                          - (rho*g*(H-y) - tau0)^((n+1)/n) ]
  plug (H - y0 <= y <= H)           : u_plug = (n/((n+1)*rho*g)) * K^(-1/n) *
                                        (tau_w - tau0)^((n+1)/n)
  upper half mirrors about y = H. Requires tau_w = rho*g*H > tau0 for any flow.

Parameters (nondimensional MFC units):
    rho  = 1.0, H = 0.1 (L_y = 0.2)
    nn   = 0.5          (shear-thinning flow index)
    g_x  = 0.1          -> tau_w = rho*g*H = 1.0e-2
    tau0 = 3.5e-3       -> y0 = tau0/(rho*g) = 3.5e-2 = 0.35 H  (clear plug + shear)
    K    = 1.5e-2       (consistency index)
    hb_m = 1.0e4        (sharp Papanastasiou yield regularization)
    mu_min = 1e-6, mu_max = 0.3  (caps the plug viscosity, which diverges as the shear
          rate -> 0 in the plug; mu_max = 0.3 is ~6x the wall effective viscosity, so
          the plug is effectively rigid. The explicit viscous timestep dt ~ dy^2
          rho/mu_max is the binding cost: the previous mu_max = 1.0 (the Bingham scale)
          drove dt to ~3e-6 and also relaxed to steady state too slowly to finish in
          minutes; mu_max = 0.3 gives dt ~ 1e-5 and a faster transient while still
          resolving a clear plug. A much smaller mu_max lets the near-plug region
          shear too freely and shrinks the plug.)
    pres = 10 -> sound speed ~3.74; u_plug ~ 4e-3 => Mach ~1e-3 (low Mach)
    grid: m = 24 (x, periodic; >= 24 for a 2-rank y-split WENO5 decomposition),
          n = 63 (y resolution of plug + sheared layer), L_x = L_y = 0.2
    cfl_adap_dt, cfl_target = 0.3, t_stop = 0.4 (the velocity profile is steady to
          ~1% between the last two saves; the regularized plug is not perfectly flat,
          so the >=95% u_max band is the meaningful plug width, matching y0)

Compare with compare_analytic.py.
"""

import json

# Channel / fluid parameters
L_x = 0.2
L_y = 0.2
rho = 1.0
pres = 10.0
K = 1.5e-2
nn = 0.5
tau0 = 3.5e-3
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
            "t_save": 0.04,
            "t_stop": 0.4,
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
            # Fluids Physical Parameters: single Herschel-Bulkley (n != 1, tau0 > 0) fluid
            "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": 1.0 / K,
            "fluid_pp(1)%non_newtonian": "T",
            "fluid_pp(1)%K": K,
            "fluid_pp(1)%nn": nn,
            "fluid_pp(1)%tau0": tau0,
            "fluid_pp(1)%hb_m": 1.0e4,
            "fluid_pp(1)%mu_min": 1e-6,
            "fluid_pp(1)%mu_max": 0.3,
        }
    )
)
