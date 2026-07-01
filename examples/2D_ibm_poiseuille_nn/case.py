#!/usr/bin/env python3
"""
2D power-law Poiseuille channel whose no-slip walls are IMMERSED BOUNDARIES,
validating the IBM + non-Newtonian viscosity interaction (the per-stencil-sample
Herschel-Bulkley viscosity used by the IBM ghost-point and force machinery).

Two rectangular IB slabs form the channel: the flow gap is y in [0.05, 0.25]
(half-height H = 0.1, centerline y = 0.15) inside a y-domain [0, 0.3]. Each slab
extends half outside the domain (and beyond the domain in x) so its only surface
seen by the flow is the flat face bounding the gap; the domain BCs behind the
slabs are no-slip walls. A constant body acceleration g_x drives the flow in the
gap; x is periodic. The steady gap profile is compared against the closed-form
power-law solution (see examples/2D_poiseuille_thickening_nn):

    u(y) = (n/(n+1)) * (rho*g/K)^(1/n) * ( H^((n+1)/n) - |y - y_c|^((n+1)/n) )

Fluid/forcing parameters match examples/2D_poiseuille_thickening_nn (K = 5e-2,
nn = 1.5, g_x = 5e-2, mu_wall = K^(1/n)*(rho*g*H)^((n-1)/n) ~ 0.0232,
mu_max = 0.035 ~ 1.5*mu_wall so the clamp stays inactive). The gap is resolved
with 64 cells (dy = 0.003125), the same resolution as the BC-walled template.

Modes (environment variable IBM_NN_MODE, default "powerlaw"):
    powerlaw  : committed config; nn = 1.5 power-law fluid, cfl_adap_dt to
                t_stop = 0.9 (~2 wall-viscous diffusion times). Validation B.
    newtonian : constant-viscosity fluid, mu = 0.02 (Re(1) = 50), FIXED
                dt = 6e-5 to t = 0.3. Validation A run (i).
    nn1       : degenerate non-Newtonian fluid, nn = 1.0, tau0 = 0, K = 0.02
                (analytically the same fluid as "newtonian"), same fixed dt.
                Validation A run (ii). Fields must match run (i) to round-off.

See README.md and compare_analytic.py / check_equivalence.py.
"""

import json
import os

MODE = os.environ.get("IBM_NN_MODE", "powerlaw")
assert MODE in ("powerlaw", "newtonian", "nn1"), f"bad IBM_NN_MODE: {MODE}"

# Channel / fluid parameters
L_x = 0.2
L_y = 0.3  # domain height; flow gap is y in [0.05, 0.25]
# Slab full height: most of it lies outside the domain. The centroids sit just
# INSIDE the domain (not exactly on the boundary, where no rank would own the
# IB and its ib_state force record would never be written) while the gap faces
# stay exactly at y = 0.05 and y = 0.25 (on cell boundaries).
slab = 0.099
rho = 1.0
pres = 10.0
g_x = 5.0e-2

case = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    "x_domain%beg": 0.0,
    "x_domain%end": L_x,
    "y_domain%beg": 0.0,
    "y_domain%end": L_y,
    "m": 24,
    "n": 95,
    "p": 0,
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
    # x: periodic channel; y: no-slip walls behind the IB slabs
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -16,
    "bc_y%end": -16,
    "viscous": "T",
    # Constant body acceleration in +x
    "bf_x": "T",
    "g_x": g_x,
    "k_x": 0.0,
    "w_x": 0.0,
    "p_x": 0.0,
    # Immersed-boundary channel walls: two rectangle slabs. Each extends beyond
    # the domain in x and half outside the domain in y, so the only IB surface
    # the flow sees is the flat face at y = 0.05 (slab 1) / y = 0.25 (slab 2).
    # mass = 0 zeroes the bf_x*mass bookkeeping term in the reported IB force,
    # so ib_state files hold the pure pressure+viscous surface integration.
    "ib": "T",
    "num_ibs": 2,
    "ib_state_wrt": "T",
    "patch_ib(1)%geometry": 3,
    "patch_ib(1)%x_centroid": 0.5 * L_x,
    "patch_ib(1)%y_centroid": 0.05 - 0.5 * slab,  # gap face at y = 0.05
    "patch_ib(1)%length_x": 2.0 * L_x,
    "patch_ib(1)%length_y": slab,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%mass": 0.0,
    "patch_ib(2)%geometry": 3,
    "patch_ib(2)%x_centroid": 0.5 * L_x,
    "patch_ib(2)%y_centroid": 0.25 + 0.5 * slab,  # gap face at y = 0.25
    "patch_ib(2)%length_x": 2.0 * L_x,
    "patch_ib(2)%length_y": slab,
    "patch_ib(2)%slip": "F",
    "patch_ib(2)%mass": 0.0,
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
    # Fluids Physical Parameters (EOS shared by all modes)
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}

if MODE == "powerlaw":
    K = 5.0e-2
    case |= {
        "cfl_adap_dt": "T",
        "cfl_target": 0.3,
        "n_start": 0,
        "t_save": 0.09,
        "t_stop": 0.9,
        "fluid_pp(1)%Re(1)": 1.0 / K,
        "fluid_pp(1)%non_newtonian": "T",
        "fluid_pp(1)%K": K,
        "fluid_pp(1)%nn": 1.5,
        "fluid_pp(1)%tau0": 0.0,
        "fluid_pp(1)%hb_m": 1000.0,
        "fluid_pp(1)%mu_min": 1e-6,
        "fluid_pp(1)%mu_max": 0.035,
    }
else:
    # Equivalence pair: identical fixed dt so the two runs share an exact
    # timestep history and can be compared at matched physical time.
    mu = 0.02
    case |= {
        "dt": 6.0e-5,
        "t_step_start": 0,
        "t_step_stop": 5000,
        "t_step_save": 1000,
        "fluid_pp(1)%Re(1)": 1.0 / mu,
    }
    if MODE == "nn1":
        case |= {
            "fluid_pp(1)%non_newtonian": "T",
            "fluid_pp(1)%K": mu,
            "fluid_pp(1)%nn": 1.0,
            "fluid_pp(1)%tau0": 0.0,
            "fluid_pp(1)%hb_m": 1000.0,
            "fluid_pp(1)%mu_min": 1e-6,
            "fluid_pp(1)%mu_max": 0.03,
        }

print(json.dumps(case))
