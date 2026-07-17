#!/usr/bin/env python3
# Real-world 2D JWL detonation with full products-air mixing.
#
# A PETN charge (Kuhl et al., UCRL-PROC-225822 (2006) Appendix JWL fit, rho0 = 1000
# kg/m^3 -- Kuhl's own pressed-booster reference density, radius 5 cm) sits at rest at
# ambient pressure in free air at STP. A small central hot spot ignites a
# self-propagating JWL++ (jwl_reactive, Souers 2000) detonation: dl/dt = G p^b (1-l)
# releases the detonation energy E0/rho0 as the reaction front sweeps outward at the
# PETN CJ speed. The reaction is gated by Y*(1-lambda), so ONLY the
# unreacted explosive burns -- the surrounding air never reacts.
#
# Two features distinguish this from the shipped program-burn examples/2D_jwl_detonation:
#
#   1. Genuine reactive burn with ZND structure. fluid_pp(1)%jwl_delta_e applies the
#      Garno et al. (2020) Eq. 17 energy offset e_eff = e + Y*(1-lambda)*delta_e, so the
#      UNREACTED charge sits on a stiffer Hugoniot than the products (the Y factor keeps the
#      offset on the explosive, leaving the surrounding air untouched). delta_e is chosen
#      (below) so the unreacted PETN at (rho0, e=0) sits exactly at ambient pressure. The
#      detonation therefore carries a von Neumann pressure spike decaying to the CJ state
#      -- not the monotonic energy-source profile of a program burn. (The sub-mm reaction
#      zone is only partially resolved at this engineering resolution; the fully resolved
#      spike-vs-analytic validation lives in examples/1D_jwl_znd_detonation.)
#
#   2. Full products-air mixing. Once the detonation reaches the charge edge the products
#      (Y=1) expand violently into the air (Y=0), and the whole 0 < Y < 1 range of the
#      weighted-composition mixture closure is exercised across the contact -- the
#      Richtmyer-Meshkov-unstable products/air interface is where the mixture EOS does its
#      real work, rather than only the pure endpoints.
#
# Uniform 400x400 mesh (grid stretching is known to go NaN with strong JWL reactive
# sources -- see the note in examples/2D_jwl_detonation).
import json

from mfc.jwl_products import AIR, PETN_KUHL, ambient_fluid, jwl_fluid, znd_delta_e

# --- PETN JWL products (Kuhl et al. 2006 fit; see toolchain/mfc/jwl_products.py) ---
jwl_rho0 = PETN_KUHL["rho0"]

# --- Reactant/product energy offset (Garno Eq. 17) ---
# Placed so unreacted PETN (lambda=0) at (rho0, e=0) sits at ambient pressure, i.e.
# on a stiffer Hugoniot than the products -> resolvable ZND (von Neumann) structure.
p_amb = 101325.0
jwl_delta_e = znd_delta_e(PETN_KUHL, p_amb)

# --- JWL++ rate: dl/dt = jwl_G * p^jwl_b_exp * (1 - lambda) ---
# Calibrated in 1D (examples/1D_jwl_znd_detonation-style planar runs, see that
# README) against PETN_KUHL's measured, G-invariant CJ eigenvalue: D_CJ = 6245
# m/s, P_CJ = 16.7 GPa (front speed and post-reaction plateau pressure held
# fixed across a 3.3x change in G -- the operational definition of the CJ
# state). This G gives a ~1.25 mm reaction zone at 1D resolution (dx=10 um);
# at this case's 1 mm mesh the zone is sub-cell, same engineering-resolution
# caveat as the original TNT case (P10 Souers: recalibrate per grid).
jwl_G = 1.0e-13
jwl_b_exp = 2.0

# --- Air (ideal gas, STP) ---
air_rho = AIR["rho0"]
air_gamma_mfc = AIR["gamma"]
air_Cv = AIR["cv"]

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": -0.2,
            "x_domain%end": 0.2,
            "y_domain%beg": -0.2,
            "y_domain%end": 0.2,
            "m": 399,
            "n": 399,
            "p": 0,
            "dt": 1.5e-8,
            "t_step_start": 0,
            "t_step_stop": 4000,
            "t_step_save": 250,
            "num_patches": 3,
            "model_eqns": 2,
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # ParaView viz arrays: pressure/velocity/alpha come from prim_vars_wrt; these
            # add the JWL fields (T, Y_products, lambda), sound speed, z-vorticity, and a
            # numerical Schlieren of the detonation front.
            "jwl_wrt": "T",
            "c_wrt": "T",
            "omega_wrt(3)": "T",
            "schlieren_wrt": "T",
            "schlieren_alpha(1)": 1.0,
            "schlieren_alpha(2)": 1.0,
            "fd_order": 2,
            # Self-propagating JWL++ reactive burn (lambda initialized to 0 everywhere).
            "jwl_reactive": "T",
            "jwl_G": jwl_G,
            "jwl_b_exp": jwl_b_exp,
            # Patch 1: ambient air background.
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 0.4,
            "patch_icpp(1)%length_y": 0.4,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p_amb,
            "patch_icpp(1)%alpha_rho(1)": 1.63e-5,
            "patch_icpp(1)%alpha_rho(2)": air_rho * (1.0 - 1.0e-8),
            "patch_icpp(1)%alpha(1)": 1.0e-8,
            "patch_icpp(1)%alpha(2)": 0.99999999,
            # Patch 2: PETN charge, radius 5 cm, unreacted (ambient pressure).
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%radius": 0.05,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": p_amb,
            "patch_icpp(2)%alpha_rho(1)": jwl_rho0 * (1.0 - 1.0e-8),
            "patch_icpp(2)%alpha_rho(2)": air_rho * 1.0e-8,
            "patch_icpp(2)%alpha(1)": 0.99999999,
            "patch_icpp(2)%alpha(2)": 1.0e-8,
            # Patch 3: central hot spot (still PETN material) igniting the burn.
            "patch_icpp(3)%geometry": 2,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%alter_patch(2)": "T",
            "patch_icpp(3)%x_centroid": 0.0,
            "patch_icpp(3)%y_centroid": 0.0,
            "patch_icpp(3)%radius": 0.008,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%pres": 25.0e9,
            "patch_icpp(3)%alpha_rho(1)": jwl_rho0 * (1.0 - 1.0e-8),
            "patch_icpp(3)%alpha_rho(2)": air_rho * 1.0e-8,
            "patch_icpp(3)%alpha(1)": 0.99999999,
            "patch_icpp(3)%alpha(2)": 1.0e-8,
            # Direct initiation: a pressure-only hot spot was found (in 1D trials)
            # to be sub-critical for PETN_KUHL's lower CJ pressure -- it fizzles
            # rather than transitioning to a self-sustained detonation. Seeding the
            # booster as already-reacted (rxn_val = 1) gives robust direct initiation.
            "patch_icpp(3)%rxn_val": 1.0,
            # Fluid 1: PETN JWL products (+ ambient references); Fluid 2: air.
            **jwl_fluid(1, PETN_KUHL, AIR, delta_e=jwl_delta_e),
            **ambient_fluid(2, AIR),
        }
    )
)
