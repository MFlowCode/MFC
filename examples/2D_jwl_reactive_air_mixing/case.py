#!/usr/bin/env python3
# Real-world 2D JWL detonation with full products-air mixing.
#
# A full-density cylindrical TNT charge (rho0 = 1630 kg/m^3, radius 5 cm) sits at
# rest at ambient pressure in free air at STP. A small central hot spot ignites a
# self-propagating JWL++ (jwl_reactive, Souers 2000) detonation: dl/dt = G p^b (1-l)
# releases the detonation energy E0/rho0 as the reaction front sweeps outward at the
# TNT CJ speed (~6.9 km/s). The reaction is gated by Y*(1-lambda), so ONLY the
# unreacted explosive burns -- the surrounding air never reacts.
#
# Two features distinguish this from the shipped program-burn examples/2D_jwl_detonation:
#
#   1. Genuine reactive burn with ZND structure. fluid_pp(1)%jwl_delta_e applies the
#      Garno et al. (2020) Eq. 17 energy offset e_eff = e + Y*(1-lambda)*delta_e, so the
#      UNREACTED charge sits on a stiffer Hugoniot than the products (the Y factor keeps the
#      offset on the explosive, leaving the surrounding air untouched). delta_e is chosen
#      (below) so the unreacted TNT at (rho0, e=0) sits exactly at ambient pressure. The
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
import math

# --- TNT JWL products (same fit as examples/2D_jwl_detonation) ---
jwl_A = 3.712e11
jwl_B = 3.231e9
jwl_R1 = 4.15
jwl_R2 = 0.95
jwl_omega = 0.30
jwl_rho0 = 1630.0
jwl_E0 = 1.0089e10  # E0 = rho0 * Q_detonation
jwl_Cv = 613.5

# --- Reactant/product energy offset (Garno Eq. 17) ---
# Placed so unreacted TNT (lambda=0) at (rho0, e=0) sits at ambient pressure, i.e.
# on a stiffer Hugoniot than the products -> resolvable ZND (von Neumann) structure.
p_amb = 101325.0
F1 = jwl_A * (1.0 - jwl_omega / jwl_R1) * math.exp(-jwl_R1) + jwl_B * (1.0 - jwl_omega / jwl_R2) * math.exp(-jwl_R2)
jwl_delta_e = (p_amb - F1) / (jwl_omega * jwl_rho0)  # ~ -1.29e7 J/kg

# --- JWL++ rate: dl/dt = jwl_G * p^jwl_b_exp * (1 - lambda) ---
# G set so the reaction zone spans a few 1 mm cells at TNT's ~21 GPa CJ pressure
# (e-fold time ~0.4 us -> zone ~3 mm), giving a sharp but resolved detonation front.
jwl_G = 5.0e-15
jwl_b_exp = 2.0

# --- Air (ideal gas, STP) ---
air_rho = 1.225
air_gamma_mfc = 2.5  # 1/(1.4 - 1)
air_Cv = 717.5

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
            # Patch 2: full-density TNT charge, radius 5 cm, unreacted (ambient pressure).
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
            # Patch 3: central hot spot (still TNT material) igniting the burn.
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
            # Fluid 1: TNT JWL products.
            "fluid_pp(1)%eos": 2,
            "fluid_pp(1)%gamma": 2.5,
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%cv": jwl_Cv,
            "fluid_pp(1)%jwl_A": jwl_A,
            "fluid_pp(1)%jwl_B": jwl_B,
            "fluid_pp(1)%jwl_R1": jwl_R1,
            "fluid_pp(1)%jwl_R2": jwl_R2,
            "fluid_pp(1)%jwl_omega": jwl_omega,
            "fluid_pp(1)%jwl_rho0": jwl_rho0,
            "fluid_pp(1)%jwl_E0": jwl_E0,
            "fluid_pp(1)%jwl_air_e0": 2.5575e5,
            "fluid_pp(1)%jwl_air_rho0": air_rho,
            "fluid_pp(1)%jwl_delta_e": jwl_delta_e,
            # Fluid 2: air.
            "fluid_pp(2)%eos": 1,
            "fluid_pp(2)%gamma": air_gamma_mfc,
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%cv": air_Cv,
        }
    )
)
