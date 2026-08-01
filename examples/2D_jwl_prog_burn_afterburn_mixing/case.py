#!/usr/bin/env python3
# 2D free-air TNT detonation combining a kinematic program burn with mixing-rate
# afterburn, run through full products-air mixing.
#
# A full-density cylindrical TNT charge (rho0 = 1630 kg/m^3, radius 5 cm) sits at
# rest at ambient pressure in free air at STP. prog_burn ignites at the charge
# center and deposits the primary detonation energy jwl_E0/jwl_rho0 as a kinematic
# front sweeping outward at the TNT CJ speed (~6.93 km/s, radius = pb_D_cj*t
# exactly, by construction - unlike jwl_reactive the front is prescribed, not
# self-propagating). Once the charge is consumed the products (Y=1) expand
# violently into the surrounding air (Y=0), exercising the full 0 < Y < 1 range
# of the composition-weighted mixture closure across the Richtmyer-Meshkov-
# unstable products/air interface (see examples/2D_jwl_reactive_air_mixing for
# the jwl_reactive analogue of this same mixing physics).
#
# jwl_afterburn (mixing-rate model, jwl_ab_model = 1) is layered on top: as
# products (Y -> 1) contact unburned air (Y -> 0) at the interface, the advected
# afterburn progress b relaxes toward 1 with rate db/dt = (1-b)(1-Y)/jwl_ab_tau,
# releasing a SECOND energy budget jwl_q_ab (secondary combustion of detonation
# products with atmospheric oxygen - physically distinct from and released after
# the primary jwl_E0 the front already deposited). The clamp caps the release at
# Y*jwl_q_ab per cell (m_jwl_sources.fpp), so afterburn energy only appears where
# products have actually mixed into air (0 < Y < 1), not in the pure-air far
# field or the still-unburned charge interior. prog_burn and jwl_afterburn are
# NOT mutually exclusive (only jwl_reactive + prog_burn are);
# together they model the real two-stage TNT energy release: primary detonation
# (JWL fit energy) then afterburn (fireball combustion) as products mix with air.
#
# Uniform 400x400 mesh (grid stretching is known to go NaN with strong JWL
# sources - see the note in examples/2D_jwl_detonation).
import json

# --- TNT JWL products (same fit as examples/2D_jwl_detonation) ---
jwl_A = 3.712e11
jwl_B = 3.231e9
jwl_R1 = 4.15
jwl_R2 = 0.95
jwl_omega = 0.30
jwl_rho0 = 1630.0
jwl_E0 = 1.0089e10  # E0 = rho0 * Q_detonation
jwl_Cv = 613.5

# --- Program burn: kinematic front at the TNT CJ speed, ignited at the charge center ---
pb_D_cj = 6930.0
pb_width = 0.003  # 3 mm ~ 3 cells at dx = 1 mm; pb_D_cj*dt must stay <= pb_width

# --- Afterburn: mixing-rate model, gated to the products-air contact zone ---
# jwl_ab_tau set so the afterburn budget releases over a few microseconds behind
# the front, comparable to the time the products/air interface takes to develop
# after the charge is consumed (~ charge_radius / D_cj ~ 7 us).
jwl_q_ab = 3.0e6  # J/kg, secondary (fireball) combustion energy on top of jwl_E0
jwl_ab_tau = 3.0e-6

# --- Air (ideal gas, STP) ---
air_rho = 1.225
air_gamma_mfc = 2.5  # 1/(1.4 - 1)
air_Cv = 717.5
p_amb = 101325.0

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
            "dt": 5.0e-8,
            "t_step_start": 0,
            "t_step_stop": 600,
            "t_step_save": 50,
            "num_patches": 2,
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
            # Program burn: center initiation.
            "prog_burn": "T",
            "pb_D_cj": pb_D_cj,
            "pb_width": pb_width,
            "pb_x_det": 0.0,
            "pb_y_det": 0.0,
            # Afterburn: mixing-rate model, active only where products (Y) have
            # reached the interface with air.
            "jwl_afterburn": "T",
            "jwl_ab_model": 1,
            "jwl_q_ab": jwl_q_ab,
            "jwl_ab_tau": jwl_ab_tau,
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
            # Fluid 2: air.
            "fluid_pp(2)%eos": 1,
            "fluid_pp(2)%gamma": air_gamma_mfc,
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%cv": air_Cv,
        }
    )
)
