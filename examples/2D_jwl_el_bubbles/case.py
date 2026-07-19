#!/usr/bin/env python3
# Inert products blast into a dispersed Euler-Lagrange bubble cloud. A full-density
# TNT-products charge on the y=0 axis is released instantaneously (no program burn,
# no reaction) and drives a shock into ambient air; the shock reaches a small cloud
# of sub-grid Lagrangian bubbles placed downstream. The bubbles are seeded and driven
# through the JWL closure without an EOS branch: the seed pressure comes from
# s_jwl_mix_state_er (m_bubbles_EL.fpp s_add_bubbles) and the Keller-Miksis liquid
# sound speed from s_jwl_mix_sound_speed, while the driving pressure is q_prim(E),
# already built through the JWL closure in cons to prim.
#
# The bubbles sit well downstream (x about 0.05, roughly 0.11 m from the charge) so the
# blast has decayed into the sub-singular Keller-Miksis envelope before it arrives;
# adap_dt subcycles the stiff radial ODE inside the fluid step. Four bubbles in
# y-symmetric pairs (y = +/- 0.008) keep a correct run mirror symmetric about y = 0.
# Needs input/lag_bubbles.dat, so this is an example, not a golden-suite test.
# Non-reactive, extrapolation BCs (-3, since CBC is prohibited with eos_jwl).
import json

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": -0.1,
            "x_domain%end": 0.1,
            "y_domain%beg": -0.1,
            "y_domain%end": 0.1,
            "m": 149,
            "n": 149,
            "p": 0,
            "dt": 2.5e-8,
            "t_step_start": 0,
            "t_step_stop": 1200,
            "t_step_save": 100,
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
            "parallel_io": "T",
            "lag_db_wrt": "T",
            # Euler-Lagrange bubble cloud (sub-grid). Keller-Miksis with full heat/mass transfer.
            "bubbles_lagrange": "T",
            "bubble_model": 2,
            "thermal": 3,
            "polytropic": "F",
            "adap_dt": "T",
            "lag_params%nBubs_glb": 4,
            "lag_params%solver_approach": 2,
            "lag_params%cluster_type": 2,
            "lag_params%pressure_corrector": "T",
            "lag_params%smooth_type": 1,
            "lag_params%heatTransfer_model": "T",
            "lag_params%massTransfer_model": "T",
            "lag_params%epsilonb": 1.0,
            "lag_params%valmaxvoid": 0.9,
            "lag_params%write_bubbles": "F",
            "lag_params%write_bubbles_stats": "F",
            # Bubble gas and interface properties (dimensional SI).
            "bub_pp%R0ref": 1.0,
            "bub_pp%p0ref": 101325.0,
            "bub_pp%rho0ref": 1.225,
            "bub_pp%T0ref": 298.0,
            "bub_pp%ss": 0.069,
            "bub_pp%pv": 2350.0,
            "bub_pp%vd": 2.5e-5,
            "bub_pp%mu_l": 1.0e-3,
            "bub_pp%gam_v": 1.333,
            "bub_pp%gam_g": 1.4,
            "bub_pp%M_v": 18.0,
            "bub_pp%M_g": 28.0,
            "bub_pp%k_v": 0.02,
            "bub_pp%k_g": 0.025,
            "bub_pp%cp_v": 2100.0,
            "bub_pp%cp_g": 1000.0,
            "bub_pp%R_v": 461.9,
            "bub_pp%R_g": 296.9,
            # Patch 1: ambient air at rest, fills the domain.
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 0.2,
            "patch_icpp(1)%length_y": 0.2,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 101325.0,
            "patch_icpp(1)%alpha_rho(1)": 1.63e-5,
            "patch_icpp(1)%alpha_rho(2)": 1.22499998775,
            "patch_icpp(1)%alpha(1)": 1.0e-8,
            "patch_icpp(1)%alpha(2)": 0.99999999,
            # Patch 2: full-density products charge on the axis, held at a post-CJ
            # pressure so it expands instantaneously (inert).
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -0.055,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%radius": 0.02,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 8.5e9,
            "patch_icpp(2)%alpha_rho(1)": 1629.99998,
            "patch_icpp(2)%alpha_rho(2)": 1.225e-8,
            "patch_icpp(2)%alpha(1)": 0.99999999,
            "patch_icpp(2)%alpha(2)": 1.0e-8,
            # Fluid 1: TNT JWL products.
            "fluid_pp(1)%eos": 2,
            "fluid_pp(1)%gamma": 2.5,
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%cv": 613.5,
            "fluid_pp(1)%jwl_A": 3.712e11,
            "fluid_pp(1)%jwl_B": 3.231e9,
            "fluid_pp(1)%jwl_R1": 4.15,
            "fluid_pp(1)%jwl_R2": 0.95,
            "fluid_pp(1)%jwl_omega": 0.30,
            "fluid_pp(1)%jwl_rho0": 1630.0,
            "fluid_pp(1)%jwl_E0": 1.0089e10,
            "fluid_pp(1)%jwl_air_e0": 2.5575e5,
            "fluid_pp(1)%jwl_air_rho0": 1.225,
            # Fluid 2: air (also the bubble gas EOS reference).
            "fluid_pp(2)%eos": 1,
            "fluid_pp(2)%gamma": 2.5,
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%cv": 717.5,
        }
    )
)
