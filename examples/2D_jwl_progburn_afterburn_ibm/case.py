#!/usr/bin/env python3
# Afterburn combustion around a free-flight immersed body: program burn -> JWL products ->
# products/air afterburn mixing -> blast shock -> moving IB cylinder.
#
# Same radial program-burn detonation as 2D_jwl_progburn_ibm_freeflight, with the mixing-rate
# afterburn model layered on top (jwl_ab_model = 1): as products mix with air behind the shock,
# progress b relaxes toward 1 and releases the secondary combustion budget Y*jwl_q_ab.
#
# This is the regression anchor for reactive/afterburn JWL + IBM support. The afterburn progress
# b (eqn_idx%abn) is a material scalar, so IBM ghost cells take the image-point value (zero
# normal gradient at a rigid wall, the same constant-extrapolation treatment as alpha) and
# moving-IB fresh cells take the same inverse-distance neighbour average as the other fields
# (both in m_ibm.fpp s_ibm_correct_state). The reaction sources skip solid cells entirely
# (m_jwl_sources.fpp), so no phantom energy is released inside the body.
#
# The setup is mirror-symmetric about y = 0, so a correct run keeps b, rho, and the body on the
# axis to round-off; b stays in [0, 1]; and the body ends slightly faster than the afterburn-off
# case since the secondary release strengthens the blast.
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
            "ib_state_wrt": "T",
            # Program burn: charge center-initiated on the y=0 axis, radial front.
            "prog_burn": "T",
            "jwl_afterburn": "T",
            "jwl_ab_model": 1,
            "jwl_q_ab": 3.0e6,
            "jwl_ab_tau": 3.0e-6,
            "pb_D_cj": 6930.0,
            "pb_width": 0.008,
            "pb_x_det": -0.055,
            "pb_y_det": 0.0,
            # Free-flight rigid cylinder: starts at rest on the y=0 axis, mass finite so the
            # blast accelerates it. moving_ibm = 2 is two-way coupling (fluid force -> body).
            "ib": "T",
            "num_ibs": 1,
            "fd_order": 2,
            "patch_ib(1)%geometry": 2,
            "patch_ib(1)%x_centroid": 0.0,
            "patch_ib(1)%y_centroid": 0.0,
            "patch_ib(1)%radius": 0.015,
            "patch_ib(1)%slip": "F",
            "patch_ib(1)%moving_ibm": 2,
            "patch_ib(1)%vel(1)": 0.0,
            "patch_ib(1)%vel(2)": 0.0,
            "patch_ib(1)%angles(1)": 0.0,
            "patch_ib(1)%angles(2)": 0.0,
            "patch_ib(1)%angles(3)": 0.0,
            "patch_ib(1)%angular_vel(1)": 0.0,
            "patch_ib(1)%angular_vel(2)": 0.0,
            "patch_ib(1)%angular_vel(3)": 0.0,
            "patch_ib(1)%mass": 5.0e-3,
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
            # Patch 2: TNT-products charge, reduced density so the unburned state sits at ambient
            # pressure; the burn front deposits the detonation energy. Centered on the axis.
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -0.055,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%radius": 0.02,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 101325.0,
            "patch_icpp(2)%alpha_rho(1)": 99.999999,
            "patch_icpp(2)%alpha_rho(2)": 1.225e-8,
            "patch_icpp(2)%alpha(1)": 0.99999999,
            "patch_icpp(2)%alpha(2)": 1.0e-8,
            # Fluid 1: TNT JWL products (verbatim from examples/2D_jwl_detonation).
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
            # Fluid 2: air.
            "fluid_pp(2)%eos": 1,
            "fluid_pp(2)%gamma": 2.5,
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%cv": 717.5,
        }
    )
)
