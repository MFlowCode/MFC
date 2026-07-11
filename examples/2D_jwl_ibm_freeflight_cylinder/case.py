#!/usr/bin/env python3
# Free-flight rigid cylinder driven by a JWL-products blast (Phase 1c validation).
#
# Purpose: exercise the two-way-coupled moving immersed boundary (moving_ibm = 2) in a JWL
# field. A high-pressure TNT-products driver slab on the left launches a shock rightward. The
# shock strikes a rigid cylinder that starts at rest, and the pressure load accelerates the
# body downstream. The blast load is the pressure-gradient integral of the primitive pressure
# q_prim(E) (s_compute_ib_forces in m_ibm.fpp), which is already built through the JWL closure
# in the cons to prim conversion, so a JWL cell feeds the body a JWL-consistent force with no
# separate EOS branch needed. Newton's law then advances the body (m_time_steppers.fpp
# s_propagate_immersed_boundaries): vel += rk*dt*force/mass.
#
# The configuration is mirror-symmetric about y = 0 (driver spans the full y extent, cylinder
# centered on the axis). A correct run therefore keeps the body on the axis: it translates in
# +x only, with zero net y-force and zero net torque, and the flow stays symmetric about y = 0
# to round-off. The 2D torque channel is a known base-code limitation (tau_z is not accumulated
# in 2D), which this symmetric setup makes irrelevant since the physical torque is zero.
#
# As the body moves, its trailing (-x) edge uncovers solid cells; those fresh cells are
# repopulated through the JWL closure by the moving-IB branch of s_ibm_correct_state. Non
# reactive, extrapolation BCs (-3, since CBC is prohibited with eos_jwl).
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
            "t_step_stop": 600,
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
            # Patch 1: ambient products at rest (low pressure), fills the domain.
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 0.2,
            "patch_icpp(1)%length_y": 0.2,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 101325.0,
            "patch_icpp(1)%alpha_rho(1)": 99.999999,
            "patch_icpp(1)%alpha_rho(2)": 1.225e-8,
            "patch_icpp(1)%alpha(1)": 0.99999999,
            "patch_icpp(1)%alpha(2)": 1.0e-8,
            # Patch 2: high-pressure products driver slab on the left (x in [-0.1, -0.03]).
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -0.065,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%length_x": 0.07,
            "patch_icpp(2)%length_y": 0.2,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 5.0e8,
            "patch_icpp(2)%alpha_rho(1)": 299.999997,
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
