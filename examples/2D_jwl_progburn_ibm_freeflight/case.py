#!/usr/bin/env python3
# End-to-end detonation chain into a free-flight immersed body:
#   program burn -> full JWL products -> products/air mixing -> blast shock -> moving IB particle.
#
# This replaces the static high-pressure driver slab used in 2D_jwl_ibm_freeflight_cylinder
# with a real program burn. A reduced-density TNT-products charge (fluid 1 at about 100 kg/m3,
# sitting at ambient pressure until it is lit) is center-initiated on the y=0 axis. A radial
# lighting-time front at pb_D_cj releases the detonation energy (prog_burn). The products expand,
# drive a shock into the ambient air, and the shock strikes a rigid cylinder that starts at rest.
# The pressure load accelerates the body downstream.
#
# Why each stage is trustworthy here:
#   burn: the front is kinematic, radius = pb_D_cj*(t - t_det); pb_D_cj*dt << pb_width so no
#         cells are skipped, and pb_width spans about 6 cells so the band is resolved.
#   JWL + mixing: the composition-weighted closure is exact at Y=0 (air) and Y=1 (products);
#         HLLC keeps the products/air contact and Y sharp.
#   shock -> body: the blast load is the pressure-gradient integral of q_prim(E)
#         (s_compute_ib_forces), already built through the JWL closure in the cons to prim
#         conversion, so a JWL cell feeds the body a JWL-consistent force with no EOS branch.
#         Newton's law then advances the body (m_time_steppers s_propagate_immersed_boundaries).
#
# The charge (centroid -0.055, radius 0.02) and the cylinder (centroid 0, radius 0.015) do not
# overlap, so the burn completes entirely in fluid and the front never enters the solid. prog_burn
# carries no advected reaction-progress variable, so the IBM ghost/fresh-cell rebuild (which fills
# cont, adv, mom, E) is sufficient and the ib + reactive prohibition does not apply.
#
# The setup is mirror-symmetric about y = 0 (charge and detonation point on the axis, cylinder
# centered on the axis). A correct run keeps the body on the axis: it translates in +x only, with
# zero net y-force and zero net torque, and the flow stays symmetric about y = 0 to round-off.
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
            "ib_state_wrt": "T",
            # Program burn: charge center-initiated on the y=0 axis, radial front.
            "prog_burn": "T",
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
