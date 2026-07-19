#!/usr/bin/env python3
# Inert products blast onto a free-flight rigid cylinder. A full-density
# TNT-products charge on the y=0 axis is released instantaneously (no program
# burn, no reaction); the shock strikes a rigid cylinder that starts at rest and
# accelerates downstream under the pressure load. This exercises the moving-IB
# fresh-cell repopulation: cells uncovered as the body translates are rebuilt by
# inverse-distance extrapolation and closed through the JWL EOS, instead of being
# left at the near-vacuum solid seed. Charge and cylinder both sit on y=0, so a
# correct run keeps the body on the axis: +x translation only, zero net y-force
# and torque, flow symmetric about y=0 to round-off.
# moving_ibm = 2 is two-way coupling (fluid force -> body).
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
            "t_step_save": 200,
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
            # Free-flight rigid cylinder: starts at rest on the y=0 axis, finite mass
            # so the blast accelerates it.
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
            # Patch 1: ambient air at rest.
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
            # Fluid 2: air.
            "fluid_pp(2)%eos": 1,
            "fluid_pp(2)%gamma": 2.5,
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%cv": 717.5,
        }
    )
)
