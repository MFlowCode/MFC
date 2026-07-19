#!/usr/bin/env python3
# 3D inert spherical products burst into air. A full-density TNT-products sphere
# (radius 3 cm) is initialized already detonated, held at a post-CJ pressure, and
# released instantaneously into ambient air. There is no program burn or reaction
# (that is the Tier-2 reactive layer); this case exercises the JWL products/air
# closure and the HLLC/cons-prim path in 3D as the products expand and drive a
# spherical blast into the air. It is the exascale-structure smoke test: the same
# case runs on-device under the GPU build with no code change.
#
# Sized for a ~5 min run on 8 CPU cores (run with --mpi -n 8): 128^3 mesh, 400 steps,
# parallel_io. Non-reactive, extrapolation BCs (-3, since CBC is prohibited with eos_jwl).
import json

L = 0.2  # cube 0.2 m on a side ([-0.1, 0.1] m)

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": -0.5 * L,
            "x_domain%end": 0.5 * L,
            "y_domain%beg": -0.5 * L,
            "y_domain%end": 0.5 * L,
            "z_domain%beg": -0.5 * L,
            "z_domain%end": 0.5 * L,
            "m": 127,
            "n": 127,
            "p": 127,
            "dt": 1.9e-8,
            "t_step_start": 0,
            "t_step_stop": 400,
            "t_step_save": 80,
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
            "bc_z%beg": -3,
            "bc_z%end": -3,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: ambient air (cuboid filling the domain).
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": L,
            "patch_icpp(1)%length_z": L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 101325.0,
            "patch_icpp(1)%alpha_rho(1)": 1.63e-5,
            "patch_icpp(1)%alpha_rho(2)": 1.22499998775,
            "patch_icpp(1)%alpha(1)": 1.0e-8,
            "patch_icpp(1)%alpha(2)": 0.99999999,
            # Patch 2: full-density products sphere held at a post-CJ pressure so it
            # expands instantaneously (inert).
            "patch_icpp(2)%geometry": 8,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%z_centroid": 0.0,
            "patch_icpp(2)%radius": 0.03,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
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
