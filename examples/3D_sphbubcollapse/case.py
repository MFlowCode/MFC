#!/usr/bin/env python3
import json

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            # For these computations, the bubble is placed at the (0,0,0)
            # domain origin. Only one octant of a full spherical bubble is
            # used, with reflecting BC's at octant boundaries acting
            # as the other bubble regions.
            "x_domain%beg": 0.0e00,
            "x_domain%end": 4.0e-03 / 1.0e-03,
            "y_domain%beg": 0.0e00,
            "y_domain%end": 4.0e-03 / 1.0e-03,
            "z_domain%beg": 0.0e00,
            "z_domain%end": 4.0e-03 / 1.0e-03,
            # Grid stretching is used in the all coordinate directions
            # to minimize computational costs. The grid is coarsened
            # away from the bubble / origin
            "stretch_x": "T",
            "a_x": 4.0e00,
            "x_a": -1.5e-03 / 1.0e-03,
            "x_b": 1.5e-03 / 1.0e-03,
            "stretch_y": "T",
            "a_y": 4.0e00,
            "y_a": -1.5e-03 / 1.0e-03,
            "y_b": 1.5e-03 / 1.0e-03,
            "stretch_z": "T",
            "a_z": 4.0e00,
            "z_a": -1.5e-03 / 1.0e-03,
            "z_b": 1.5e-03 / 1.0e-03,
            "cyl_coord": "F",
            "m": 99,
            "n": 99,
            "p": 99,
            "dt": 0.2e-09 / 1.0e-03,
            "t_step_start": 0,
            "t_step_stop": 133300,
            "t_step_save": 100,
            # Simulation Algorithm Parameters
            # Only two patches are necessary, the background liquid and the
            # gas bubble
            "num_patches": 2,
            # Use the 6 equation model
            "model_eqns": 3,
            # 6 equations model does not need the K \div(u) term
            "alt_soundspeed": "F",
            # Two fluids: water and air
            "num_fluids": 2,
            # time step
            "mpp_lim": "T",
            # Correct errors when computing speed of sound
            "mixture_err": "T",
            # Use TVD RK3 for time marching
            "time_stepper": 3,
            # Use WENO5
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "avg_state": 2,
            # Use the mapped WENO weights to maintain monotinicity
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            # Use the HLLC  Riemann solver
            "riemann_solver": 2,
            "wave_speeds": 1,
            # We use reflective boundary conditions at octant edges and
            # non-reflective boundary conditions at the domain edges
            "bc_x%beg": -2,
            "bc_x%end": -6,
            "bc_y%beg": -2,
            "bc_y%end": -6,
            "bc_z%beg": -2,
            "bc_z%end": -6,
            # Formatted Database Files Structure Parameters
            # Export primitive variables in double precision with parallel
            # I/O to minimize I/O computational time during large simulations
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: High pressured water
            # Specify the cubic water background grid geometry
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 80.0e-03 / 1.0e-03,
            "patch_icpp(1)%y_centroid": 80.0e-03 / 1.0e-03,
            "patch_icpp(1)%z_centroid": 80.0e-03 / 1.0e-03,
            "patch_icpp(1)%length_x": 160.0e-03 / 1.0e-03,
            "patch_icpp(1)%length_y": 160.0e-03 / 1.0e-03,
            "patch_icpp(1)%length_z": 160.0e-03 / 1.0e-03,
            # Specify the patch primitive variables
            "patch_icpp(1)%vel(1)": 0.0e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%vel(3)": 0.0e00,
            "patch_icpp(1)%pres": 1.0e05,
            "patch_icpp(1)%alpha_rho(1)": 1000.0e00,
            "patch_icpp(1)%alpha_rho(2)": 0.1e00,
            "patch_icpp(1)%alpha(1)": 0.9e00,
            "patch_icpp(1)%alpha(2)": 0.1e00,
            # Patch 3: Air bubble
            # Specify the spherical gas bubble grid geometry
            "patch_icpp(2)%geometry": 8,
            "patch_icpp(2)%smoothen": "T",
            "patch_icpp(2)%smooth_patch_id": 1,
            "patch_icpp(2)%smooth_coeff": 0.5e00,
            "patch_icpp(2)%x_centroid": 0.0e00,
            "patch_icpp(2)%y_centroid": 0.0e00,
            "patch_icpp(2)%z_centroid": 0.0e00,
            "patch_icpp(2)%radius": 1.0e-03 / 1.0e-03,
            "patch_icpp(2)%alter_patch(1)": "T",
            # Specify the patch primitive variables
            "patch_icpp(2)%vel(1)": 0.0e00,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%vel(3)": 0.0e00,
            "patch_icpp(2)%pres": 1.0e03,
            "patch_icpp(2)%alpha_rho(1)": 100.0e00,
            "patch_icpp(2)%alpha_rho(2)": 0.9e00,
            "patch_icpp(2)%alpha(1)": 0.1e00,
            "patch_icpp(2)%alpha(2)": 0.9e00,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 4.4e00 * 6.0e08 / (4.4e00 - 1.0e00),
            "fluid_pp(2)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(2)%pi_inf": 0.0e00,
        }
    )
)
