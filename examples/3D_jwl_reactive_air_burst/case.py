#!/usr/bin/env python3
# 3D spherical JWL++ reactive detonation bursting into air.
#
# A full-density TNT sphere (rho0 = 1630, radius 5 cm) sits at rest at ambient
# pressure in free air. A small central hot spot ignites a self-propagating
# JWL++ (jwl_reactive, Souers 2000) detonation, dl/dt = G p^b (1-l), that sweeps
# outward at the TNT CJ speed and then vents the products into the surrounding
# air -- the 3D analogue of examples/2D_jwl_reactive_air_mixing. The Y*(1-lambda)
# gate keeps only the unreacted explosive reacting; the air stays inert.
#
# Sized to run in ~8 min on 8 CPU cores (run with --mpi -n 8): 120^3 mesh (dx=1.7mm),
# 400 steps. A fast JWL++ rate (jwl_G) keeps the reaction zone sub-cell so the detonation
# self-sustains at the TNT CJ speed on this coarse (laptop-affordable) mesh.
import json

from mfc.jwl_products import AIR, TNT, ambient_fluid, jwl_fluid, znd_delta_e

# TNT products in air, from the shared preset library; jwl_rho0/air_rho feed the patch ICs.
jwl_rho0 = TNT["rho0"]
air_rho = AIR["rho0"]
p_amb = 101325.0

# Reactant/product energy offset (Garno Eq. 17): unreacted TNT at (rho0, e=0) at ambient p.
jwl_delta_e = znd_delta_e(TNT, p_amb)

# JWL++ rate: dl/dt = jwl_G * p^jwl_b_exp * (1 - lambda)
jwl_G = 5.0e-14
jwl_b_exp = 2.0

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
            "m": 119,
            "n": 119,
            "p": 119,
            "dt": 3.0e-8,
            "t_step_start": 0,
            "t_step_stop": 400,
            "t_step_save": 20,
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
            "bc_z%beg": -3,
            "bc_z%end": -3,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Self-propagating JWL++ reactive burn.
            "jwl_reactive": "T",
            "jwl_G": jwl_G,
            "jwl_b_exp": jwl_b_exp,
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
            "patch_icpp(1)%pres": p_amb,
            "patch_icpp(1)%alpha_rho(1)": jwl_rho0 * 1.0e-8,
            "patch_icpp(1)%alpha_rho(2)": air_rho * (1.0 - 1.0e-8),
            "patch_icpp(1)%alpha(1)": 1.0e-8,
            "patch_icpp(1)%alpha(2)": 0.99999999,
            # Patch 2: TNT charge sphere, unreacted, at ambient pressure.
            "patch_icpp(2)%geometry": 8,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%z_centroid": 0.0,
            "patch_icpp(2)%radius": 0.05,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": p_amb,
            "patch_icpp(2)%alpha_rho(1)": jwl_rho0 * (1.0 - 1.0e-8),
            "patch_icpp(2)%alpha_rho(2)": air_rho * 1.0e-8,
            "patch_icpp(2)%alpha(1)": 0.99999999,
            "patch_icpp(2)%alpha(2)": 1.0e-8,
            # Patch 3: central hot spot (still TNT) igniting the burn.
            "patch_icpp(3)%geometry": 8,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%alter_patch(2)": "T",
            "patch_icpp(3)%x_centroid": 0.0,
            "patch_icpp(3)%y_centroid": 0.0,
            "patch_icpp(3)%z_centroid": 0.0,
            "patch_icpp(3)%radius": 0.006,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%vel(3)": 0.0,
            "patch_icpp(3)%pres": 25.0e9,
            "patch_icpp(3)%alpha_rho(1)": jwl_rho0 * (1.0 - 1.0e-8),
            "patch_icpp(3)%alpha_rho(2)": air_rho * 1.0e-8,
            "patch_icpp(3)%alpha(1)": 0.99999999,
            "patch_icpp(3)%alpha(2)": 1.0e-8,
            # Fluid 1: TNT JWL products (+ ambient references); Fluid 2: air.
            **jwl_fluid(1, TNT, AIR, delta_e=jwl_delta_e),
            **ambient_fluid(2, AIR),
        }
    )
)
