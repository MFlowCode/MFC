#!/usr/bin/env python3
"""
2D JWL+Air Blast Tube with Moving Immersed Cylinder

This test case exercises mixed-material JWL EOS with a moving 2D immersed
cylinder. A compact TNT-like JWL charge expands into air and forms a shock
that travels down a longer channel before loading the cylinder.

Configuration notes:
- Longer domain gives the blast wave room to develop before hitting the IB.
- Smaller charge radius avoids an unrealistically huge initial product region.
- Cylinder is a 2D circular IB, which represents a long cylinder cross-section.
- Grid spacing gives about 20 cells per cylinder diameter, or 10 per radius.
"""

import json

eps = 1.0e-8  # Interface thickness parameter

# Fluid densities at reference conditions
rho_jwl = 1630.0  # Detonation products (TNT-like)
rho_air = 1.225  # Ambient air at sea level

# JWL equation of state parameters (TNT-like explosive)
jwl = {
    "A": 3.712e11,  # JWL A coefficient (Pa)
    "B": 3.231e9,  # JWL B coefficient (Pa)
    "R1": 4.15,  # JWL R1 exponent
    "R2": 0.95,  # JWL R2 exponent
    "omega": 0.30,  # Gruneisen parameter
    "rho0": rho_jwl,  # Reference density (kg/m³)
    "E0": 1.0089e10,  # Reference internal energy (J/kg)
    "air_e0": 2.5575e5,  # Air internal energy (J/kg)
    "air_rho0": rho_air,  # Air reference density
    "air_gamma": 0.4,  # Air heat capacity ratio coefficient
}

print(
    json.dumps(
        {
            "run_time_info": "T",
            # ========== Domain Configuration ==========
            "x_domain%beg": 0.0,
            "x_domain%end": 2.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 0.5,
            "m": 250,  # dx = 0.004, about 20 cells per cylinder diameter
            "n": 75,  # dy = 0.004, matching x-resolution
            "p": 0,  # 2D case
            # ========== Time Stepping with CFL Adaptation ==========
            # "dt": 1.0e-8,    # Initial small dt for safety with high gradients
            "t_step_start": 0,
            "t_step_stop": 4000,
            "t_step_save": 50,
            "cfl_adap_dt": "T",  # Enable adaptive CFL stepping
            "cfl_target": 0.35,  # Conservative for compact high-pressure blast
            "n_start": 0,
            "t_stop": 2.0e-4,  # Long enough for shock-cylinder interaction
            "t_save": 5.0e-6,  # Frequent output for blast/IB diagnostics
            # ========== Numerics & Discretization ==========
            "num_patches": 2,
            "model_eqns": 2,  # Multicomponent flow
            "num_fluids": 2,  # JWL + air
            "mpp_lim": "T",  # Multi-phase pressure limiter
            "mixture_err": "T",  # Mixture error tracking
            "alt_soundspeed": "F",
            "time_stepper": 3,  # 3rd order TVD RK
            "weno_order": 3,  # 3rd order WENO-Z (recommended for 2D)
            "weno_eps": 1.0e-16,  # WENO epsilon for robustness
            "mapped_weno": "T",  # Mapped WENO for better conservation
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,  # HLLC (better for multiphase)
            "wave_speeds": 1,  # Roe state for wave speeds
            "avg_state": 2,  # Mixture-based averaging
            "bc_x%beg": -3,  # Non-reflective BC (transmissive)
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # ========== Immersed Boundary Setup ==========
            "ib": "T",
            "num_ibs": 1,
            "viscous": "T",  # Include viscous effects near IB
            "weno_Re_flux": "T",  # WENO-based viscous fluxes
            "weno_avg": "T",  # WENO averaging at IB interface
            # ========== Output Configuration ==========
            "format": 1,  # ASCII format
            "precision": 2,  # Double precision
            "prim_vars_wrt": "T",
            "rho_wrt": "T",
            "pres_wrt": "T",  # Pressure important for blast diagnostics
            "c_wrt": "T",  # Sound speed for shock identification
            "E_wrt": "T",  # Total energy
            "ib_state_wrt": "T",  # IB response tracking
            "parallel_io": "F",
            # ========== Patch 1: Air Background ==========
            "patch_icpp(1)%geometry": 3,  # Rectangular patch
            "patch_icpp(1)%x_centroid": 1.0,
            "patch_icpp(1)%y_centroid": 0.25,
            "patch_icpp(1)%length_x": 2.0,
            "patch_icpp(1)%length_y": 0.5,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 101325.0,  # Standard atmospheric pressure
            "patch_icpp(1)%alpha_rho(1)": eps * rho_jwl,
            "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_air,
            "patch_icpp(1)%alpha(1)": eps,
            "patch_icpp(1)%alpha(2)": 1.0 - eps,
            # ========== Patch 2: JWL High-Pressure Driver ==========
            "patch_icpp(2)%geometry": 2,  # Circular patch
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.16,
            "patch_icpp(2)%y_centroid": 0.25,
            "patch_icpp(2)%radius": 0.035,  # Compact charge region
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 2.0e10,  # Strong compact JWL blast driver
            "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_jwl,
            "patch_icpp(2)%alpha_rho(2)": eps * rho_air,
            "patch_icpp(2)%alpha(1)": 1.0 - eps,
            "patch_icpp(2)%alpha(2)": eps,
            # ========== Immersed Boundary (Moving 2D Cylinder Cross-Section) ==========
            "patch_ib(1)%geometry": 2,  # Circular IB = 2D cylinder cross-section
            "patch_ib(1)%x_centroid": 0.50,  # Downstream of compact charge
            "patch_ib(1)%y_centroid": 0.25,
            "patch_ib(1)%radius": 0.04,  # Cylinder radius
            "patch_ib(1)%slip": "T",  # Slippage allowed (frictionless)
            "patch_ib(1)%moving_ibm": 2,  # Two-way coupled dynamics
            "patch_ib(1)%vel(1)": 0.0,
            "patch_ib(1)%vel(2)": 0.0,
            "patch_ib(1)%vel(3)": 0.0,
            "patch_ib(1)%angles(1)": 0.0,
            "patch_ib(1)%angles(2)": 0.0,
            "patch_ib(1)%angles(3)": 0.0,
            "patch_ib(1)%angular_vel(1)": 0.0,
            "patch_ib(1)%angular_vel(2)": 0.0,
            "patch_ib(1)%angular_vel(3)": 0.0,
            "patch_ib(1)%mass": 0.03,  # Light cylinder for visible blast response
            # ========== Fluid 1: JWL Detonation Products ==========
            "fluid_pp(1)%eos": 2,  # JWL EOS
            "fluid_pp(1)%gamma": 1.0 / 0.4,  # Heat capacity ratio for products
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%cv": 613.5,  # Specific heat at constant volume
            "fluid_pp(1)%Re(1)": 2.5e6,  # High Reynolds (inviscid dominant)
            "fluid_pp(1)%jwl_A": jwl["A"],
            "fluid_pp(1)%jwl_B": jwl["B"],
            "fluid_pp(1)%jwl_R1": jwl["R1"],
            "fluid_pp(1)%jwl_R2": jwl["R2"],
            "fluid_pp(1)%jwl_omega": jwl["omega"],
            "fluid_pp(1)%jwl_rho0": jwl["rho0"],
            "fluid_pp(1)%jwl_E0": jwl["E0"],
            "fluid_pp(1)%jwl_air_e0": jwl["air_e0"],
            "fluid_pp(1)%jwl_air_rho0": jwl["air_rho0"],
            "fluid_pp(1)%jwl_air_gamma": jwl["air_gamma"],
            # ========== Fluid 2: Ideal Gas Air ==========
            "fluid_pp(2)%eos": 1,  # Ideal gas EOS
            "fluid_pp(2)%gamma": 1.0 / 0.4,  # Air specific heat ratio
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%Re(1)": 2.5e6,  # High Reynolds
        }
    )
)
