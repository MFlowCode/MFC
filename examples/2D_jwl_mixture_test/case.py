#!/usr/bin/env python3
"""
2D JWL+Air Bottom-Wall Blast with Moving Particle Bed

This test case exercises mixed-material JWL EOS with a moving 2D immersed
particle bed. A compact TNT-like JWL charge is placed near the lower slip
wall and drives a shock into a semicircular bed of small cylindrical particles.

Configuration notes:
- Bottom boundary is a slip wall and the other open boundaries use -3.
- Particles are 2D circular IBs, representing cylinder cross-sections.
- Particle bed is arranged in an upper half-cylinder around the blast.
"""

import json
import math

eps = 1.0e-8  # Volume-fraction floor to avoid exactly pure-fluid mixture states

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

params = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": 2.0,
    "y_domain%beg": 0.0,
    "y_domain%end": 1.2,
    "m": 200,
    "n": 120,
    "p": 0,
    "t_step_start": 0,
    "t_step_stop": 5000,
    "t_step_save": 50,
    "cfl_adap_dt": "T",
    "cfl_target": 0.35,
    "n_start": 0,
    "t_stop": 2.0e-4,
    "t_save": 5.0e-6,
    "num_patches": 2,
    "model_eqns": 2,
    "num_fluids": 2,
    "mpp_lim": "T",
    "mixture_err": "T",
    "alt_soundspeed": "F",
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
    "bc_y%beg": -1,
    "bc_y%end": -3,
    "ib": "T",
    "num_ibs": 40,
    "viscous": "T",
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "collision_model": 1,
    "ib_coefficient_of_friction": 0.03,
    "collision_time": 5.0e-7,
    "coefficient_of_restitution": 0.8,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "rho_wrt": "T",
    "pres_wrt": "T",
    "c_wrt": "T",
    "E_wrt": "T",
    "ib_state_wrt": "T",
    "parallel_io": "F",
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 1.0,
    "patch_icpp(1)%y_centroid": 0.6,
    "patch_icpp(1)%length_x": 2.0,
    "patch_icpp(1)%length_y": 1.2,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 101325.0,
    "patch_icpp(1)%alpha_rho(1)": eps * rho_jwl,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_air,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": 1.0,
    "patch_icpp(2)%y_centroid": 0.06,
    "patch_icpp(2)%radius": 0.025,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": 2.0e10,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_jwl,
    "patch_icpp(2)%alpha_rho(2)": eps * rho_air,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "fluid_pp(1)%eos": 2,
    "fluid_pp(1)%gamma": 1.0 / 0.4,
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%cv": 613.5,
    "fluid_pp(1)%Re(1)": 2.5e6,
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
    "fluid_pp(2)%eos": 1,
    "fluid_pp(2)%gamma": 1.0 / 0.4,
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(2)%Re(1)": 2.5e6,
}

particle_radius = 0.01
particle_mass = 1.5e-3
bed_center = (1.0, 0.06)
particle_diameter = 2.0 * particle_radius
center_spacing = 1.2 * particle_diameter
bed_radii = [0.14 + i * center_spacing for i in range(5)]
bed_angles = [35.0, 50.714, 66.429, 82.143, 97.857, 113.571, 129.286, 145.0]

patch_id = 0
for radius in bed_radii:
    for angle_deg in bed_angles:
        patch_id += 1
        theta = math.radians(angle_deg)
        x_centroid = bed_center[0] + radius * math.cos(theta)
        y_centroid = bed_center[1] + radius * math.sin(theta)

        params[f"patch_ib({patch_id})%geometry"] = 2
        params[f"patch_ib({patch_id})%x_centroid"] = x_centroid
        params[f"patch_ib({patch_id})%y_centroid"] = y_centroid
        params[f"patch_ib({patch_id})%radius"] = particle_radius
        params[f"patch_ib({patch_id})%slip"] = "T"
        params[f"patch_ib({patch_id})%moving_ibm"] = 2
        params[f"patch_ib({patch_id})%vel(1)"] = 0.0
        params[f"patch_ib({patch_id})%vel(2)"] = 0.0
        params[f"patch_ib({patch_id})%vel(3)"] = 0.0
        params[f"patch_ib({patch_id})%angles(1)"] = 0.0
        params[f"patch_ib({patch_id})%angles(2)"] = 0.0
        params[f"patch_ib({patch_id})%angles(3)"] = 0.0
        params[f"patch_ib({patch_id})%angular_vel(1)"] = 0.0
        params[f"patch_ib({patch_id})%angular_vel(2)"] = 0.0
        params[f"patch_ib({patch_id})%angular_vel(3)"] = 0.0
        params[f"patch_ib({patch_id})%mass"] = particle_mass

print(json.dumps(params))
