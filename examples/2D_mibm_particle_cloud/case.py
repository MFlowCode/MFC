import json
import math

# 2D shock wave interacting with a bed of 20 free-floating circular particles.

gam_a = 1.4

# Shock parameters (Mach 1.5)
mach_number = 1.5
pre_shock_pressure = 1
pre_shock_density = 1.4
pre_shock_speed = 0.0
post_shock_pressure = 2.4583
post_shock_density = 2.6069
post_shock_speed = 0.6944

domain_size = 4.0
wave_front = -1.5

total_time = 1.5
num_time_steps = 2000
dt = float(total_time / num_time_steps)
num_saves = 100
steps_to_save = int(num_time_steps / num_saves)

# Soft-sphere collision parameters (from 3D_mibm_sphere_head_on_collision)
collision_time = 20.0 * dt

# Particle bed parameters
bed_x = 0.5
bed_y = 0.0
bed_lx = 2.0
bed_ly = 3.5
particle_radius = 0.15
particle_mass = 0.25
particle_min_spacing = 0.05

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -domain_size * 0.5,
            "x_domain%end": domain_size * 0.5,
            "y_domain%beg": -domain_size * 0.5,
            "y_domain%end": domain_size * 0.5,
            "cyl_coord": "F",
            "m": 256,
            "n": 256,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": num_time_steps,
            "t_step_save": steps_to_save,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "T",
            "weno_avg": "T",
            "avg_state": 2,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "bc_x%beg": -17,
            "bc_x%end": -8,
            "bc_y%beg": -15,
            "bc_y%end": -15,
            # Immersed boundaries — all circles come from the particle bed
            "ib": "T",
            "num_ibs": 0,
            "viscous": "T",
            "many_ib_patch_parallelism": "T",
            # Collision model (soft-sphere, from 3D_mibm_sphere_head_on_collision)
            "collision_model": 1,
            "coefficient_of_restitution": 0.9,
            "collision_time": collision_time,
            "ib_coefficient_of_friction": 0.1,
            # Particle bed: 20 free-floating circles placed randomly in region
            "num_particle_clouds": 1,
            "particle_cloud(1)%x_centroid": bed_x,
            "particle_cloud(1)%y_centroid": bed_y,
            "particle_cloud(1)%z_centroid": 0.0,
            "particle_cloud(1)%length_x": bed_lx,
            "particle_cloud(1)%length_y": bed_ly,
            "particle_cloud(1)%length_z": 0.0,
            "particle_cloud(1)%num_particles": 20,
            "particle_cloud(1)%radius": particle_radius,
            "particle_cloud(1)%mass": particle_mass,
            "particle_cloud(1)%min_spacing": particle_min_spacing,
            "particle_cloud(1)%moving_ibm": 2,
            "particle_cloud(1)%seed": 42,
            "particle_cloud(1)%packing_method": 1,
            # Output
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "E_wrt": "T",
            "ib_state_wrt": "F",
            "parallel_io": "T",
            # IC Patch 1: post-shock region (left of wave front)
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * wave_front - 0.25 * domain_size,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 0.5 * domain_size + wave_front,
            "patch_icpp(1)%length_y": domain_size,
            "patch_icpp(1)%vel(1)": post_shock_speed,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": post_shock_pressure,
            "patch_icpp(1)%alpha_rho(1)": post_shock_density,
            "patch_icpp(1)%alpha(1)": 1.0,
            # IC Patch 2: pre-shock region (right of wave front)
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.5 * wave_front + 0.25 * domain_size,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%length_x": 0.5 * domain_size - wave_front,
            "patch_icpp(2)%length_y": domain_size,
            "patch_icpp(2)%vel(1)": pre_shock_speed,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": pre_shock_pressure,
            "patch_icpp(2)%alpha_rho(1)": pre_shock_density,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluid properties: air
            "fluid_pp(1)%gamma": 1.0 / (gam_a - 1.0),
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(1)%Re(1)": 2500000,
        }
    )
)
