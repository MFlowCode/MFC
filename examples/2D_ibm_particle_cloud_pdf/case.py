import json
import math

# 2D Mach 1.5 shock passing a bed of free immersed-boundary particles placed with an analytic pdf: the placement density grows
# as the square of height, as if the particles were falling from the top of the bed when the shock arrives.

heat_capacity_ratio = 1.4

# Mach 1.5 shock, pre-shock sound speed 1, Anderson (2003) Modern Compressible Flow ch. 7
pre_shock_pressure = 1.0
pre_shock_density = 1.4
pre_shock_speed = 0.0
post_shock_pressure = 2.4583
post_shock_density = 2.6069
post_shock_speed = 0.6944

# Particles: 10 cells per diameter, ten times denser than the gas
domain_size = 16.0
particle_diameter = domain_size / 150.0
particle_density = 10.0 * pre_shock_density
particle_reynolds_number = 20.0  # on post-shock density and speed
viscosity = post_shock_density * post_shock_speed * particle_diameter / particle_reynolds_number

# Bed x in [10D, 50D] spanning the channel one diameter off each wall, filled to a 10% area fraction
bed_length_x = 40.0 * particle_diameter
bed_length_y = domain_size - 2.0 * particle_diameter
bed_x_centroid = 30.0 * particle_diameter
volume_fraction = 0.1
num_particles = round(volume_fraction * bed_length_x * domain_size / (math.pi * (0.5 * particle_diameter) ** 2))

# Grid and time step, CFL 0.33 on post-shock u + c
num_cells = round(10.0 * domain_size / particle_diameter)
cell_size = domain_size / num_cells
wave_front = 5.0 * particle_diameter
post_shock_sound_speed = math.sqrt(heat_capacity_ratio * post_shock_pressure / post_shock_density)
time_step = 0.33 * cell_size / (post_shock_speed + post_shock_sound_speed)
num_time_steps = 9800
steps_to_save = 100

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": domain_size,
            "y_domain%beg": -0.5 * domain_size,
            "y_domain%end": 0.5 * domain_size,
            "m": num_cells - 1,
            "n": num_cells - 1,
            "p": 0,
            "dt": time_step,
            "t_step_start": 0,
            "t_step_stop": num_time_steps,
            "t_step_save": steps_to_save,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": "5eq",
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": "rk3",
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "T",
            "weno_avg": "T",
            "avg_state": "arithmetic",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": "hll",
            "wave_speeds": "direct",
            "fd_order": 2,
            "viscous": "T",
            "bc_x%beg": -17,
            "bc_x%end": -8,
            "bc_y%beg": -15,
            "bc_y%end": -15,
            # Immersed boundaries: every particle comes from the particle cloud
            "ib": "T",
            "num_ibs": 0,
            "collision_model": 1,
            "coefficient_of_restitution": 0.9,
            "collision_time": 20.0 * time_step,
            "ib_coefficient_of_friction": 0.1,
            # Particle cloud: box bed, pdf zero at its bottom edge and largest at its top (compiled into pre_process)
            "num_particle_clouds": 1,
            "particle_cloud(1)%particle_type": "IB",
            "particle_cloud(1)%cloud_geometry": "box",
            "particle_cloud(1)%packing_method": 1,
            "particle_cloud(1)%seed": 42,
            "particle_cloud(1)%x_centroid": bed_x_centroid,
            "particle_cloud(1)%y_centroid": 0.0,
            "particle_cloud(1)%z_centroid": 0.0,
            "particle_cloud(1)%length_x": bed_length_x,
            "particle_cloud(1)%length_y": bed_length_y,
            "particle_cloud(1)%num_particles": num_particles,
            "particle_cloud(1)%radius": 0.5 * particle_diameter,
            "particle_cloud(1)%min_spacing": 0.2 * particle_diameter,
            "particle_cloud(1)%pdf": "((y - (yc - 0.5*ly)) / ly)**2",
            "particle_cloud(1)%density": particle_density,
            "particle_cloud(1)%moving_ibm": 2,
            # Output
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "ib_state_wrt": "T",
            "parallel_io": "T",
            # IC Patch 1: post-shock region (left of wave front)
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * wave_front,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": wave_front,
            "patch_icpp(1)%length_y": domain_size,
            "patch_icpp(1)%vel(1)": post_shock_speed,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": post_shock_pressure,
            "patch_icpp(1)%alpha_rho(1)": post_shock_density,
            "patch_icpp(1)%alpha(1)": 1.0,
            # IC Patch 2: pre-shock region (right of wave front)
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.5 * (wave_front + domain_size),
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%length_x": domain_size - wave_front,
            "patch_icpp(2)%length_y": domain_size,
            "patch_icpp(2)%vel(1)": pre_shock_speed,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": pre_shock_pressure,
            "patch_icpp(2)%alpha_rho(1)": pre_shock_density,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluid properties: air
            "fluid_pp(1)%gamma": 1.0 / (heat_capacity_ratio - 1.0),
            "fluid_pp(1)%eos": "ideal_gas",
            "fluid_pp(1)%Re(1)": 1.0 / viscosity,
        }
    )
)
