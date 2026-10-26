import json
import math

# 3D Mach 1.5 shock passing a rotated cuboid cloud of rigid Euler-Lagrange particles with lognormal radii, placed by
# rejection sampling inside the cuboid patch shape.

gam_a = 1.4

# Shock parameters (Mach 1.5, pre-shock sound speed 1)
pre_shock_pressure = 1.0
pre_shock_density = 1.4
pre_shock_speed = 0.0
post_shock_pressure = 2.4583
post_shock_density = 2.6069
post_shock_speed = 0.6944

# Particles: median a fifth of a cell across, ten times denser than the gas
domain_size = 16.0
particle_diameter = domain_size / 150.0
particle_radius = 0.5 * particle_diameter  # median radius
radius_sigma = 0.3  # standard deviation of ln(radius / particle_radius), clipped at +/- 3 sigma
particle_density = 10.0 * pre_shock_density
particle_reynolds_number = 20.0  # on post-shock density and speed
viscosity = post_shock_density * post_shock_speed * particle_diameter / particle_reynolds_number

# Cuboid cloud 6 x 3 x 3, rotated about each axis, filled to a 1% volume fraction
cloud_lengths = (6.0, 3.0, 3.0)
cloud_angles = (0.3, 0.4, 0.5)  # radians about x, y, z
cloud_x = 0.4 * domain_size
volume_fraction = 0.01
mean_particle_volume = 4.0 / 3.0 * math.pi * particle_radius**3 * math.exp(4.5 * radius_sigma**2)  # lognormal third moment
num_particles = round(volume_fraction * math.prod(cloud_lengths) / mean_particle_volume)

# Grid and time step, CFL 0.33 on post-shock u + c
num_cells = round(0.2 * domain_size / particle_diameter)
dx = domain_size / num_cells
wave_front = 5.0 * particle_diameter
post_shock_sound_speed = math.sqrt(gam_a * post_shock_pressure / post_shock_density)
dt = 0.33 * dx / (post_shock_speed + post_shock_sound_speed)
num_time_steps = 400
steps_to_save = 20

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
            "z_domain%beg": -0.5 * domain_size,
            "z_domain%end": 0.5 * domain_size,
            "m": num_cells - 1,
            "n": num_cells - 1,
            "p": num_cells - 1,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": num_time_steps,
            "t_step_save": steps_to_save,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": "5eq",
            "alt_soundspeed": "F",
            "num_fluids": 2,
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
            "bc_z%beg": -15,
            "bc_z%end": -15,
            # Euler-Lagrange particles: rigid (fixed radius), inertial, Stokes drag
            "bubbles_lagrange": "T",
            "bubble_model": 0,
            "thermal": 3,
            "polytropic": "F",
            "lag_params%nBubs_glb": num_particles,
            "lag_params%vel_model": 2,
            "lag_params%drag_model": 2,
            "lag_params%pressure_force": "T",
            "lag_params%solver_approach": 2,
            "lag_params%cluster_type": 1,
            "lag_params%pressure_corrector": "F",
            "lag_params%smooth_type": 1,
            "lag_params%heatTransfer_model": "F",
            "lag_params%massTransfer_model": "F",
            "lag_params%epsilonb": 1.0,
            "lag_params%valmaxvoid": 0.9,
            "lag_params%write_bubbles": "F",
            # Particle cloud: rotated cuboid patch shape, lognormal radii, each mass from density
            "num_particle_clouds": 1,
            "particle_cloud(1)%particle_type": "EL",
            "particle_cloud(1)%cloud_geometry": "cuboid",
            "particle_cloud(1)%packing_method": 1,
            "particle_cloud(1)%seed": 42,
            "particle_cloud(1)%x_centroid": cloud_x,
            "particle_cloud(1)%y_centroid": 0.0,
            "particle_cloud(1)%z_centroid": 0.0,
            "particle_cloud(1)%length_x": cloud_lengths[0],
            "particle_cloud(1)%length_y": cloud_lengths[1],
            "particle_cloud(1)%length_z": cloud_lengths[2],
            "particle_cloud(1)%angles(1)": cloud_angles[0],
            "particle_cloud(1)%angles(2)": cloud_angles[1],
            "particle_cloud(1)%angles(3)": cloud_angles[2],
            "particle_cloud(1)%num_particles": num_particles,
            "particle_cloud(1)%radius": particle_radius,
            "particle_cloud(1)%radius_sigma": radius_sigma,
            "particle_cloud(1)%min_spacing": 0.2 * particle_diameter,
            "particle_cloud(1)%density": particle_density,
            # Output
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "lag_db_wrt": "T",
            "lag_pos_wrt": "T",
            "lag_vel_wrt": "T",
            "parallel_io": "T",
            # IC Patch 1: post-shock region (left of wave front)
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.5 * wave_front,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_x": wave_front,
            "patch_icpp(1)%length_y": domain_size,
            "patch_icpp(1)%length_z": domain_size,
            "patch_icpp(1)%vel(1)": post_shock_speed,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": post_shock_pressure,
            "patch_icpp(1)%alpha_rho(1)": post_shock_density,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%alpha_rho(2)": 0.0,
            "patch_icpp(1)%alpha(2)": 0.0,
            # IC Patch 2: pre-shock region (right of wave front)
            "patch_icpp(2)%geometry": 9,
            "patch_icpp(2)%x_centroid": 0.5 * (wave_front + domain_size),
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%z_centroid": 0.0,
            "patch_icpp(2)%length_x": domain_size - wave_front,
            "patch_icpp(2)%length_y": domain_size,
            "patch_icpp(2)%length_z": domain_size,
            "patch_icpp(2)%vel(1)": pre_shock_speed,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": pre_shock_pressure,
            "patch_icpp(2)%alpha_rho(1)": pre_shock_density,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%alpha_rho(2)": 0.0,
            "patch_icpp(2)%alpha(2)": 0.0,
            # Fluid properties: air carrier and the particle gas the Euler-Lagrange model requires (never on the grid)
            "fluid_pp(1)%gamma": 1.0 / (gam_a - 1.0),
            "fluid_pp(1)%eos": "ideal_gas",
            "fluid_pp(1)%Re(1)": 1.0 / viscosity,
            "fluid_pp(2)%gamma": 1.0 / (gam_a - 1.0),
            "fluid_pp(2)%eos": "ideal_gas",
            "fluid_pp(2)%Re(1)": 1.0 / viscosity,
            # Particle gas properties: nondimensional reference states
            "bub_pp%R0ref": 1.0,
            "bub_pp%p0ref": 1.0,
            "bub_pp%rho0ref": 1.0,
            "bub_pp%T0ref": 1.0,
            "bub_pp%ss": 1.0e-12,
            "bub_pp%pv": 1.0e-6,
            "bub_pp%vd": 1.0e-6,
            "bub_pp%mu_l": viscosity,
            "bub_pp%gam_v": 1.33,
            "bub_pp%gam_g": gam_a,
            "bub_pp%M_v": 18.0,
            "bub_pp%M_g": 28.97,
            "bub_pp%k_v": 1.0e-3,
            "bub_pp%k_g": 1.0e-3,
            "bub_pp%cp_v": 1.0,
            "bub_pp%cp_g": 1.0,
            "bub_pp%R_v": 1.0 / particle_density,
            "bub_pp%R_g": pre_shock_pressure / particle_density,
        }
    )
)
