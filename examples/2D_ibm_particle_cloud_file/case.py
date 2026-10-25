import json

# 2D shock wave passing a static bed of circular IBs whose centres are read from centres.txt.

gam_a = 1.4

# Shock parameters (Mach 1.5)
pre_shock_pressure = 1
pre_shock_density = 1.4
post_shock_pressure = 2.4583
post_shock_density = 2.6069
post_shock_speed = 0.6944

domain_size = 4.0
wave_front = -1.5

total_time = 1.0
num_time_steps = 1000
dt = float(total_time / num_time_steps)
num_saves = 50
steps_to_save = int(num_time_steps / num_saves)

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
            "m": 200,
            "n": 200,
            "p": 0,
            "dt": dt,
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
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "avg_state": "arithmetic",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            "riemann_solver": "hllc",
            "wave_speeds": "direct",
            "fd_order": 2,
            "bc_x%beg": -17,
            "bc_x%end": -8,
            "bc_y%beg": -15,
            "bc_y%end": -15,
            # Immersed boundaries: every circle comes from the file-based particle cloud
            "ib": "T",
            "num_ibs": 0,
            "num_particle_clouds": 1,
            "particle_cloud(1)%particle_type": "IB",
            "particle_cloud(1)%cloud_geometry": 21,
            "particle_cloud(1)%filepath": "centres.txt",
            "particle_cloud(1)%num_particles": 6,
            "particle_cloud(1)%radius": 0.2,
            "particle_cloud(1)%mass": 1.0,
            # Output
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "ib_state_wrt": "T",
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
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": pre_shock_pressure,
            "patch_icpp(2)%alpha_rho(1)": pre_shock_density,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluid properties: air
            "fluid_pp(1)%gamma": 1.0 / (gam_a - 1.0),
            "fluid_pp(1)%eos": "ideal_gas",
        }
    )
)
