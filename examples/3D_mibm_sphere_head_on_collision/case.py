import json
import math

Mu = 1.84e-05
gam_a = 1.4

# lead-up-properties
lead_distance = 0.1
velocity = 1.0
dt = 0.5e-3
collision_time = 20. * dt

# simulation runs long enough to collide and travel about lead distance away again
simulation_time = 2. * (lead_distance / velocity) + collision_time
num_time_steps = int(simulation_time / dt)
num_saves = 10
t_step_save = int(num_time_steps / num_saves)

# parerticle properties
radius = 1.0
collision_angle_degrees = 45.
collision_anlge_radians = collision_angle_degrees * math.pi / 180.
collision_point_y = radius * math.sin(collision_anlge_radians)
collision_point_x = radius * math.cos(collision_anlge_radians)
domain_size = 5. * radius

# Configuring case dictionary
print(json.dumps({
      # Logistics
        "run_time_info": "T",
        # Computational Domain Parameters
        "x_domain%beg": -0.5 * domain_size,
        "x_domain%end": 0.5 * domain_size,
        "y_domain%beg": -0.5 * domain_size,
        "y_domain%end": 0.5 * domain_size,
        "z_domain%beg": -0.5 * domain_size,
        "z_domain%end": 0.5 * domain_size,
        "cyl_coord": "F",
        "m": 100,
        "n": 100,
        "p": 100,
        "dt": dt,
        "t_step_start": 0,
        "t_step_stop": num_time_steps,
        "t_step_save": t_step_save,
        # Simulation Algorithm Parameters
        "num_patches": 1,
        # Use the 5 equation model
        "model_eqns": 2,
        "alt_soundspeed": "F",
        # One fluids: air
        "num_fluids": 1,
        "mpp_lim": "F",
        # Correct errors when computing speed of sound
        "mixture_err": "T",
        # Use TVD RK3 for time marching
        "time_stepper": 3,
        # Use WENO5
        "weno_order": 5,
        "weno_eps": 1.0e-16,
        "weno_avg": "T",
        "avg_state": 2,
        "mapped_weno": "T",
        "null_weights": "F",
        "mp_weno": "T",
        "riemann_solver": 2,
        "wave_speeds": 1,
        # We use ghost-cell
        "bc_x%beg": -3,
        "bc_x%end": -3,
        "bc_y%beg": -3,
        "bc_y%end": -3,
        "bc_z%beg": -3,
        "bc_z%end": -3,
        # Set IB to True and add 1 patch
        "ib": "T",
        "num_ibs": 2,
        # Formatted Database Files Structure Parameters
        "format": 1,
        "precision": 2,
        "prim_vars_wrt": "T",
        "E_wrt": "T",
        "parallel_io": "T",
        # Patch: Constant Tube filled with air
        # Specify the cylindrical air tube grid geometry
        "patch_icpp(1)%geometry": 9,
        "patch_icpp(1)%x_centroid": 0.0,
        "patch_icpp(1)%y_centroid": 0.0,
        "patch_icpp(1)%z_centroid": 0.0,
        "patch_icpp(1)%length_x": 0.5 * domain_size,
        "patch_icpp(1)%length_y": 0.5 * domain_size,
        "patch_icpp(1)%length_z": 0.5 * domain_size,
        # Specify the patch primitive variables
        "patch_icpp(1)%vel(1)": 0.0e00,
        "patch_icpp(1)%vel(2)": 0.0e00,
        "patch_icpp(1)%vel(3)": 0.0e00,
        "patch_icpp(1)%pres": 1.0e00,
        "patch_icpp(1)%alpha_rho(1)": 1.0e00,
        "patch_icpp(1)%alpha(1)": 1.0e00,
        # Patch: Cylinder Immersed Boundary
        "patch_ib(1)%geometry": 8,
        "patch_ib(1)%x_centroid": 2. * collision_point_x + lead_distance,
        "patch_ib(1)%y_centroid": 2. * collision_point_y,
        "patch_ib(1)%z_centroid": 0.0,
        "patch_ib(1)%radius": radius,
        "patch_ib(1)%slip": "F",
        "patch_ib(1)%mass": 100.0,
        "patch_ib(1)%vel(1)": -1.0,
        "patch_ib(1)%moving_ibm": 2,
        # Patch 2
        "patch_ib(2)%geometry": 8,
        "patch_ib(2)%x_centroid": 0.0,
        "patch_ib(2)%y_centroid": 0.0,
        "patch_ib(2)%z_centroid": 0.0,
        "patch_ib(2)%radius": radius,
        "patch_ib(2)%slip": "F",
        "patch_ib(2)%mass": 100.0,
        "patch_ib(2)%moving_ibm": 2,
        # Collisions!
        "collision_model": 1,
        "ib_coefficient_of_friction": 0.0988,
        "collision_time": collision_time,
        "coefficient_of_restitution": 0.98,
        # Fluids Physical Parameters
        "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00), # 2.50(Not 1.40)
        "fluid_pp(1)%pi_inf": 0,
        }
    )
)