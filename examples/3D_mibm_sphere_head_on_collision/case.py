import json
import math

# This case file is based on verification cases run in
# "Analytical solution for the problem of frictional-elastic collisions
# of spherical particles using the linear model"
# by Francesco Paolo Di Maio and Alberto Di Renzo

# This can be used to check collision rebound angles to recreate
# figure 4 of that  paper

Mu = 1.84e-05
gam_a = 1.4

# lead-up-properties
velocity = 3.9
dt = 5.0e-6
collision_time = 20.0 * dt

# parerticle properties
radius = 5.0
collision_angle_degrees = 30.0
collision_angle_radians = collision_angle_degrees * math.pi / 180.0
domain_size = 4.0 * radius
lead_distance = 0.2 * radius

# simulation runs long enough to collide and travel about lead distance away again
simulation_time = 2.0 * (lead_distance / velocity) + collision_time
num_time_steps = int(simulation_time / dt)
num_saves = 10
t_step_save = int(num_time_steps / num_saves)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -0.5 * domain_size,
            "x_domain%end": 0.5 * domain_size,
            "y_domain%beg": 0.0,
            "y_domain%end": domain_size,
            "z_domain%beg": -0.5 * domain_size,
            "z_domain%end": 0.5 * domain_size,
            "cyl_coord": "F",
            "m": 60,
            "n": 60,
            "p": 60,
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
            "bc_y%beg": -15,
            "bc_y%end": -15,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            # Set IB to True and add 1 patch
            "ib": "T",
            "num_ibs": 1,
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
            "patch_icpp(1)%y_centroid": 0.5 * domain_size,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_x": domain_size,
            "patch_icpp(1)%length_y": domain_size,
            "patch_icpp(1)%length_z": domain_size,
            # Specify the patch primitive variables
            "patch_icpp(1)%vel(1)": 0.0e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%vel(3)": 0.0e00,
            "patch_icpp(1)%pres": 1.0e00,
            "patch_icpp(1)%alpha_rho(1)": 1.0e00,
            "patch_icpp(1)%alpha(1)": 1.0e00,
            # Patch: Sphere Immersed Boundary
            "patch_ib(1)%geometry": 8,
            "patch_ib(1)%x_centroid": -1.0 * lead_distance * math.sin(collision_angle_radians),  # get a lead up distance to the collision
            "patch_ib(1)%y_centroid": radius + lead_distance * math.sin(collision_angle_radians),
            "patch_ib(1)%z_centroid": 0.0,
            "patch_ib(1)%radius": radius,
            "patch_ib(1)%slip": "F",
            "patch_ib(1)%mass": 1.0e6,  # arbitrarily high mass to ignore fluid
            "patch_ib(1)%vel(1)": velocity * math.sin(collision_angle_radians),
            "patch_ib(1)%vel(2)": -velocity * math.cos(collision_angle_radians),
            "patch_ib(1)%moving_ibm": 2,
            # Collisions
            "collision_model": 1,  # soft-sphere collision model
            "ib_coefficient_of_friction": 0.092,
            "collision_time": collision_time,
            "coefficient_of_restitution": 0.98,  # almost perfectly elastic
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50(Not 1.40)
            "fluid_pp(1)%pi_inf": 0,
        }
    )
)
