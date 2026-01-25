import json
import math

# This case is a recreation of the case from "Moving overlapping grids with adaptive mesh refinement for high-speed reactive and non-reactive flow"
# by William D. Henshaw and Donald W. Schwendeman

# fluid parameters
gam_a = 1.4

# domain size and speed
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

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            # For these computations, the cylinder is placed at the (0,0,0)
            # domain origin.
            # axial direction
            "x_domain%beg": -domain_size * 0.5,
            "x_domain%end": domain_size * 0.5,
            # r direction
            "y_domain%beg": -domain_size * 0.5,
            "y_domain%end": domain_size * 0.5,
            "cyl_coord": "F",
            "m": 1000,
            "n": 1000,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": num_time_steps,  # 10000,
            "t_step_save": steps_to_save,
            # Simulation Algorithm Parameters
            # Only one patches are necessary, the air tube
            "num_patches": 2,
            # Use the 5 equation model
            "model_eqns": 2,
            "alt_soundspeed": "F",
            # One fluids: air
            "num_fluids": 1,
            # time step
            "mpp_lim": "F",
            # Correct errors when computing speed of sound
            "mixture_err": "T",
            # Use TVD RK3 for time marching
            "time_stepper": 3,
            # Use WENO5
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
            # We use ghost-cell
            "bc_x%beg": -17,
            "bc_x%end": -8,
            "bc_y%beg": -15,
            "bc_y%end": -15,
            # Set IB to True and add 1 patch
            "ib": "T",
            "num_ibs": 1,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "E_wrt": "T",
            "parallel_io": "T",
            # Patch: Constant Tube filled with air
            # Specify the cylindrical air tube grid geometry
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(2)%geometry": 3,
            # patch locations
            "patch_icpp(1)%x_centroid": 0.5 * wave_front + 0.25 * domain_size,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 0.5 * domain_size - wave_front,
            "patch_icpp(1)%length_y": domain_size,
            "patch_icpp(2)%x_centroid": 0.5 * wave_front - 0.25 * domain_size,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%length_x": 0.5 * domain_size + wave_front,
            "patch_icpp(2)%length_y": domain_size,
            # Specify the patch primitive variables
            "patch_icpp(1)%vel(1)": pre_shock_speed,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": pre_shock_pressure,
            "patch_icpp(1)%alpha_rho(1)": pre_shock_density,
            "patch_icpp(1)%alpha(1)": 1.0e00,
            "patch_icpp(2)%vel(1)": post_shock_speed,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": post_shock_pressure,
            "patch_icpp(2)%alpha_rho(1)": post_shock_density,
            "patch_icpp(2)%alpha(1)": 1.0e00,
            # Patch: Cylinder Immersed Boundary
            "patch_ib(1)%geometry": 2,
            "patch_ib(1)%x_centroid": -0.5,
            "patch_ib(1)%y_centroid": 0.0,
            "patch_ib(1)%radius": 0.5,
            "patch_ib(1)%slip": "T",
            "patch_ib(1)%moving_ibm": 2,
            "patch_ib(1)%vel(1)": 0,
            "patch_ib(1)%vel(2)": 0,
            "patch_ib(1)%vel(3)": 0,
            "patch_ib(1)%angles(1)": 0.0,  # x-axis rotation in radians
            "patch_ib(1)%angles(2)": 0.0,  # y-axis rotation
            "patch_ib(1)%angles(3)": 0.0,  # z-axis rotation
            "patch_ib(1)%angular_vel(1)": 0.0,  # x-axis rotational velocity in radians per second
            "patch_ib(1)%angular_vel(2)": 0.0,  # y-axis rotation
            "patch_ib(1)%angular_vel(3)": 0.0,  # z-axis rotation
            "patch_ib(1)%mass": 0.5,  # z-axis rotation
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50(Not 1.40)
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(1)%Re(1)": 2500000,
        }
    )
)
