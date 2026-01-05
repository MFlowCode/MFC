import json
import math

Mu = 0.01 * 1e-4 # kenimatic viscosity of 0.01 cm^2/s
gam_a = 1.4

total_time = 0.35
# total_time = 2.0
num_time_steps = 100000
dt = float(total_time / num_time_steps)
num_saves = 1000
steps_per_save = int(num_time_steps / num_saves)

length = 0.4e-2
avg_speed = 5.5*length / total_time

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            # For these computations, the cylinder is placed at the (0,0,0)
            # domain origin.
            "x_domain%beg": 0.0e-2,
            "x_domain%end": length,
            "y_domain%beg": 0.0e-2,
            "y_domain%end": 7.*length,
            "cyl_coord": "F",
            "m": 100,
            "n": 700,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": num_time_steps,  # 3000
            "t_step_save": steps_per_save,  # 10
            # Simulation Algorithm Parameters
            # Only one patches are necessary, the air tube
            "num_patches": 1,
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
            "bc_x%beg": -15,
            "bc_x%end": -15,
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
            # Specify the rectangle patch geometry
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": length/2.,
            # Uniform medium density, centroid is at the center of the domain
            "patch_icpp(1)%y_centroid": 7.*length/2.,
            "patch_icpp(1)%length_x": length,
            "patch_icpp(1)%length_y": 7. * length,
            # Specify the patch primitive variables
            "patch_icpp(1)%vel(1)": 0.00e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": f"1.0 + 9.81*{7. * length} - 9.81*y", # Set up a linear pressure gradient to start with # 1.0e00 at the top,
            "patch_icpp(1)%alpha_rho(1)": 1.0e00,
            "patch_icpp(1)%alpha(1)": 1.0e00,
            # Patch: Cylinder Immersed Boundary
            "patch_ib(1)%geometry": 6,
            "patch_ib(1)%moving_ibm": 2,
            "patch_ib(1)%x_centroid": length/2.,
            "patch_ib(1)%y_centroid": 6. * length,
            "patch_ib(1)%length_x": length / 4.,
            "patch_ib(1)%length_y": length / 8,
            "patch_ib(1)%vel(2)": 0.0e00,
            "patch_ib(1)%angles(3)": math.pi / 4.,
            "patch_ib(1)%mass": 1.1 * math.pi * length * length / 128., # density of 1.1 times the volume of the ellipse to give density ratio of 1.1
            "patch_ib(1)%slip": "F",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50(Not 1.40)
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(1)%Re(1)": 1. / Mu,
            # Body Forces
            "bf_y": "T",
            "k_y": 0.0,
            "w_y": 0.0,
            "p_y": 0.0,
            "g_y": -9.81,  # gravity
        }
    )
)
