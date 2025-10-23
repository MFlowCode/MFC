import json
import math

Mu = 1.84e-05
gam_a = 1.4
gam_b = 1.1

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            # axial direction
            "x_domain%beg": 0.0e00,
            "x_domain%end": 6.0e-03,
            # r direction
            "y_domain%beg": 0.0e00,
            "y_domain%end": 6.0e-03,
            "cyl_coord": "F",
            "m": 200,
            "n": 200,
            "p": 0,
            "dt": 0.57e-5 / 2,
            "t_step_start": 0,
            "t_step_stop": 7200,  # 3000
            "t_step_save": 30,  # 10
            # Simulation Algorithm Parameters
            # Only one patches are necessary, the air tube
            "num_patches": 1,
            # Use the 5 equation model
            "model_eqns": 2,
            # 6 equations model does not need the K \div(u) term
            "alt_soundspeed": "F",
            # One fluids: air
            "num_fluids": 2,
            # time step
            "mpp_lim": "F",
            # Correct errors when computing speed of sound
            "mixture_err": "T",
            # Use TVD RK3 for time marching
            "time_stepper": 3,
            # Reconstruct the primitive variables to minimize spurious
            # Use WENO5
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "T",
            "avg_state": 2,
            # Use the mapped WENO weights to maintain monotinicity
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            # Use the HLLC  Riemann solver
            "riemann_solver": 2,
            "wave_speeds": 1,
            # We use reflective boundary conditions at octant edges and
            # non-reflective boundary conditions at the domain edges
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
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
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 3.0e-03,
            # Uniform medium density, centroid is at the center of the domain
            "patch_icpp(1)%y_centroid": 3.0e-03,
            "patch_icpp(1)%length_x": 6.0e-03,
            "patch_icpp(1)%length_y": 6.0e-03,
            # Specify the patch primitive variables
            "patch_icpp(1)%vel(1)": 0.1e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": 1.0e00,
            "patch_icpp(1)%alpha_rho(1)": 0.8e00,
            "patch_icpp(1)%alpha(1)": 0.8e00,
            "patch_icpp(1)%alpha_rho(2)": 0.2e00,
            "patch_icpp(1)%alpha(2)": 0.2e00,
            # Patch: Airfoil Immersed Boundary
            "patch_ib(1)%geometry": 4,
            "patch_ib(1)%x_centroid": 1.0e-03,
            "patch_ib(1)%y_centroid": 3.0e-03,
            "patch_ib(1)%c": 1.0e-03,
            "patch_ib(1)%t": 0.15,
            "patch_ib(1)%p": 0.4,
            "patch_ib(1)%m": 0.02,
            "patch_ib(1)%angles(3)": -0.5235987756,  # 30 degrees clockwise rotation, in radians
            # Fluids Physical Parameters
            # Use the same stiffness as the air bubble
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50 (Not 1.40)
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(2)%gamma": 1.0e00 / (gam_b - 1.0e00),  # 2.50 (Not 1.40)
            "fluid_pp(2)%pi_inf": 0,
        }
    )
)
