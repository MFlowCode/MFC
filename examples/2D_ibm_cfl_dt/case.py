import json
import math

Mu = 1.84e-05
gam_a = 1.4

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
            "x_domain%beg": 0.0e00,
            "x_domain%end": 8.0e-03,
            # r direction
            "y_domain%beg": 0.0e00,
            "y_domain%end": 6.0e-03,
            "cyl_coord": "F",
            "m": 199,
            "n": 139,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.3,
            "n_start": 0,
            "t_save": 8e-5,
            "t_stop": 8e-3,
            # Simulation Algorithm Parameters
            # Only one patches are necessary, the air tube
            "num_patches": 1,
            # Use the 5 equation model
            "model_eqns": 2,
            "alt_soundspeed": "F",
            # One fluids: air
            "num_fluids": 1,
            # No need to ensure the volume fractions sum to unity at the end of each
            # time step
            "mpp_lim": "F",
            # Correct errors when computing speed of sound
            "mixture_err": "T",
            # Use TVD RK3 for time marching
            "time_stepper": 3,
            # Use WENO5
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            # 'weno_Re_flux'                 : 'T',
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
            "bc_y%beg": -2,
            "bc_y%end": -2,
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
            "c_wrt": "T",
            # Patch: Constant Tube filled with air
            # Specify the cylindrical air tube grid geometry
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 4.0e-03,
            # Uniform medium density, centroid is at the center of the domain
            "patch_icpp(1)%y_centroid": 3.0e-03,
            "patch_icpp(1)%length_x": 8.0e-03,
            "patch_icpp(1)%length_y": 6.0e-03,
            # Specify the patch primitive variables
            "patch_icpp(1)%vel(1)": 5e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": 1.0e00,
            "patch_icpp(1)%alpha_rho(1)": 1.0e00,
            "patch_icpp(1)%alpha(1)": 1.0e00,
            # Patch: Cylinder Immersed Boundary
            "patch_ib(1)%geometry": 2,
            "patch_ib(1)%x_centroid": 1.5e-03,
            "patch_ib(1)%y_centroid": 3e-03,
            "patch_ib(1)%radius": 0.4e-03,
            "patch_ib(1)%slip": "F",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50(Not 1.40)
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(1)%Re(1)": 250000,
        }
    )
)
