import json
import math

gam = 1.4

D = 0.1
N = 400

tf = 0.01
saveFreq = tf / 200

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            # x direction
            "x_domain%beg": -5.0 * D,
            "x_domain%end": 5.0 * D,
            # y direction
            "y_domain%beg": -2.5 * D,
            "y_domain%end": 2.5 * D,
            # z direction
            "z_domain%beg": -2.5 * D,
            "z_domain%end": 2.5 * D,
            "cyl_coord": "F",
            "m": 2 * N,
            "n": N,
            "p": N,
            "cfl_adap_dt": "T",
            "cfl_target": 0.4,
            "n_start": 0,
            "t_save": saveFreq,
            "t_stop": tf,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            # Use the 5 equation model
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            # Use IGR5
            "igr": "T",
            "igr_pres_lim": "T",
            "igr_order": 5,
            "alf_factor": 10,
            "igr_iter_solver": 2,
            "num_igr_iters": 15,
            "num_igr_warm_start_iters": 50,
            # We use ghost-cell extrapolation
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            # Set IB to True
            "ib": "T",
            "num_ibs": 1,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "E_wrt": "T",
            "parallel_io": "T",
            # Main Patch: 10D:5D:5D rectangular prism centered at the (0, 0, 0)
            # HCID 390 smooths the x-velocity profile in a small region around
            # the sphere using a tanh profile.
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%hcid": 390,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_x": 10 * D,
            "patch_icpp(1)%length_y": 5 * D,
            "patch_icpp(1)%length_z": 5 * D,
            "patch_icpp(1)%vel(1)": 527.2e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%vel(3)": 0.0e00,
            "patch_icpp(1)%pres": 10918.2549,
            "patch_icpp(1)%alpha_rho(1)": 0.2199,
            "patch_icpp(1)%alpha(1)": 1.0e00,
            # Patch: Sphere Immersed Boundary
            "patch_ib(1)%geometry": 8,
            "patch_ib(1)%x_centroid": -3.0 * D,
            "patch_ib(1)%y_centroid": 0.0,
            "patch_ib(1)%z_centroid": 0.0,
            "patch_ib(1)%radius": D / 2,
            "patch_ib(1)%slip": "T",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam - 1.0e00),
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(1)%Re(1)": 650000,
        }
    )
)
