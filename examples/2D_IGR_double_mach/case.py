import json
import math

# Numerical setup
Nx = 799
Ny = int((Nx + 1) / 4) - 1
dx = 4.0 / (1.0 * (Nx + 1))

c = 1.0
cfl = 0.1
mydt = cfl * dx / (11.0 * c)
Nt = int(0.2 / (1 * mydt))

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0,
            "x_domain%end": 4,
            "y_domain%beg": 0,
            "y_domain%end": 1,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": mydt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": int(Nt / 100),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 1,
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -17,
            "bc_x%end": -3,
            "bc_y%beg": -2,
            "bc_y%end": -3,
            "igr": "T",
            "igr_order": 3,
            "igr_iter_solver": 1,
            "num_igr_iters": 5,
            "num_igr_warm_start_iters": 150,
            "alf_factor": 10,
            "double_mach": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "omega_wrt(3)": "F",
            "fd_order": 2,
            # Patch 1
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 2.0,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 4.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.4,
            "patch_icpp(1)%hcid": 285,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
