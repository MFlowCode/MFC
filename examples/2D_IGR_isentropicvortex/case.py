import json
import math

# Parameters
epsilon = "5d0"
alpha = "1d0"
gamma = "1.4"


# Numerical setup
Nx = 149
dx = 20.0 / (1.0 * (Nx + 1))

c = 1.4**0.5
C = 0.3
mydt = C * dx / c
Nt = int(20 / (0.1 * mydt))

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -10,
            "x_domain%end": 10,
            "y_domain%beg": -10,
            "y_domain%end": 10,
            "stretch_x": "F",
            "stretch_y": "F",
            "loops_x": 2,
            "loops_y": 2,
            "a_x": 1.03,
            "a_y": 1.03,
            "x_a": -1.5,
            "y_a": -1.5,
            "x_b": 1.5,
            "y_b": 1.5,
            "m": Nx,
            "n": Nx,
            "p": 0,
            "dt": mydt,
            "t_step_start": 0,
            "t_step_stop": int(30 * Nt),
            "t_step_save": int(Nt),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "F",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 5,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "igr": "T",
            "igr_order": 5,
            "igr_iter_solver": 1,
            "num_igr_iters": 5,
            "num_igr_warm_start_iters": 150,
            "alf_factor": 2,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "omega_wrt(3)": "T",
            "fd_order": 2,
            # Patch 1
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0,
            "patch_icpp(1)%y_centroid": 0,
            "patch_icpp(1)%length_x": 20.0,
            "patch_icpp(1)%length_y": 20.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 0.0,
            "patch_icpp(1)%alpha_rho(1)": 0.0,
            "patch_icpp(1)%hcid": 283,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
