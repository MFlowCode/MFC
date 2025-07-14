import math
import json

# Parameters
epsilon = "5d0"
alpha = "1d0"
gamma = "1.4"

# Initial conditions
vel1_i = "0d0"
vel2_i = "0d0"
T_i = "1d0"
alpha_rho1_i = "1d0"
pres_i = "1d0"

# Perturbations
vel1 = f"{vel1_i} + (y - yc)*({epsilon}/(2d0*pi))*" + f"exp({alpha}*(1d0 - (x - xc)**2d0 - (y - yc)**2d0))"
vel2 = f"{vel2_i} - (x - xc)*({epsilon}/(2d0*pi))*" + f"exp({alpha}*(1d0 - (x - xc)**2d0 - (y - yc)**2d0))"
alpha_rho1 = (
    f"{alpha_rho1_i}*(1d0 - ({alpha_rho1_i}/{pres_i})*({epsilon}/(2d0*pi))*"
    + f"({epsilon}/(8d0*{alpha}*({gamma} + 1d0)*pi))*"
    + f"exp(2d0*{alpha}*(1d0 - (x - xc)**2d0"
    + f"- (y - yc)**2d0))"
    + f")**{gamma}"
)
pres = (
    f"{pres_i}*(1d0 - ({alpha_rho1_i}/{pres_i})*({epsilon}/(2d0*pi))*"
    + f"({epsilon}/(8d0*{alpha}*({gamma} + 1d0)*pi))*"
    + f"exp(2d0*{alpha}*(1d0 - (x - xc)**2d0"
    + f"- (y - yc)**2d0))"
    + f")**({gamma} + 1d0)"
)

# Numerical setup
Nx = 399
dx = 10.0 / (1.0 * (Nx + 1))

c = 1.4**0.5
C = 0.3
mydt = C * dx / c * 0.01
Nt = 1 / mydt

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -3,
            "x_domain%end": 3,
            "y_domain%beg": -3,
            "y_domain%end": 3,
            "stretch_x": "T",
            "stretch_y": "T",
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
            "t_step_stop": 10,
            "t_step_save": 1,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 3,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -4,
            "bc_x%end": -4,
            "bc_y%beg": -4,
            "bc_y%end": -4,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "omega_wrt(3)": "T",
            "fd_order": 2,
            # Patch 1
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0,
            "patch_icpp(1)%y_centroid": 0,
            "patch_icpp(1)%length_x": 10.0,
            "patch_icpp(1)%length_y": 10.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 0.0,
            "patch_icpp(1)%alpha_rho(1)": 0.0,
            "patch_icpp(1)%hcid": 280,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
