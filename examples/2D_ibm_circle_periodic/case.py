import json
import math
import numpy as np

gam_a = 1.4

# sphere diameter
D = 0.1

# circle density
rho_c = 200

# domain length
L = 10 * D

# mach number
M = 1.2

# reynolds number
Re = 400.0

# pressure
P = 101325

# density
rho = 1.0

# fluid velocity
v1 = M * np.sqrt(gam_a * P / rho)

# dynamic viscosity
mu = rho * v1 * D / Re

dt = 1.0e-06
Nt = 5000 
t_save = 20 

Nx = 511  # to fully resolve requires ~ 40-60 cells across sphere diameter
Ny = Nx
Nz = 0

# immersed boundary dictionary
ib_dict = {}
ib_dict.update(
    {
        f"patch_ib({1})%geometry": 2,
        f"patch_ib({1})%x_centroid": 0.5 * L,
        f"patch_ib({1})%y_centroid": 0.5 * L,
        f"patch_ib({1})%radius": D / 2,
        f"patch_ib({1})%slip": "F",
        f"patch_ib({1})%moving_ibm": 0,
        f"patch_ib({1})%mass": rho_c * np.pi * (D / 2.0) ** 2,
    }
)

# Configuring case dictionary
case_dict = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    # x direction
    "x_domain%beg": -5.0 * D,
    "x_domain%end": 5.0 * D,
    # y direction
    "y_domain%beg": -5.0 * D,
    "y_domain%end": 5.0 * D,
    "cyl_coord": "F",
    "m": Nx,
    "n": Ny,
    "p": Nz,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": t_save,
    # Simulation Algorithm Parameters
    # Only one patch is necessary for one fluid
    "num_patches": 1,
    # Use the 5 equation model
    "model_eqns": 2,
    # 6 equations model does not need the K \div(u) term
    "alt_soundspeed": "F",
    # One fluids: air
    "num_fluids": 1,
    # time step
    "mpp_lim": "F",
    # Correct errors when computing speed of sound
    "mixture_err": "T",
    # Use TVD RK3 for time marching
    "time_stepper": 3,
    # Reconstruct the primitive variables to minimize spurious
    # Use WENO5
    "weno_order": 5,
    "weno_eps": 1.0e-14,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": 2,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    # periodic bc
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    # Set IB to True and add 1 patch for every sphere
    "ib": "T",
    "num_ibs": 1,
    "viscous": "T",
    # Formatted Database Files Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "E_wrt": "T",
    "parallel_io": "T",
    # Patch: square filled with air
    "patch_icpp(1)%geometry": 3,
    # Uniform properties, centroid is at the center of the domain
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": 10 * D,
    "patch_icpp(1)%length_y": 10 * D,
    # Specify the patch primitive variables
    "patch_icpp(1)%vel(1)": v1,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": P,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0e00,
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50(Not 1.40)
    "fluid_pp(1)%pi_inf": 0,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    # ibs wrap around domain
    "periodic_ibs": "T",
}

case_dict.update(ib_dict)

print(json.dumps(case_dict))
