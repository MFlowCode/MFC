import json
import math
import numpy as np

Mu = 1.84e-05
gam_a = 1.4
R = 287.0

D = 0.1

P = 101325 # Pa
rho = 1.225 # kg/m^3

T = P/(rho*R)

M = 1.2
Re = 1500.0
v1 = M*(gam_a*P/rho)**(1.0/2.0)

mu = rho*v1*D/Re # dynamic viscosity for current case

#print('mu: ', mu)
#print('v1: ', v1)
#print('rho: ', rho)
#print('Kn = ' + str( np.sqrt(np.pi*gam_a/2)*(M/Re) )) # Kn < 0.01 = continuum flow

dt = 4.0E-06
Nt = 1000
t_save = 10

Nx = 63
Ny = 63
Nz = 63

# immersed boundary dictionary
ib_dict = {}
ib_dict.update({
    f"patch_ib({1})%geometry": 8,
    f"patch_ib({1})%x_centroid": 0.0,
    f"patch_ib({1})%y_centroid": 0.0,
    f"patch_ib({1})%z_centroid": 0.0,
    f"patch_ib({1})%radius": D / 2,
    f"patch_ib({1})%slip": "F",
    })

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
    # z direction
    "z_domain%beg": -5.0 * D,
    "z_domain%end": 5.0 * D,
    "cyl_coord": "F",
    "m": Nx,
    "n": Ny,
    "p": Nz,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,  # 3000
    "t_step_save": t_save,  # 10
    "t_step_stat_start": 100,
    # Simulation Algorithm Parameters
    # Only one patches are necessary, the air tube
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
    "low_Mach": 1,
    "wave_speeds": 1,
    # periodic bc
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "bc_z%beg": -1,
    "bc_z%end": -1,
    # Set IB to True and add 1 patch
    "ib": "T",
    "num_ibs": 1,
    "viscous": "T",
    # Formatted Database Files Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "E_wrt": "T",
    "q_filtered_wrt": "T",
    "parallel_io": "T",
    # Patch: Constant Tube filled with air
    # Specify the cylindrical air tube grid geometry
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    # Uniform medium density, centroid is at the center of the domain
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": 10 * D,
    "patch_icpp(1)%length_y": 10 * D,
    "patch_icpp(1)%length_z": 10 * D,
    # Specify the patch primitive variables
    "patch_icpp(1)%vel(1)": v1,
    "patch_icpp(1)%vel(2)": 0.0e00,
    "patch_icpp(1)%vel(3)": 0.0e00,
    "patch_icpp(1)%pres": P,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0e00,
    # Patch: Sphere Immersed Boundary
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50(Not 1.40)
    "fluid_pp(1)%pi_inf": 0,
    "fluid_pp(1)%Re(1)": 1.0 / mu,

    # new case additions
    "periodic_forcing": "T",
    "periodic_ibs": "T",
    "volume_filtering_momentum_eqn": "T",

    "u_inf_ref": v1,
    "rho_inf_ref": rho,
    "T_inf_ref": T,

    "store_levelset": "F",
    "slab_domain_decomposition": "T", 
    "compute_autocorrelation": "T",
    }

case_dict.update(ib_dict)

print(json.dumps(case_dict))
