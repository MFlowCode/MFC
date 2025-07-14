import math
import json

# Numerical setup
Nx = 99
Ny = 99
dx = 8.0 / (1.0 * (Nx + 1))

alf_st = 0.4

p_inf = 101325
rho_inf = 1
gam = 1.4

c = math.sqrt(gam * (p_inf) / rho_inf)
cfl = 0.3
mydt = cfl * dx / c
Tfinal = 80 * 1 / c
Nt = int(Tfinal / mydt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -4,
            "x_domain%end": 4,
            "y_domain%beg": -4,
            "y_domain%end": 4,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": mydt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": int(Nt / 100),
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
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
            "bc_x%beg": -8,
            "bc_x%end": -8,
            "bc_y%beg": -8,
            "bc_y%end": -8,
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
            "patch_icpp(1)%length_x": 8.0,
            "patch_icpp(1)%length_y": 8.0,
            "patch_icpp(1)%vel(1)": 0,
            "patch_icpp(1)%vel(2)": 0,
            "patch_icpp(1)%pres": p_inf,
            "patch_icpp(1)%alpha_rho(1)": rho_inf,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%x_centroid": 0,
            "patch_icpp(2)%y_centroid": 0,
            "patch_icpp(2)%radius": 1.0,
            "patch_icpp(2)%vel(1)": 0,
            "patch_icpp(2)%vel(2)": 0,
            "patch_icpp(2)%pres": 0.0,
            "patch_icpp(2)%alpha_rho(1)": 0.0,
            "patch_icpp(2)%hcid": 281,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%alter_patch(1)": "T",
            # CBC Inflow / Outflow
            "bc_x%grcbc_in": "F",
            "bc_x%grcbc_out": "T",
            "bc_x%grcbc_vel_out": "F",
            "bc_x%pres_out": p_inf,
            "bc_y%grcbc_in": "F",
            "bc_y%grcbc_out": "T",
            "bc_y%grcbc_vel_out": "F",
            "bc_y%pres_out": p_inf,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
