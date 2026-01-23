import math
import json

# Numerical setup
Nx = 250
Ny = 125
dx = 40.0 / (1.0 * (Nx + 1))

M_inf = 0.3
Mv = 0.1

p_inf = 101325
rho_inf = 1
gam = 1.4

c = math.sqrt(gam * (p_inf) / rho_inf)
u_inf = M_inf * c
cfl = 0.3
mydt = cfl * dx / c
Tfinal = 150 * 1 / u_inf
Nt = int(Tfinal / mydt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -20,
            "x_domain%end": 20,
            "y_domain%beg": -10,
            "y_domain%end": 10,
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
            "bc_x%beg": -7,
            "bc_x%end": -8,
            "bc_y%beg": -6,
            "bc_y%end": -6,
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
            "patch_icpp(1)%length_x": 40.0,
            "patch_icpp(1)%length_y": 20.0,
            "patch_icpp(1)%vel(1)": u_inf,
            "patch_icpp(1)%vel(2)": 0,
            "patch_icpp(1)%pres": 101325,
            "patch_icpp(1)%alpha_rho(1)": 1,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%x_centroid": 0,
            "patch_icpp(2)%y_centroid": 0,
            "patch_icpp(2)%radius": 1.0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%hcid": 282,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%alter_patch(1)": "T",
            # CBC Inflow / Outflow
            "bc_x%grcbc_in": "T",
            "bc_x%grcbc_out": "T",
            "bc_x%grcbc_vel_out": "T",
            "bc_x%vel_in(1)": u_inf,
            "bc_x%vel_in(2)": 0,
            "bc_x%vel_in(3)": 0,
            "bc_x%pres_in": p_inf,
            "bc_x%alpha_rho_in(1)": rho_inf,
            "bc_x%alpha_in(1)": 1,
            "bc_x%vel_out(1)": u_inf,
            "bc_x%vel_out(2)": 0,
            "bc_x%vel_out(3)": 0,
            "bc_x%pres_out": p_inf,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
