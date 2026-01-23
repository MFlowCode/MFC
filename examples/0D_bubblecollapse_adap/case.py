#!/usr/bin/env python
import math
import json

# FLUID PROPERTIES
# Water
gam_l = 7.1  # [1]
pi_inf_l = 306.0e06  # [N/m2]
rho_l = 1.0e03  # [kg/m3]
mu_l = 1.002e-03  # [kg/m/s]
ss = 0.07275  # [kg/s2]
pv = 2.3388e03  # [N/m2]

# Vapor
gam_v = 1.33  # [1]
M_v = 18.02  # [g/mol]
mu_v = 0.8816e-05  # [kg/m/s]
k_v = 0.019426  # [W/m/K]

# Air
gam_g = 1.4  # [1]
M_g = 28.97  # [g/mol]
mu_g = 1.8e-05  # [kg/m/s]
k_g = 0.02556  # [W/m/K]

# BUBBLES
R0ref = 50.0e-06  # [m]
p0ref = 8236.0  # [N/m2]
rho0ref = rho_l  # [kg/m3]
vf0 = 1e-5  # [1]
nb = 1

# External pressure
p0ext = 1000 * p0ref  # [N/m2]

# REFERENCE VALUES
x0 = R0ref  # [m]
p0 = p0ref  # [N/m2]
rho0 = rho0ref  # [kg/m3]
u0 = math.sqrt(p0 / rho0)  # [m/s]
t0 = x0 / u0  # [s]
tc = 0.915 * math.sqrt(rho_l * R0ref**2 / (p0ext - pv))  # [s]

# DOMAIN
Nx = 30
L = 20.0e-03  # [m]

# TIME STEPS
Tfinal = 1.5 * tc / t0
Nt = int(5e2 + 1)
t_save = 1
dt = 0.0001

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -0.5 * L / x0,
            "x_domain%end": 0.5 * L / x0,
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": t_save,
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
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "fd_order": 1,
            # Patch 1 _ Background
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%length_x": L / x0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": p0ext / p0,
            "patch_icpp(1)%alpha_rho(1)": (1.0 - vf0),
            "patch_icpp(1)%alpha(1)": vf0,
            "patch_icpp(1)%r0": R0ref / x0,
            "patch_icpp(1)%v0": 0.0,
            # Bubbles
            "bubbles_euler": "T",
            "bubble_model": 2,
            # Nondimensional numbers
            "bub_pp%R0ref": R0ref / x0,
            "bub_pp%p0ref": p0ref / p0,
            "bub_pp%rho0ref": rho0ref / rho0,
            "bub_pp%ss": ss / (rho0 * x0 * u0 * u0),
            "bub_pp%pv": pv / p0,
            "bub_pp%mu_l": mu_l / (rho0 * x0 * u0),
            "bub_pp%gam_g": gam_g,
            "adv_n": "T",
            "adap_dt": "T",
            "adap_dt_tol": 1e-4,
            "adap_dt_max_iters": 10000,
            "polytropic": "T",
            "thermal": 1,
            "polydisperse": "F",
            "nb": nb,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_l - 1.0e00),
            "fluid_pp(1)%pi_inf": gam_l * (pi_inf_l / p0) / (gam_l - 1.0),
        }
    )
)
