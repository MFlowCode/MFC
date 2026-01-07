import json
import math

# FLUID PROPERTIES
R_uni = 8314.0  # [J/mol/K]
# Water
gam_l = 7.1  # [1]
pi_inf_l = 3.43e5  # [N/m2]
rho_l = 1.0e03  # [kg/m3]
mu_l = 1.002e-03  # [kg/m/s]
ss = 0.07275  # [kg/s2]
pv = 2.3388e03  # [N/m2]
vd = 0.242e-4  # []

# Vapor
gam_v = 1.33  # [1]
M_v = 18.02  # [g/mol]
mu_v = 0.8816e-05  # [kg/m/s]
k_v = 0.019426  # []
cp_v = 2.1e3  # [J/kg/K]
R_v = R_uni / M_v  # []

# Air
gam_g = 1.4  # [1]
M_g = 28.97  # [g/mol]
mu_g = 1.8e-05  # [kg/m/s]
k_g = 0.02556  # []
cp_g = 1.0e3  # [J/kg/K]
R_g = R_uni / M_g  # []

# Bubble
R0ref = 10.0e-06  # [m]
p0ref = 101325.0  # [N/m2]
rho0ref = rho_l  # [kg/m3]
T0ref = 293.15  # [K]
vf0 = 1e-3  # [1]
nb = 1

# External condition
p0ext = 101325.0  # [N/m2]

# Reference values
x0 = R0ref  # [m]
p0 = p0ref  # [N/m2]
rho0 = rho0ref  # [kg/m3]
T0 = T0ref  # [K]
u0 = math.sqrt(p0 / rho0)  # [m/s]
t0 = x0 / u0  # [s]

###
cact = math.sqrt(gam_l * (p0 + pi_inf_l) / ((1 - vf0) * rho0))
cfl = 0.3
Nx = 400
Ny = 200
dx = 6.0e-03 / (x0 * float(Nx))
dt = cfl * dx * u0 / cact

Min = 1.2
beta = (gam_l + 1) * Min**2 / (Min**2 * (gam_l - 1 + 2 * vf0) + 2 * (1 - vf0))
vel = Min * cact * (beta - 1) / beta
delta = (1 - vf0) + gam_l * Min**2 * (beta - 1) * (1 + pi_inf_l / p0) / beta

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "F",
            # Computational Domain Parameters
            "x_domain%beg": 0.0e00,
            "x_domain%end": 6.0e-03 / x0,
            "y_domain%beg": 0.0e00,
            "y_domain%end": 3.0e-03 / x0,
            "cyl_coord": "F",
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": 1000,
            "t_step_save": 10,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "T",
            "avg_state": 2,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "bc_x%beg": -6,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # Set IB to True and add 1 patch
            "ib": "T",
            "num_ibs": 1,
            # Formatted Database Files Structure Parameters
            # Export primitive variables in double precision with parallel
            # I/O to minimize I/O computational time during large simulations
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "fd_order": 1,
            "omega_wrt(3)": "T",
            "parallel_io": "T",
            # Ambient State
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 3.0e-03 / x0,
            "patch_icpp(1)%y_centroid": 1.50e-03 / x0,
            "patch_icpp(1)%length_x": 6.0e-03 / x0,
            "patch_icpp(1)%length_y": 3.0e-03 / x0,
            "patch_icpp(1)%alpha_rho(1)": (1.0 - vf0),
            "patch_icpp(1)%alpha(1)": vf0,
            "patch_icpp(1)%vel(1)": 1.5,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": 1.0e00,
            "patch_icpp(1)%r0": 1.0,
            "patch_icpp(1)%v0": 0.0e00,
            # Shocked State
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.5e-03 / x0,
            "patch_icpp(2)%y_centroid": 1.50e-03 / x0,
            "patch_icpp(2)%length_x": 1.0e-03 / x0,
            "patch_icpp(2)%length_y": 3.0e-03 / x0,
            "patch_icpp(2)%alpha_rho(1)": beta,
            "patch_icpp(2)%alpha(1)": beta * vf0,
            "patch_icpp(2)%vel(1)": vel / u0,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%pres": delta,
            "patch_icpp(2)%r0": 1.0,
            "patch_icpp(2)%v0": 0.0e00,
            "patch_icpp(2)%alter_patch(1)": "T",
            # Patch: Cylinder Immersed Boundary
            "patch_ib(1)%geometry": 4,
            "patch_ib(1)%x_centroid": 1.5e-03 / x0,
            "patch_ib(1)%y_centroid": 1.5e-03 / x0,
            "patch_ib(1)%c": 1.0e-03 / x0,
            "patch_ib(1)%t": 0.15,
            "patch_ib(1)%p": 0.4,
            "patch_ib(1)%m": 0.02,
            "patch_ib(1)%slip": "F",
            "patch_ib(1)%theta": 15,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_l - 1.0e00),
            "fluid_pp(1)%pi_inf": gam_l * (pi_inf_l) / (gam_l - 1.0),
            # Bubbles
            "bubbles_euler": "T",
            "bubble_model": 2,
            "polytropic": "T",
            "polydisperse": "F",
            "poly_sigma": 0.3,
            "thermal": 3,
            # Bubble parameters
            "bub_pp%R0ref": R0ref / x0,
            "bub_pp%p0ref": p0ref / p0,
            "bub_pp%rho0ref": rho0ref / rho0,
            "bub_pp%ss": ss / (rho0 * x0 * u0 * u0),
            "bub_pp%pv": pv / p0,
            "bub_pp%mu_l": mu_l / (rho0 * x0 * u0),
            "bub_pp%gam_g": gam_g,
        }
    )
)
