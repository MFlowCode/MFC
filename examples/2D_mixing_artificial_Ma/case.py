#!/usr/bin/env python3
import math
import json

# FLUID PROPERTIES (WATER, VAPOR & AIR)
# Water
rho_w = 1.0e03  # [kg/m3] density of water
gamma_w = 7.1  # [1] specific heat ratio
pi_inf_w = 3.06e08  # [N/m2] water stiffness

# REFERENCE VALUES
x_ref = 0.002475  # [m] initial vorticity thickness
u_ref = 3.4343  # [m/s] upper stream velocity

# NON-DIMENSIONAL NUMBERS
Re = 50.0  # [1] Reynolds number based on the upper stream
# velocity and the initial vorticity thickness

# SPEED OF SOUND
pres = 8236.0  # [N/m2] pressure of water
c_w = math.sqrt(gamma_w * (pres + pi_inf_w) / rho_w)

# MODIFIED PROPERTIES FOR ARTIFICIAL MACH NUMBER
Ma_t = 0.1  # [1] target Mach number
pi_fac = (rho_w * u_ref**2 / (gamma_w * Ma_t**2) - pres) / pi_inf_w
pi_inf_w = pi_inf_w * pi_fac
c_w = math.sqrt(gamma_w * (pres + pi_inf_w) / rho_w)

# SIMULATION SETUP
# Domain size
Lx = 59.0
Ly = 59.0

# Number of grid cells
Nx = 191
Ny = 191

# Grid spacing
dx = Lx / float(Nx + 1)
dy = Ly / float(Ny + 1)

# Time advancement
cfl = 5e-1
T = 20.0
dt = cfl * dx / (1.0 + c_w / u_ref)
Ntfinal = int(T / dt)
Ntstart = 0
Nfiles = 20
t_save = int(math.ceil((Ntfinal - 0) / float(Nfiles)))
Nt = t_save * Nfiles
t_step_start = Ntstart
t_step_stop = int(Nt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": Lx,
            "y_domain%beg": -Ly / 2.0,
            "y_domain%end": Ly / 2.0,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": dt,
            "t_step_start": t_step_start,
            "t_step_stop": t_step_stop,
            "t_step_save": t_save,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "T",
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -6,
            "bc_y%end": -6,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "cons_vars_wrt": "T",
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "fd_order": 1,
            "omega_wrt(3)": "T",
            # Patch 1
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": Lx / 2.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%length_y": Ly,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%vel(1)": 1.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": pres / (rho_w * u_ref**2),
            # Mixing layer
            "mixlayer_vel_profile": "T",
            "mixlayer_vel_coef": 1.0,
            # Artificial Mach number
            "pi_fac": pi_fac,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gamma_w - 1.0e00),
            "fluid_pp(1)%pi_inf": gamma_w * (pi_inf_w / (rho_w * u_ref**2)) / (gamma_w - 1.0e00),
            "fluid_pp(1)%Re(1)": Re,
        }
    )
)
