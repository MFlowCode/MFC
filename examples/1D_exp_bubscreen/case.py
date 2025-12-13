#!/usr/bin/env python3
import math
import json

# FLUID PROPERTIES
R_uni = 8314.0  # [J/kmol/K]
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
cp_v = 2.1e3  # [J/kg/K]
R_v = R_uni / M_v  # [J/kg/K]

# Air
gam_g = 1.4  # [1]
M_g = 28.97  # [g/mol]
mu_g = 1.8e-05  # [kg/m/s]
k_g = 0.02556  # [W/m/K]
cp_g = 1.0e3  # [J/kg/K]
R_g = R_uni / M_g  # [J/kg/K]

# Bubble
R0ref = 17.0e-05  # [m]
p0ref = 101325.0  # [N/m2]
rho0ref = rho_l  # [kg/m3]
T0ref = 293.15  # [K]
vf0 = 0.0023  # [1]
nb = 1

# External condition
p0ext = p0ref  # [N/m2]

# Pressure wave
pa = 0.001e6  # [N/m2]
freq = 300000.0  # [Hz]
cphysical = 1475.0  # [m/s]
wavelength = cphysical / freq  # [m]

# Reference values
x0 = R0ref  # [m]
p0 = p0ref  # [N/m2]
rho0 = rho0ref  # [kg/m3]
T0 = T0ref  # [K]
u0 = math.sqrt(p0 / rho0)  # [m/s]
t0 = x0 / u0  # [s]

# CFL number should be < 1  for numerical stability
cfl = 0.2
Nx = 500
Ldomain = 0.3
L = Ldomain / x0
dx = L / float(Nx)
dt = cfl * dx / (cphysical / u0)

Lpulse = 0.3 * Ldomain
Tpulse = Lpulse / cphysical
Tfinal = 0.3 * 0.25 * 10.0 * Tpulse * u0 / x0
Nt = int(Tfinal / dt)

Nfiles = 100.0
Nout = int(math.ceil(Nt / Nfiles))
Nt = int(Nout * Nfiles)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -0.15 / x0,
            "x_domain%end": 0.15 / x0,
            "stretch_x": "F",
            "cyl_coord": "F",
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": Nout,
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
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -8,
            "bc_x%end": -8,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "fd_order": 1,
            #'schlieren_wrt'                :'T',
            "probe_wrt": "T",
            "num_probes": 1,
            "probe(1)%x": 0.0,
            # Patch 1 _ Background
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%length_x": 0.3 / x0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": p0ext / p0,
            "patch_icpp(1)%alpha_rho(1)": (1.0 - 1.0e-12) * (1.0e03 / rho0),
            "patch_icpp(1)%alpha(1)": 1.0e-12,
            "patch_icpp(1)%r0": 1.0,
            "patch_icpp(1)%v0": 0.0e00,
            # Patch 2 Screen
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%length_x": 0.08 / x0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": p0ext / p0,
            "patch_icpp(2)%alpha_rho(1)": (1.0 - vf0) * 1.0e03 / rho0,
            "patch_icpp(2)%alpha(1)": vf0,
            "patch_icpp(2)%r0": 1.0,
            "patch_icpp(2)%v0": 0.0e00,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_l - 1.0e00),
            "fluid_pp(1)%pi_inf": gam_l * (pi_inf_l / p0) / (gam_l - 1.0),
            # Bubbles
            "bubbles_euler": "T",
            "bubble_model": 2,
            "polytropic": "T",
            "polydisperse": "F",
            "nb": 1,
            # Bubble parameters
            "bub_pp%R0ref": R0ref / x0,
            "bub_pp%p0ref": p0ref / p0,
            "bub_pp%rho0ref": rho0ref / rho0,
            "bub_pp%ss": ss / (rho0 * x0 * u0 * u0),
            "bub_pp%pv": pv / p0,
            "bub_pp%mu_l": mu_l / (rho0 * x0 * u0),
            "bub_pp%gam_g": gam_g,
            # Acoustic source
            "acoustic_source": "T",
            "num_source": 1,
            "acoustic(1)%support": 1,
            "acoustic(1)%loc(1)": -0.05 / x0,
            "acoustic(1)%npulse": 1,
            "acoustic(1)%dir": 1.0,
            "acoustic(1)%pulse": 1,
            "acoustic(1)%mag": pa / p0,
            "acoustic(1)%wavelength": wavelength / x0,
        }
    )
)
