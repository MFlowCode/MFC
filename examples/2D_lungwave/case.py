#!/usr/bin/env python3
import math
import json

pi = 3.141592653589
# material parameters

# material1 :: gas
# patterson 2018

gammag = 1.4  # unitless
Bg = 0  # pascals
rhog = 1.18  # kg/m^3
c_g = 347.2  # m/s
G_g = 0  # pa

# material2 :: lung

gammal = 5.5
Bl = 492.0e06
rhol = 996.0
c_l = 1648.7
G_l = 1e3

# primitive variables
patmos = 101325.0  # pa

# problem specific variable
lambda_wave = 1e-3

# define pulse
P_amp = 10.0e6
P_len = 45  # length of the impulse
theta = -math.pi / 2  # direction of propagation

# non-dim

# define characteristic density, length, time, stress material
rho_char = rhog
length_char = lambda_wave
c_char = c_g
time_char = length_char / c_char
stress_char = rho_char * c_char * c_char / gammag

# non-dim the properties
rhog_n = rhog / rho_char
c_g_n = c_g / c_char
rhol_n = rhol / rho_char
c_l_n = c_l / c_char
Bg_n = Bg / stress_char
Bl_n = Bl / stress_char
G_g_n = G_g / stress_char
G_l_n = G_l / stress_char
patmos_n = patmos / stress_char
P_amp_n = P_amp / stress_char

# geometry
dlengx = 1.0
dlengy = 20.0
Nx = 200
Ny = dlengy * Nx
dx = dlengx / Nx
dy = dlengy / Ny
alphal_back = 1.0
alphag_back = 0.0
alphal_lung = 0.0
alphag_lung = 1.0

interface_amp = 0.5

# time stepping requirements
time_end = 2.5
cfl = 0.5

dt = cfl * dx / c_l
Nt = int(time_end / dt)
Nframes = 50000
tstart = 0
tstop = Nt
tsave = int(Nt / Nframes)

# interface profile
interface_amp = 0.5

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": dlengx,
            "y_domain%beg": -dlengy / 2.0,
            "y_domain%end": dlengy / 2.0,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "dt": dt,
            "t_step_start": tstart,
            "t_step_stop": tstop,
            "t_step_save": tsave,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
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
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -6,
            "bc_y%end": -6,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Monopole setting
            "acoustic_source": "T",  # creating an acoustic wave
            "num_source": 1,  # place in the middle and expand
            "acoustic(1)%support": 2,
            "acoustic(1)%pulse": 1,  # sine  wave
            "acoustic(1)%npulse": 1,  # 1 pulse
            "acoustic(1)%wavelength": dlengx,
            "acoustic(1)%mag": 10.0 * patmos_n,  # magnitude
            "acoustic(1)%length": 1 * dlengx,  # impulse length
            "acoustic(1)%loc(1)": dlengx / 2,  # x_center of the domain
            "acoustic(1)%loc(2)": 5.0 * dlengx,  # upper boundary of the domain
            "acoustic(1)%dir": -math.pi / 2,  # direction: -pi/2
            # Patch 1: Background
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": dlengx / 2,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": dlengx,
            "patch_icpp(1)%length_y": dlengy,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": patmos_n,
            "patch_icpp(1)%alpha_rho(1)": rhol_n * alphal_back,
            "patch_icpp(1)%alpha_rho(2)": rhog_n * alphag_back,
            "patch_icpp(1)%alpha(1)": alphal_back,
            "patch_icpp(1)%alpha(2)": alphag_back,
            # Patch 2: Lung
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%hcid": 205,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": dlengx / 2.0,
            "patch_icpp(2)%y_centroid": -dlengy / 4.0,
            "patch_icpp(2)%length_x": dlengx,
            "patch_icpp(2)%length_y": dlengy / 2.0 + 2,
            "patch_icpp(2)%a(2)": interface_amp,
            "patch_icpp(2)%vel(1)": 0.0e00,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": patmos_n,
            "patch_icpp(2)%alpha_rho(1)": rhol_n * alphal_lung,
            "patch_icpp(2)%alpha_rho(2)": rhog_n * alphag_lung,
            "patch_icpp(2)%alpha(1)": alphal_lung,
            "patch_icpp(2)%alpha(2)": alphag_lung,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gammal - 1.0e00),
            "fluid_pp(1)%pi_inf": gammal * Bl_n / (gammal - 1.0e00),
            "fluid_pp(2)%gamma": 1.0e00 / (gammag - 1.0e00),
            "fluid_pp(2)%pi_inf": gammag * Bg_n / (gammag - 1.0e00),
        }
    )
)
