#!/usr/bin/env python3
import math
import json

# Single bubble collapse
# Description: A planar acoustic wave interacts with a single bubble
# in water. The background field is modeled in using an Eulerian framework,
# while the bubbles are tracked using a Lagrangian framework.

# Reference values for nondimensionalization
x0 = 1.0e-03  # length - m
rho0 = 1.0e03  # density - kg/m3
c0 = 1475.0  # speed of sound - m/s
T0 = 298  # temperature - K
p0 = rho0 * c0 * c0  # pressure - Pa

# Host properties (water + glicerol)
gamma_host = 2.7466  # Specific heat ratio
pi_inf_host = 792.02e06  # Stiffness - Pa
mu_host = 0.006  # Dynamic viscosity - Pa.s
c_host = 1475.0  # speed of sound - m/s
rho_host = 1000  # density kg/m3
T_host = T0  # temperature K

# Lagrangian bubble's properties
R_uni = 8314  # Universal gas constant - J/kmol/K
MW_g = 28.0  # Molar weight of the gas - kg/kmol
MW_v = 18.0  # Molar weight of the vapor - kg/kmol
gamma_g = 1.33  # Specific heat ratio of the gas
gamma_v = 1.33  # Specific heat ratio of the vapor
pv = 2500  # Vapor pressure of the host - Pa
cp_g = 1.0e3  # Specific heat of the gas - J/kg/K
cp_v = 2.1e3  # Specific heat of the vapor - J/kg/K
k_g = 0.025  # Thermal conductivity of the gas - W/m/K
k_v = 0.02  # Thermal conductivity of the vapor - W/m/K
diffVapor = 2.5e-5  # Diffusivity coefficient of the vapor - m2/s
sigBubble = 0.07  # Surface tension of the bubble - N/m
mu_gas = 1.48e-5

# Acoustic source properties
patm = 1.0e05  # Atmospheric pressure - Pa
pamp = 1.32e05  # Amplitude of the acoustic source - Pa
freq = 21.4e03  # Source frequency - Hz
wlen = c_host / freq  # Wavelength - m

# Domain and time set up

xb = -6.0e-3  # Domain boundaries - m (x direction)
xe = 6.0e-3
yb = -3.0e-3  # Domain boundaries - m (y direction)
ye = 3.0e-3
zb = -3.0e-3  # Domain boundaries - m (z direction)
ze = 3.0e-3

Nx = 2 * 60  # number of elements into x direction
Ny = 60  # number of elements into y direction
Nz = 60  # number of elements into z direction

dt = 4.0e-08  # time-step - sec
stopTime = 60.0e-06  # stop time - sec
saveTime = 30.0e-06  # save time - sec

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics ================================================
            "run_time_info": "T",
            # Computational Domain Parameters ==========================
            "x_domain%beg": xb / x0,
            "x_domain%end": xe / x0,
            "y_domain%beg": yb / x0,
            "y_domain%end": ye / x0,
            "z_domain%beg": zb / x0,
            "z_domain%end": ze / x0,
            "stretch_x": "F",
            "stretch_y": "F",
            "stretch_z": "F",
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "dt": round(dt * c0 / x0, 6),
            "n_start": 0,
            "t_save": saveTime * (c0 / x0),
            "t_stop": stopTime * (c0 / x0),
            # Simulation Algorithm Parameters ==========================
            "model_eqns": 2,
            "num_fluids": 2,
            "num_patches": 1,
            "mpp_lim": "F",
            "viscous": "T",
            "time_stepper": 4,  # 4th/5th RKCK
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -6,
            "bc_x%end": -6,
            "bc_y%beg": -6,
            "bc_y%end": -6,
            "bc_z%beg": -6,
            "bc_z%end": -6,
            # Formatted Database Files Structure Parameters ============
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: Water (left) ====================================
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_x": 2 * (xe - xb) / x0,
            "patch_icpp(1)%length_y": 2 * (ye - yb) / x0,
            "patch_icpp(1)%length_z": 2 * (ze - zb) / x0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": patm / p0,
            "patch_icpp(1)%alpha_rho(1)": rho_host / rho0,
            "patch_icpp(1)%alpha_rho(2)": 0.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%alpha(2)": 0.0,
            # Acoustic source ==========================================
            "acoustic_source": "T",
            "num_source": 1,
            "acoustic(1)%support": 3,
            "acoustic(1)%pulse": 1,
            "acoustic(1)%npulse": 10,
            "acoustic(1)%mag": -pamp / p0,
            "acoustic(1)%wavelength": wlen / x0,
            "acoustic(1)%length": 2 * (ze - zb) / x0,
            "acoustic(1)%height": 2 * (ye - yb) / x0,
            "acoustic(1)%loc(1)": -2.0e-03 / x0,
            "acoustic(1)%loc(2)": 0.0,
            "acoustic(1)%loc(3)": 0.0,
            "acoustic(1)%dir": 0.0,
            "acoustic(1)%delay": 0.0,
            # Lagrangian Bubbles ===========================
            "bubbles_lagrange": "T",
            "bubble_model": 2,  # Keller-Miksis model
            "rkck_adap_dt": "T",  # Activate adaptive time stepper
            "rkck_tolerance": 1.0e-05,
            "lag_params%nBubs_glb": 1,
            "lag_params%solver_approach": 2,  # Two-way coupled
            "lag_params%cluster_type": 2,
            "lag_params%pressure_corrector": "T",
            "lag_params%smooth_type": 1,
            "lag_params%heatTransfer_model": "T",
            "lag_params%massTransfer_model": "F",
            "lag_params%epsilonb": 1.0,
            "lag_params%valmaxvoid": 0.9,
            "lag_params%write_bubbles": "F",
            "lag_params%write_bubbles_stats": "F",
            "lag_params%c0": c0,
            "lag_params%rho0": rho0,
            "lag_params%T0": T0,
            "lag_params%x0": x0,
            "lag_params%diffcoefvap": diffVapor,
            "lag_params%Thost": T_host,
            # Fluids Physical Parameters ===============================
            # Host medium
            "fluid_pp(1)%gamma": 1.0 / (gamma_host - 1.0),
            "fluid_pp(1)%pi_inf": gamma_host * (pi_inf_host / p0) / (gamma_host - 1.0),
            "fluid_pp(1)%Re(1)": 1.0 / (mu_host / (rho0 * c0 * x0)),
            "fluid_pp(1)%mul0": mu_host,
            "fluid_pp(1)%ss": sigBubble,
            "fluid_pp(1)%pv": pv,
            "fluid_pp(1)%gamma_v": gamma_v,
            "fluid_pp(1)%M_v": MW_v,
            "fluid_pp(1)%k_v": k_v,
            "fluid_pp(1)%cp_v": cp_v,
            # Bubble gas state
            "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
            "fluid_pp(2)%pi_inf": 0.0e00,
            "fluid_pp(2)%Re(1)": 1.0 / (mu_gas / (rho0 * c0 * x0)),
            "fluid_pp(2)%gamma_v": gamma_g,
            "fluid_pp(2)%M_v": MW_g,
            "fluid_pp(2)%k_v": k_g,
            "fluid_pp(2)%cp_v": cp_g,
        }
    )
)
