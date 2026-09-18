#!/usr/bin/env python3
# 1D Fourier conduction verification: uniform pressure, zero velocity, periodic.
# Ideal gas with T = p / ((Gamma - 1) * rho * cv), so setting
#     rho(x) = RHO0 / (1 + A*sin(2*pi*x/L))
# gives exactly T(x) = T0 * (1 + A*sin(2*pi*x/L)) with T0 = p / ((Gamma-1)*RHO0*cv).
# At t = 0 velocity is zero and pressure is uniform, so every Euler flux vanishes and
#     d(rho*E)/dt = k * d2T/dx2 = -k * T0 * A * (2*pi/L)**2 * sin(2*pi*x/L)
# exactly. One time step therefore measures the conduction term in isolation.
import json
import os

NX = int(os.environ.get("NX", "100"))

GAM = 1.4
CV = 1.0
RHO0 = 1.0
P0 = 1.0
AMP = 0.1
K_THERM = 1.0e-3
L = 1.0
DT = 1.0e-8

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "m": NX - 1,
            "n": 0,
            "p": 0,
            "dt": DT,
            "t_step_start": 0,
            "t_step_stop": 1,
            "t_step_save": 1,
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,  # periodic
            "bc_x%end": -1,
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": P0,
            "patch_icpp(1)%alpha_rho(1)": f"{RHO0} / (1.0 + {AMP} * sin(2.0 * pi * x / {L}))",
            "patch_icpp(1)%alpha(1)": 1.0,
            "fluid_pp(1)%eos": "ideal_gas",
            "fluid_pp(1)%gamma": 1.0 / (GAM - 1.0),
            "fluid_pp(1)%cv": CV,
            "fluid_pp(1)%k_therm": K_THERM,
            "parallel_io": "T",
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
        }
    )
)
