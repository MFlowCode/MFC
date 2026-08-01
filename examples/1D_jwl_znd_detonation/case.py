#!/usr/bin/env python3
import json
import math

# LX-10-0 JWL products parameters from Garno et al., J. Appl. Phys. 128,
# 195903 (2020), Table I. The case passes Q (J/kg) and MFC derives
# jwl_E0 = rho0*Q internally.
jwl_A = 880.2e9
jwl_B = 17.437e9
jwl_R1 = 4.60
jwl_R2 = 1.20
jwl_omega = 0.30
jwl_rho0 = 1860.0
jwl_Cv = 1000.0
jwl_Q = 5.59e6

# Reactant/product energy offset (Garno Eq. 17): chosen so unreacted
# (lambda = 0) explosive at (rho0, e = 0) sits at ambient pressure, which
# places it on a stiffer Hugoniot than the products and gives the reactive
# burn a resolved ZND structure (von Neumann spike -> CJ state).
p_amb = 101325.0
F1 = jwl_A * (1.0 - jwl_omega / jwl_R1) * math.exp(-jwl_R1) + jwl_B * (1.0 - jwl_omega / jwl_R2) * math.exp(-jwl_R2)
jwl_delta_e = (p_amb - F1) / (jwl_omega * jwl_rho0)

# JWL++ rate constants: dlambda/dt = jwl_G * p^jwl_b_exp * (1 - lambda).
# G is set for a ~0.5 mm reaction zone at the ~40 GPa post-spike pressures
# of this material (e-fold time ~100 ns), so the zone spans ~50 cells here.
jwl_G = 6.0e-15
jwl_b_exp = 2.0

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": 0.02,
            "m": 1999,
            "n": 0,
            "p": 0,
            "dt": 2.0e-10,
            "t_step_start": 0,
            "t_step_stop": 9500,
            "t_step_save": 950,
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "recon_type": 2,
            "muscl_order": 2,
            "muscl_lim": 2,
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: high-pressure hot spot igniting the unreacted charge
            # (jwl_reactive initializes lambda = 0 everywhere).
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.0005,
            "patch_icpp(1)%length_x": 0.001,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": 45.0e9,
            "patch_icpp(1)%alpha_rho(1)": jwl_rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: unreacted explosive at rest, ambient pressure.
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.0105,
            "patch_icpp(2)%length_x": 0.019,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": 2.0e5,
            "patch_icpp(2)%alpha_rho(1)": jwl_rho0,
            "patch_icpp(2)%alpha(1)": 1.0,
            "fluid_pp(1)%eos": 2,
            "fluid_pp(1)%gamma": 1.0 / 0.4,
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%cv": jwl_Cv,
            "fluid_pp(1)%jwl_A": jwl_A,
            "fluid_pp(1)%jwl_B": jwl_B,
            "fluid_pp(1)%jwl_R1": jwl_R1,
            "fluid_pp(1)%jwl_R2": jwl_R2,
            "fluid_pp(1)%jwl_omega": jwl_omega,
            "fluid_pp(1)%jwl_rho0": jwl_rho0,
            "fluid_pp(1)%jwl_Q": jwl_Q,
            "fluid_pp(1)%jwl_air_e0": 2.5575e5,
            "fluid_pp(1)%jwl_air_rho0": 1.225,
            "fluid_pp(1)%jwl_delta_e": jwl_delta_e,
            "jwl_reactive": "T",
            "jwl_G": jwl_G,
            "jwl_b_exp": jwl_b_exp,
        }
    )
)
