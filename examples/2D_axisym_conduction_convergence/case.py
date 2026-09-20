#!/usr/bin/env python3
# Axisymmetric Fourier conduction verification: uniform pressure, zero velocity, cylindrical (x, r).
# Ideal gas with T = p / ((Gamma - 1) * rho * cv), so setting
#     rho(r) = RHO0 / (1 + A*(1 - r**2/R**2) + B*(1 - r**4/R**4))
# gives exactly T(r) = T0 * (1 + A*(1 - r**2/R**2) + B*(1 - r**4/R**4)) with
# T0 = p / ((Gamma-1)*RHO0*cv), whose cylindrical Laplacian is
#     (1/r) d/dr (r dT/dr) = -4*A*T0/R**2 - 16*B*T0*r**2/R**4,
# finite on the axis. At B = 0 (the default) that is a constant, so every cell including the
# axis cell must return the same number: a wrong sign or factor in the k = 0 axis source stands
# out against a flat field, where a sinusoid would hide it. B > 0 makes the answer vary with r
# and turns the same case into an ordinary second-order convergence test.
import json
import os

NR = int(os.environ.get("NR", "100"))

GAM = 1.4
CV = 1.0
RHO0 = 1.0
P0 = 1.0
AMP = float(os.environ.get("AMP", "0.4"))
BAMP = float(os.environ.get("BAMP", "0.0"))
K_THERM = float(os.environ.get("K_THERM", "1.0e-2"))
R = 1.0
LX = 1.0
NX = 30
DT = 1.0e-8

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": LX,
            "y_domain%beg": 0.0,
            "y_domain%end": R,
            "cyl_coord": "T",
            "m": NX - 1,
            "n": NR - 1,
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
            "bc_x%beg": -1,  # periodic: T has no x variation, so the x flux is identically zero
            "bc_x%end": -1,
            "bc_y%beg": -2,  # axis
            "bc_y%end": -2,  # wall; only the outermost cells see its mirrored ghost temperature
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * LX,
            "patch_icpp(1)%y_centroid": 0.5 * R,
            "patch_icpp(1)%length_x": LX,
            "patch_icpp(1)%length_y": R,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": P0,
            "patch_icpp(1)%alpha_rho(1)": (f"{RHO0} / (1.0 + {AMP} * (1.0 - y * y / ({R} * {R}))" f" + {BAMP} * (1.0 - y * y * y * y / ({R} * {R} * {R} * {R})))"),
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
