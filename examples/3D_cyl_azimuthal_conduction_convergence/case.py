#!/usr/bin/env python3
# Azimuthal Fourier conduction verification: 3D cylindrical (x, r, theta), uniform pressure, zero velocity.
# In 3D cylindrical MFC the third coordinate is theta and dz is in RADIANS, so the azimuthal term
# carries two metric factors, (1/r**2) d2T/dtheta2. This case is the only one that exercises them.
#
# Ideal gas with T = p / ((Gamma - 1) * rho * cv), so setting
#     rho(r, theta) = RHO0 / (1 + A * (r/R)**2 * cos(theta))
# gives exactly T = T0 * (1 + A * (r/R)**2 * cos(theta)) with T0 = p / ((Gamma-1)*RHO0*cv).
# That profile is single-valued and smooth on the axis, and its cylindrical Laplacian is
#     (1/r) d/dr (r dT/dr) + (1/r**2) d2T/dtheta2 = 4*A*T0*cos(theta)/R**2 - A*T0*cos(theta)/R**2
#                                                 = 3*A*T0*cos(theta)/R**2,
# independent of r. The radial half is exact on the uniform 3D cylindrical r-grid (T is quadratic
# in r, so the two-point face gradient and the two-point average are both exact), which leaves the
# azimuthal second difference as the sole source of error: any mistake in the 1/r**2 metric shows up
# as an r-dependent residual rather than as a converging one.
import json
import os

GAM = 1.4
CV = 1.0
RHO0 = 1.0
P0 = 1.0
AMP = float(os.environ.get("AMP", "0.4"))
K_THERM = float(os.environ.get("K_THERM", "1.0"))
R = 2.0
LX = 1.0
NX = int(os.environ.get("NX", "32"))  # weno_order 5 needs m + 1 >= 26; T is uniform in x, so the x flux is zero
NR = int(os.environ.get("NR", "32"))
NP = int(os.environ.get("NP", "32"))
DT = 1.0e-8
TWO_PI = 2.0 * 3.141592653589793

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": LX,
            "y_domain%beg": 0.0,
            "y_domain%end": R,
            "z_domain%beg": 0.0,
            "z_domain%end": TWO_PI,
            "cyl_coord": "T",
            "m": NX - 1,
            "n": NR - 1,
            "p": NP - 1,
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
            "bc_y%beg": -14,  # axis
            "bc_y%end": -2,  # wall; only the outermost cell sees its mirrored ghost temperature
            "bc_z%beg": -1,  # periodic in theta: the azimuthal direction has no boundary error at all
            "bc_z%end": -1,
            "patch_icpp(1)%geometry": 10,  # cylinder about the x-axis, covering the whole domain
            "patch_icpp(1)%x_centroid": 0.5 * LX,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_x": LX,
            "patch_icpp(1)%radius": R,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": P0,
            "patch_icpp(1)%alpha_rho(1)": f"{RHO0} / (1.0 + {AMP} * (y * y / ({R} * {R})) * cos(z))",
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
