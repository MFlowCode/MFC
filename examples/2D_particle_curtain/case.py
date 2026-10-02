#!/usr/bin/env python3
# Planar shock (driven by a Dirichlet inflow) hitting a dense particle curtain between slip walls
# (Euler-Lagrange, two-way coupling). The particles come from gen_particles.py.
import argparse
import json
import math
import os

import gen_particles

parser = argparse.ArgumentParser(prog="2D_particle_curtain", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
args = parser.parse_args()

# Write input/particles.dat only when ./mfc.sh run runs pre_process (validate, test and docs also load this file)
if args.mfc.get("command") == "run" and "pre_process" in args.mfc.get("targets", []):
    gen_particles.write_particles(os.path.dirname(os.path.abspath(__file__)))

# Air (ideal gas)
gamma = 1.4
R_air = 287.0  # J/(kg K)


def normal_shock(M, p1, rho1):
    """Pressure, density and lab-frame velocity behind a normal shock of Mach number M moving into still gas (p1, rho1)."""
    p2 = p1 * (1.0 + 2.0 * gamma / (gamma + 1.0) * (M**2 - 1.0))
    rho2 = rho1 * (gamma + 1.0) * M**2 / ((gamma - 1.0) * M**2 + 2.0)
    u2 = 2.0 / (gamma + 1.0) * math.sqrt(gamma * p1 / rho1) * (M - 1.0 / M)
    return p2, rho2, u2


# Ambient air (Wagner et al. 2012) and the state behind the incident shock
M_shock = 1.66
p_amb, T_amb = 82700.0, 296.4
rho_amb = p_amb / (R_air * T_amb)
p_post, rho_post, u_post = normal_shock(M_shock, p_amb, rho_amb)

# Grid: x in [-0.2, 0.3], y in [0, 0.015] (slip walls at y = 0 and y = 0.015)
xb, xe, yb, ye = -0.2, 0.3, 0.0, 0.015
Nx, Ny = 2000, 30
dx = (xe - xb) / Nx

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational domain
            "x_domain%beg": xb,
            "x_domain%end": xe,
            "y_domain%beg": yb,
            "y_domain%end": ye,
            "m": Nx - 1,
            "n": Ny - 1,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 0.4,
            "n_start": 0,
            "t_stop": 5.0e-4,
            "t_save": 5.0e-5,
            # Simulation algorithm
            "model_eqns": 2,
            "num_fluids": 1,
            "num_patches": 2,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -17,
            "bc_x%end": -8,
            "bc_y%beg": -15,
            "bc_y%end": -15,
            # Output
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "lag_db_wrt": "T",
            "lag_voidfrac_wrt": "T",
            # Patch 1: ambient air
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * (xb + xe),
            "patch_icpp(1)%y_centroid": 0.5 * (yb + ye),
            "patch_icpp(1)%length_x": xe - xb,
            "patch_icpp(1)%length_y": ye - yb,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p_amb,
            "patch_icpp(1)%alpha_rho(1)": rho_amb,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: post-shock state on the left (x < -5 mm); the Dirichlet inflow holds it
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -0.1025,
            "patch_icpp(2)%y_centroid": 0.5 * (yb + ye),
            "patch_icpp(2)%length_x": 0.195,
            "patch_icpp(2)%length_y": ye - yb,
            "patch_icpp(2)%vel(1)": u_post,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": p_post,
            "patch_icpp(2)%alpha_rho(1)": rho_post,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluid: air
            "fluid_pp(1)%eos": "ideal_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%cv": 717.5,
            # Lagrangian particles
            "particles_lagrange": "T",
            "fd_order": 2,
            "particle_pp%rho0ref_particle": 2520.0,
            "particle_params%input_path": "input/particles.dat",
            "particle_params%nparticles_glb": gen_particles.n_particles,
            "particle_params%solver_approach": 2,
            "particle_params%qs_force": 3,
            "particle_params%pressure_gradient_force": "T",
            "particle_params%added_mass_force": 1,
            "particle_params%mu_ref(1)": 1.716e-5,
            "particle_params%suth(1)": 110.4,
            "particle_params%interpolation_order": 2,
            "particle_params%epsilonb": 1.0,
            "particle_params%valmaxvoid": 0.9,
            "particle_params%charwidth": gen_particles.charwidth,
        },
        indent=4,
    )
)
