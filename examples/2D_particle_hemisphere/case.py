#!/usr/bin/env python3
# Mach 10 cylindrical blast wave hitting a half-ring of solid particles (Euler-Lagrange, two-way coupling).
# The particles come from gen_particles.py.
import argparse
import json
import math
import os

import gen_particles

parser = argparse.ArgumentParser(prog="2D_particle_hemisphere", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
args = parser.parse_args()

# Write input/particles.dat only when ./mfc.sh run runs pre_process (validate, test and docs also load this file)
if args.mfc.get("command") == "run" and "pre_process" in args.mfc.get("targets", []):
    gen_particles.write_particles(os.path.dirname(os.path.abspath(__file__)))

# Air (ideal gas)
gamma = 1.4
R_air = 287.05  # J/(kg K)


def normal_shock(M, p1, rho1):
    """Pressure, density and lab-frame velocity behind a normal shock of Mach number M moving into still gas (p1, rho1)."""
    p2 = p1 * (1.0 + 2.0 * gamma / (gamma + 1.0) * (M**2 - 1.0))
    rho2 = rho1 * (gamma + 1.0) * M**2 / ((gamma - 1.0) * M**2 + 2.0)
    u2 = 2.0 / (gamma + 1.0) * math.sqrt(gamma * p1 / rho1) * (M - 1.0 / M)
    return p2, rho2, u2


# Ambient air, and a driver at rest at the pressure and density behind a Mach M_blast shock
M_blast = 10.0
p_amb, T_amb = 82.7e3, 297.0
rho_amb = p_amb / (R_air * T_amb)
p_blast, rho_blast, _ = normal_shock(M_blast, p_amb, rho_amb)
r_blast = 0.0115  # m

# Grid: x in [-0.1, 0.1], y in [0, 0.1] (slip wall at y = 0)
dx = 250e-6
xb, xe, yb, ye = -0.1, 0.1, 0.0, 0.1
Nx, Ny = round((xe - xb) / dx), round((ye - yb) / dx)

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
            "t_stop": 5.0e-5,
            "t_save": 5.0e-6,
            # Simulation algorithm
            "model_eqns": 2,
            "num_fluids": 1,
            "num_patches": 2,
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "riemann_solver": 1,  # HLL: HLLC shows a carbuncle on the blast front along the x = 0 grid line
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -8,
            "bc_x%end": -8,
            "bc_y%beg": -15,
            "bc_y%end": -8,
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
            # Patch 2: high-pressure blast driver at the origin
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%radius": r_blast,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": p_blast,
            "patch_icpp(2)%alpha_rho(1)": rho_blast,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluid: air
            "fluid_pp(1)%eos": "ideal_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%cv": 717.5,
            # Lagrangian particles
            "particles_lagrange": "T",
            "fd_order": 2,
            "particle_pp%rho0ref_particle": 1738.0,
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
            "particle_params%charwidth": dx,
        },
        indent=4,
    )
)
