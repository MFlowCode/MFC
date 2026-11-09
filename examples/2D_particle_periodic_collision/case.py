#!/usr/bin/env python3
# Four identical head-on collisions of dense Euler-Lagrange particles (four-way coupling) in a doubly periodic box: one in the
# middle, one across the x boundary, one across the y boundary, and one across the corner (top left with bottom right).
# The particles are dense enough that the still fluid barely acts on them, so every pair should feel the same contact force.
import argparse
import json
import math
import os

parser = argparse.ArgumentParser(prog="2D_particle_periodic_collision", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
args = parser.parse_args()

radius = 0.005  # diameter below the cell size, as the one-cell collision halo requires
particle_density = 1.0e8  # relative to the fluid; the contact acceleration does not depend on it, the fluid force falls as 1/density
speed = 0.1  # of each particle toward its partner
gap = 0.010146  # between the surfaces at t = 0; contact at t = 0.05073 avoids the RK stage times t, t + dt/2 and t + dt
collision_time = 5.0e-3
dt = collision_time / 50.0
Nx = Ny = 64

# (contact point, unit normal from particle 1 to particle 2); the box is [0, 1]^2, so a contact on a face or corner wraps
pairs = [
    ((0.5, 0.5), (1.0, 0.0)),  # middle, along x
    ((0.0, 0.75), (-1.0, 0.0)),  # across the x boundary
    ((0.25, 0.0), (0.0, -1.0)),  # across the y boundary, the x pair rotated 90 degrees
    ((0.0, 1.0), (-1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0))),  # top left with bottom right, across the corner
]


def particle_rows():
    """Two rows (x, y, z, u, v, w, radius) per pair: particle 1 moves along the normal, particle 2 against it."""
    offset = radius + 0.5 * gap
    rows = []
    for (cx, cy), (nx, ny) in pairs:
        for sign in (-1.0, 1.0):
            x = (cx + sign * offset * nx) % 1.0
            y = (cy + sign * offset * ny) % 1.0
            rows.append((x, y, 0.0, -sign * speed * nx, -sign * speed * ny, 0.0, radius))
    return rows


# Write input/particles.dat only when ./mfc.sh run runs pre_process
if args.mfc.get("command") == "run" and "pre_process" in args.mfc.get("targets", []):
    case_dir = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(os.path.join(case_dir, "input"), exist_ok=True)
    with open(os.path.join(case_dir, "input", "particles.dat"), "w") as f:
        for row in particle_rows():
            f.write(" ".join(f"{value:.16e}" for value in row) + "\n")

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "F",
            # Computational domain
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "m": Nx - 1,
            "n": Ny - 1,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": 800,
            "t_step_save": 800,
            # Simulation algorithm
            "model_eqns": 2,
            "num_fluids": 1,
            "num_patches": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            # Output
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: still fluid
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            # Lagrangian particles: no drag or added mass, and dense, so only contacts change their velocity
            "particles_lagrange": "T",
            "fd_order": 2,
            "particle_pp%rho0ref_particle": particle_density,
            "particle_params%input_path": "input/particles.dat",
            "particle_params%nparticles_glb": 2 * len(pairs),
            "particle_params%coupling": "four_way",
            "particle_params%vel_model": 2,
            "particle_params%qs_force": 0,
            "particle_params%pressure_gradient_force": "T",
            "particle_params%added_mass_force": 0,
            "particle_params%interpolation_order": 2,
            "particle_params%valmaxvoid": 0.9,
            "particle_params%charwidth": 2.0 * radius,
            "particle_params%write_particles": "T",
            # Soft-sphere contacts
            "collision_model": 1,
            "coefficient_of_restitution": 0.8,
            "collision_time": collision_time,
            "coefficient_of_friction": 0.0,
        },
        indent=4,
    )
)
