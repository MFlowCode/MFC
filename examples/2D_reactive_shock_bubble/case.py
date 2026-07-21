#!/usr/bin/env python3
# 2D shock-bubble interaction (Haas & Sturtevant, JFM 1987). A planar shock in inert N2
# strikes a cylindrical bubble of light H2/O2 (~0.4x the ambient density). The misaligned
# pressure (across the shock) and density (across the bubble) gradients deposit baroclinic
# vorticity, winding the light bubble into a forward-jetting vortex pair -- the canonical
# Richtmyer-Meshkov roll-up. Multi-species (chemistry framework) so the bubble and ambient
# differ in composition and density; at the cool ambient temperature the mixture does not
# ignite, so this is the clean hydrodynamic roll-up (set --react T and a hot ambient for the
# reactive variant).
import argparse
import json
import sys

import cantera as ct

parser = argparse.ArgumentParser(prog="2D_reactive_shock_bubble")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--scale", type=float, default=1.0, help="Grid multiplier.")
parser.add_argument("--ndim", type=int, default=2, choices=(2, 3), help="Spatial dimensions (3 = spherical bubble -> vortex ring).")
parser.add_argument("--mach", type=float, default=2.5, help="Incident shock Mach number (in the N2 ambient).")
parser.add_argument("--ctfile", default="h2o2.yaml", help="Cantera mechanism (h2o2_xe.yaml adds inert Xe for the heavy reactive bubble).")
parser.add_argument("--xbub", default="H2:2,O2:1", help="Bubble mole-fraction composition. Literature deflagration: 'H2:0.30,O2:0.15,XE:0.55'.")
parser.add_argument("--pamb", type=float, default=101325.0, help="Ambient pressure [Pa]. Literature: ~0.25-0.5 atm.")
parser.add_argument("--tamb", type=float, default=300.0, help="Ambient temperature [K].")
parser.add_argument("--cfl", type=float, default=0.3, help="CFL number for the fixed timestep.")
parser.add_argument("--react", default="F", choices=("T", "F"), help="Enable reactions (default off: at cool ambient the light bubble rolls up without igniting).")
parser.add_argument("--substeps", type=int, default=4, help="alpha-QSS reaction sub-steps per flow step (floor when adaptive).")
parser.add_argument("--adap", default="F", choices=("T", "F"), help="Per-rank adaptive sub-step count (each rank sizes nsub from its own peak stiffness).")
parser.add_argument("--substeps-max", type=int, default=16, help="alpha-QSS sub-step ceiling when --adap T.")
parser.add_argument("--tend", type=float, default=1.0e-4, help="Physical end time [s].")
args = parser.parse_args()

ctfile = args.ctfile
Pamb = args.pamb
M = args.mach

# Inert N2 ambient and an H2/O2 bubble at the same p, T (composition sets its density).
ambient = ct.Solution(ctfile)
ambient.TPX = args.tamb, Pamb, "N2:1"
bubble = ct.Solution(ctfile)
bubble.TPX = args.tamb, Pamb, args.xbub

# Normal-shock (Rankine-Hugoniot) jump in the ambient, ideal gas with N2's cp/cv.
g = ambient.cp_mass / ambient.cv_mass
a1 = (g * Pamb / ambient.density) ** 0.5
p2 = Pamb * (1.0 + 2.0 * g / (g + 1.0) * (M**2 - 1.0))
rho_ratio = (g + 1.0) * M**2 / ((g - 1.0) * M**2 + 2.0)
rho2 = ambient.density * rho_ratio
u2 = M * a1 * (1.0 - 1.0 / rho_ratio)  # post-shock gas velocity (lab frame, +x)
a2 = a1 * (p2 / Pamb / rho_ratio) ** 0.5
print(
    f"ambient N2 rho={ambient.density:.4f} a1={a1:.0f} | bubble rho={bubble.density:.4f} "
    f"(ratio {bubble.density / ambient.density:.2f}) | shock M={M} u2={u2:.0f} T2/T1={p2 / Pamb / rho_ratio:.2f}",
    file=sys.stderr,
)

is_3d = args.ndim == 3
Ly = 0.03
Lz = Ly
if is_3d:
    Lx = 2.5 * Ly
    Ny = Nz = int(120 * args.scale)
    Nx = int(2.5 * Ny)
else:
    Lx = 3.0 * Ly
    Ny = int(360 * args.scale)
    Nx = 3 * Ny
    Nz = 0
dx = Lx / Nx
R = 0.007 if not is_3d else 0.25 * Ly  # bubble radius
x_bub = (0.30 if not is_3d else 0.35) * Lx  # upstream so the vortex has room to develop
x_shock = 0.15 * Lx  # initial incident-shock location, left of the bubble
geom_box = 9 if is_3d else 3
geom_bub = 8 if is_3d else 2

dt = args.cfl * dx / (u2 + a2)
NT = int(args.tend / dt)

case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": Ly,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": max(1, NT // 90),
    "parallel_io": "T",
    "model_eqns": "5eq",
    "num_fluids": 1,
    "num_patches": 3,
    "mpp_lim": "F",
    "mixture_err": "T",
    "weno_avg": "F",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "chemistry": "T",
    "chem_params%diffusion": "F",
    "chem_params%reactions": args.react,
    "chem_params%reaction_substeps": args.substeps,
    "chem_params%adap_substeps": args.adap,
    "chem_params%reaction_substeps_max": args.substeps_max,
    "cantera_file": ctfile,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "chem_wrt_T": "T",
    # Patch 1: inert N2 ambient at rest, whole domain
    "patch_icpp(1)%geometry": geom_box,
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%y_centroid": Ly / 2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": Pamb,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha_rho(1)": ambient.density,
    # Patch 2: post-shock N2 slab at the far left -> the incident shock
    "patch_icpp(2)%geometry": geom_box,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": x_shock / 2,
    "patch_icpp(2)%y_centroid": Ly / 2,
    "patch_icpp(2)%length_x": x_shock,
    "patch_icpp(2)%length_y": Ly,
    "patch_icpp(2)%vel(1)": u2,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p2,
    "patch_icpp(2)%alpha(1)": 1.0,
    "patch_icpp(2)%alpha_rho(1)": rho2,
    # Patch 3: light H2/O2 bubble at rest (circle in 2D, sphere in 3D)
    "patch_icpp(3)%geometry": geom_bub,
    "patch_icpp(3)%alter_patch(1)": "T",
    "patch_icpp(3)%x_centroid": x_bub,
    "patch_icpp(3)%y_centroid": Ly / 2,
    "patch_icpp(3)%radius": R,
    "patch_icpp(3)%vel(1)": 0.0,
    "patch_icpp(3)%vel(2)": 0.0,
    "patch_icpp(3)%pres": Pamb,
    "patch_icpp(3)%alpha(1)": 1.0,
    "patch_icpp(3)%alpha_rho(1)": bubble.density,
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}

for i in range(len(ambient.Y)):
    case[f"chem_wrt_Y({i + 1})"] = "T"
    case[f"patch_icpp(1)%Y({i + 1})"] = float(ambient.Y[i])  # N2
    case[f"patch_icpp(2)%Y({i + 1})"] = float(ambient.Y[i])  # shocked N2
    case[f"patch_icpp(3)%Y({i + 1})"] = float(bubble.Y[i])  # reactive bubble

if is_3d:
    case["p"] = Nz
    case.update(
        {
            "z_domain%beg": 0.0,
            "z_domain%end": Lz,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            "patch_icpp(1)%z_centroid": Lz / 2,
            "patch_icpp(1)%length_z": Lz,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(2)%z_centroid": Lz / 2,
            "patch_icpp(2)%length_z": Lz,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(3)%z_centroid": Lz / 2,
            "patch_icpp(3)%vel(3)": 0.0,
        }
    )

if __name__ == "__main__":
    print(json.dumps(case))
