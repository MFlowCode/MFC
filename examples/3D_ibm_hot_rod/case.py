#!/usr/bin/env python3
# Mach 1.5 air crossflow past a finite 3D cylinder (rod axis along x, flat end caps) with a hot isothermal wall.
# Chemistry is on for its thermodynamics (H2/O2 mechanism, no fuel). Without face_bound = 1 the WENO face states next to
# the hot caps overshoot past T_surface_max, where the thermodynamic fit gives a negative cp and a NaN sound speed.
import argparse
import json

import cantera as ct

parser = argparse.ArgumentParser(prog="3D_ibm_hot_rod")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--cells_per_D", type=int, default=8, help="Cells across the rod diameter.")
parser.add_argument("--steps", type=int, default=300, help="Number of time steps.")
parser.add_argument("--face_bound", type=int, default=1, choices=(0, 1), help="Positivity bound on face states.")
args = parser.parse_args()

D = 1.0e-3  # rod diameter [m]
Mach, Twall = 1.5, 2500.0

gas = ct.Solution("h2o2.yaml")
gas.TPX = 300.0, 101325.0, "O2:0.21,N2:0.79"
u_inf = Mach * gas.sound_speed

dx = D / args.cells_per_D
Lx, Ly = 5 * D, 4 * D  # the rod (length 3 D) spans x in [-1.5 D, 1.5 D]
m, n = round(Lx / dx) - 1, round(Ly / dx) - 1
dt = 0.1 * dx / (u_inf + gas.sound_speed)

case = {
    "run_time_info": "T",
    "x_domain%beg": -Lx / 2,
    "x_domain%end": Lx / 2,
    "y_domain%beg": -Ly / 2,
    "y_domain%end": Ly / 2,
    "z_domain%beg": -Ly / 2,
    "z_domain%end": Ly / 2,
    "m": m,
    "n": n,
    "p": n,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": args.steps,
    "t_step_save": args.steps,
    "parallel_io": "T",
    # Algorithm
    "num_patches": 1,
    "model_eqns": 2,
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "weno_avg": "T",
    "weno_Re_flux": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "fd_order": 2,
    "viscous": "T",
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "bc_z%beg": -3,
    "bc_z%end": -3,
    "face_bound": args.face_bound,
    # Chemistry (thermodynamics and diffusion; no fuel, so no reactions)
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "F",
    "cantera_file": "h2o2.yaml",
    "chem_wrt_T": "T",
    # Hot no-slip rod with flat caps
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 10,
    "patch_ib(1)%x_centroid": 0.0,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%z_centroid": 0.0,
    "patch_ib(1)%radius": D / 2,
    "patch_ib(1)%length_x": 3 * D,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%thermal_bc": 1,
    "patch_ib(1)%Twall": Twall,
    # Cold crossflow in y
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%length_z": Ly,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": u_inf,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": gas.P,
    "patch_icpp(1)%alpha_rho(1)": gas.density,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / gas.viscosity,
    # Output
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
}

for i in range(len(gas.Y)):
    case[f"patch_icpp(1)%Y({i + 1})"] = float(gas.Y[i])
    case[f"chem_wrt_Y({i + 1})"] = "T"

if __name__ == "__main__":
    print(json.dumps(case))
