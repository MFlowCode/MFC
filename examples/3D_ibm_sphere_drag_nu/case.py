#!/usr/bin/env python3
"""
Sphere in crossflow at Re_D = 20-300, Mach 0.2: drag and Nusselt number.

A no-slip isothermal sphere (immersed boundary, thermal_bc = 1) in dry air at Mach 0.2. Chemistry is on for the gas
thermodynamics and transport only (H2/O2 mechanism, no fuel, reactions off). D follows from --Re at fixed p_inf and T_inf.
--tw_ratio 1 gives a drag run with the wall at the free-stream temperature; --tw_ratio 1.1 a heated wall for the Nusselt
number. The force goes to D/ib_forces.dat (ib_force_wrt, from the face fluxes) and the heat the body takes to D/ib_heat.dat.

    ./mfc.sh run examples/3D_ibm_sphere_drag_nu/case.py -n 1 --gpu -- --Re 100 --cells_per_D 20 --xu 4 --xd 8 --hy 4 --steps 0 --tstar 40
    python3 examples/3D_ibm_sphere_drag_nu/check_sphere.py <case directory> [--hot <heated run directory>]

Defaults are a CI-sized run (Re 20, 6 cells per D, a 4.2 D cube, 30 fixed steps: the flow is not steady and the body is a
smoke-level 6 cells across; the chemistry transport makes 3D cells expensive). Table runs use --cells_per_D 20, --xu 4 --xd 8 --hy 4 (blockage 1.2 %) and the adaptive
time step (--steps 0 with --tstar); the README lists the commands.
"""

import argparse
import json
import math
import sys

import cantera as ct

parser = argparse.ArgumentParser(prog="3D_ibm_sphere_drag_nu")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--Re", type=float, default=20.0, help="Re_D = rho U D / mu at the free stream.")
parser.add_argument("--mach", type=float, default=0.2, help="Free-stream Mach number.")
parser.add_argument("--tw_ratio", type=float, default=1.0, help="Wall temperature over free-stream temperature.")
parser.add_argument("--cells_per_D", type=float, default=6.0, help="Cells across the diameter.")
parser.add_argument("--xu", type=float, default=1.6, help="Upstream extent [D].")
parser.add_argument("--xd", type=float, default=2.6, help="Downstream extent [D].")
parser.add_argument("--hy", type=float, default=2.1, help="Lateral half width in y and z [D]; the area blockage is pi/(16 hy^2).")
parser.add_argument("--steps", type=int, default=30, help="Fixed number of steps (CI-sized run); 0 selects the adaptive step with --tstar.")
parser.add_argument("--tstar", type=float, default=50.0, help="Run length t U / D when --steps is 0.")
parser.add_argument("--cfl", type=float, default=0.4, help="CFL number of the fixed step (the adaptive step uses cfl_target, set from Re and --cells_per_D).")
parser.add_argument("--frames", type=int, default=10, help="Saves of the adaptive run (a fixed-step run saves once, at the end).")
parser.add_argument("--stride", type=int, default=20, help="ib_force_stride: steps between force and heat rows.")
parser.add_argument("--ib_ip_interp", type=int, default=1, choices=(1, 2), help="Image-point interpolation (1: inverse distance, the default).")
parser.add_argument("--silo", action="store_true", help="Write density, pressure and velocity for post_process (the momentum-balance check).")
parser.add_argument("--serial_io", action="store_true", help="parallel_io = F (a build without MPI); the checker then reads p_all/.")
args = parser.parse_args()

T_inf, p_inf = 300.0, 101325.0
gas = ct.Solution("h2o2.yaml")
gas.TPX = T_inf, p_inf, "O2:0.21,N2:0.79"
rho, mu, a = gas.density, gas.viscosity, gas.sound_speed
u_inf = args.mach * a
D = args.Re * mu / (rho * u_inf)
dx = D / args.cells_per_D
Lx, Ly = (args.xu + args.xd) * D, 2.0 * args.hy * D
m, n = round(Lx / dx) - 1, round(Ly / dx) - 1
p = n
steps = args.steps

case = {
    "run_time_info": "T",
    "x_domain%beg": -args.xu * D,
    "x_domain%end": args.xd * D,
    "y_domain%beg": -args.hy * D,
    "y_domain%end": args.hy * D,
    "z_domain%beg": -args.hy * D,
    "z_domain%end": args.hy * D,
    "m": m,
    "n": n,
    "p": p,
    "cyl_coord": "F",
    "num_patches": 1,
    "num_fluids": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": 2,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "viscous": "T",
    "fd_order": 4,
    # Chemistry: thermodynamics and mixture-averaged transport only
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "F",
    "chem_params%transport_model": 1,
    "cantera_file": "h2o2.yaml",
    "chem_wrt_T": "T",
    # Free stream
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.5 * (args.xd - args.xu) * D,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_z": Ly,
    "patch_icpp(1)%vel(1)": u_inf,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p_inf,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0,
    # Subsonic characteristic inflow, outflow holding p_inf (without it the far-field pressure creeps up), extrapolated sides
    "bc_x%beg": -7,
    "bc_x%grcbc_in": "T",
    "bc_x%vel_in(1)": u_inf,
    "bc_x%vel_in(2)": 0.0,
    "bc_x%vel_in(3)": 0.0,
    "bc_x%pres_in": p_inf,
    "bc_x%alpha_rho_in(1)": rho,
    "bc_x%alpha_in(1)": 1.0,
    "bc_x%end": -8,
    "bc_x%grcbc_out": "T",
    "bc_x%pres_out": p_inf,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "bc_z%beg": -3,
    "bc_z%end": -3,
    # No-slip isothermal sphere, off the symmetry planes by fractions of dx
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 8,
    "patch_ib(1)%x_centroid": dx / 3.0,
    "patch_ib(1)%y_centroid": dx / 5.0,
    "patch_ib(1)%z_centroid": dx / 7.0,
    "patch_ib(1)%radius": 0.5 * D,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%thermal_bc": 1,
    "patch_ib(1)%Twall": args.tw_ratio * T_inf,
    "ib_ip_interp": args.ib_ip_interp,
    "ib_state_wrt": "T",
    "ib_force_wrt": "T",
    "ib_heat_wrt": "T",
    "ib_force_stride": args.stride,
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    "parallel_io": "F" if args.serial_io else "T",
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "F",
}
for i, y in enumerate(gas.Y):
    case[f"patch_icpp(1)%Y({i + 1})"] = float(y)
    case[f"chem_wrt_Y({i + 1})"] = "F"
if args.silo:
    case.update({"rho_wrt": "T", "pres_wrt": "T", "vel_wrt(1)": "T", "vel_wrt(2)": "T", "vel_wrt(3)": "T"})

if steps > 0:
    # Fixed step: acoustic CFL, and the viscous limit nu dt / dx^2 stays well below 0.1 at these Re
    dt = args.cfl * dx / (u_inf + a)
    case.update({"dt": dt, "t_step_start": 0, "t_step_stop": steps, "t_step_save": steps})
else:
    t_stop = args.tstar * D / u_inf
    case.update(
        {
            "cfl_adap_dt": "T",
            "cfl_target": min(0.5, 0.4 * args.Re / args.cells_per_D),
            "chem_params%transport_dt": "T",
            "n_start": 0,
            "t_save": t_stop / args.frames,
            "t_stop": t_stop,
        }
    )

if __name__ == "__main__":
    print(
        f"Re_D {args.Re:g}  Mach {args.mach:g}  D = {D:.4e} m  U = {u_inf:.2f} m/s  grid {m + 1} x {n + 1} x {p + 1}  "
        f"blockage {math.pi / (16.0 * args.hy**2):.4f}  Twall = {args.tw_ratio * T_inf:.1f} K",
        file=sys.stderr,
    )
    print(json.dumps(case))
