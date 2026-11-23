#!/usr/bin/env python3
"""
Inviscid gas around a fixed slip circular cylinder (immersed boundary, slip = T): mass budget and exact vortex.

--mode leak    uniform flow at Mach --mach in a fully periodic box: the gas mass must stay constant. This is the CI-sized
               default; --steps 0 runs --tr flow-through times of the box (the table of the README). --mass_relax 1 returns the
               gas mass the slip ghost cells take in the step it occurs; 0 leaves it.
--mode vortex  the steady isentropic vortex u_theta = M c0 R / r (stagnation state p0, rho0 at infinity) is an exact Euler
               solution tangent to the wall: the near-wall error of the ghost state is measured against it. The initial
               condition is an analytic patch, so this mode compiles its own pre_process (about as long as a full build): it is
               not in the test suite.

    ./mfc.sh run examples/2D_ibm_slip_cylinder/case.py -t pre_process simulation post_process -n 1 -- --steps 0 --cells_per_R 16
    python3 examples/2D_ibm_slip_cylinder/check_slip.py examples/2D_ibm_slip_cylinder --case-args '--steps 0 --cells_per_R 16'

Physical values of the tables: R = 50 um, rho = 1.2 kg/m^3 and p = 1e5 Pa (stagnation state), gamma 1.4, M 0.5 (leak) or 0.3
and 0.6 (vortex).
"""

import argparse
import json

parser = argparse.ArgumentParser(prog="2D_ibm_slip_cylinder", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--mode", default="leak", choices=("leak", "vortex"))
parser.add_argument("--cells_per_R", type=int, default=8, help="Cells per cylinder radius.")
parser.add_argument("--half", type=float, default=6.0, help="Half width of the box [R] (vortex mode: 3).")
parser.add_argument("--mach", type=float, default=0.5, help="Free-stream Mach number (leak) or Mach number at the wall (vortex).")
parser.add_argument("--steps", type=int, default=100, help="Fixed steps (CI-sized run); 0 runs --tr flow-through times (leak) or 0.4 half widths of sound travel (vortex) at CFL 0.5.")
parser.add_argument("--tr", type=float, default=1.0, help="leak mode, --steps 0: flow-through times of the box.")
parser.add_argument("--bc", type=int, default=-1, help="leak mode: boundary condition of all four faces (-1 periodic, -3 extrapolation).")
parser.add_argument("--mass_relax", type=int, default=-1, help="ib_mass_relax (-1 auto: off for a fixed body).")
parser.add_argument("--ccst", action="store_true", help="ib_slip_ccst = T.")
parser.add_argument("--ip_interp", type=int, default=1, choices=(1, 2), help="ib_ip_interp (2: TSC image-point weights).")
args = parser.parse_args()

R = 50.0e-6
gam, rho0, p0 = 1.4, 1.2, 1.0e5
c0 = (gam * p0 / rho0) ** 0.5
half = 3.0 if args.mode == "vortex" else args.half
nhalf = round(half * args.cells_per_R)
dx = R / args.cells_per_R
u_ref = args.mach * c0
vortex = args.mode == "vortex"

case = {
    "run_time_info": "T",
    "x_domain%beg": -nhalf * dx,
    "x_domain%end": nhalf * dx,
    "y_domain%beg": -nhalf * dx,
    "y_domain%end": nhalf * dx,
    "m": 2 * nhalf - 1,
    "n": 2 * nhalf - 1,
    "p": 0,
    "num_patches": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "viscous": "F",
    "fd_order": 2,
    "bc_x%beg": -3 if vortex else args.bc,
    "bc_x%end": -3 if vortex else args.bc,
    "bc_y%beg": -3 if vortex else args.bc,
    "bc_y%end": -3 if vortex else args.bc,
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": 0.0,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%radius": R,
    "patch_ib(1)%slip": "T",
    "patch_ib(1)%moving_ibm": 0,
    "ib_mass_relax": args.mass_relax,
    "ib_slip_ccst": "T" if args.ccst else "F",
    "ib_ip_interp": args.ip_interp,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": 2 * nhalf * dx,
    "patch_icpp(1)%length_y": 2 * nhalf * dx,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "F",
    "rho_wrt": "T",
    "pres_wrt": "T",
    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",
    "parallel_io": "F",
}

if vortex:
    # u_theta = M c0 R / max(r, R); T/T0 = 1 - (gamma - 1)/2 (u_theta/c0)^2; isentropic rho and p
    r = "max(sqrt(x**2 + y**2), 1e-30)"
    ut = f"({args.mach!r}*{c0!r}*{R!r}/max({r}, {R!r}))"
    tr = f"(1 - {0.5 * (gam - 1.0)!r}*({ut}/{c0!r})**2)"
    case.update(
        {
            "patch_icpp(1)%vel(1)": f"-{ut}*y/{r}",
            "patch_icpp(1)%vel(2)": f"{ut}*x/{r}",
            "patch_icpp(1)%pres": f"{p0!r}*{tr}**{gam / (gam - 1.0)!r}",
            "patch_icpp(1)%alpha_rho(1)": f"{rho0!r}*{tr}**{1.0 / (gam - 1.0)!r}",
        }
    )
    t_end = 0.4 * half * R / c0
    cfl_run = True
else:
    case.update({"patch_icpp(1)%vel(1)": u_ref, "patch_icpp(1)%vel(2)": 0.0, "patch_icpp(1)%pres": p0, "patch_icpp(1)%alpha_rho(1)": rho0})
    t_end = args.tr * 2 * half * R / u_ref
    cfl_run = args.steps == 0

if cfl_run:
    case.update({"cfl_adap_dt": "T", "cfl_target": 0.5, "n_start": 0, "t_stop": t_end, "t_save": t_end if vortex else t_end / 2})
else:
    dt = 0.5 * dx / (c0 + u_ref)
    case.update({"dt": dt, "t_step_start": 0, "t_step_stop": args.steps, "t_step_save": args.steps // 2})

if __name__ == "__main__":
    print(json.dumps(case))
