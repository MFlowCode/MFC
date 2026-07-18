#!/usr/bin/env python3
"""Temporal reacting (H2/air) mixing layer, initialized from a 1-D flamelet solve
extruded via hcid=273. See flamelet_ic.py's module docstring for the required MFC
axis convention (x = cross-stream/profile, y = streamwise/periodic/extruded).
"""

import argparse
import json
import os

import cantera as ct
import flamelet_ic

current_dir = os.path.dirname(os.path.abspath(__file__))
ctfile = "h2o2.yaml"

parser = argparse.ArgumentParser(prog="2D_reacting_mixing_layer", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--scale", type=float, default=1.0, help="Scales cross-stream grid resolution; use <1 for cheap runs.")
# Default is the cold (non-reacting, mollified) profile: the flamelet Newton/BDF solve
# runs on a fixed 101-point mixture-fraction grid independent of --scale, so it is not
# cheap to skip via a smaller physical grid -- and this file must load fast every time,
# since it's invoked by `./mfc.sh validate`/precheck (batched over all examples) and
# twice per `./mfc.sh run` (see toolchain/mfc/run/run.py + build.py). Pass --hot for the
# real, physically-converged reacting profile.
parser.add_argument("--hot", action="store_true", help="Run the full flamelet Newton/BDF solve for a physically-converged reacting profile (slow; skipped by default).")
args = parser.parse_args()

# Physical parameters: representative temporal H2/air mixing layer.
pressure = 101_325.0
temperature_ox = 500.0
temperature_fu = 300.0
fuel = "H2"
mole_fraction_ox = 0.21
mole_fraction_fu = 1.0
vort_thickness = 1.0e-3
mach_c = 0.3
num_iter = 5

# Grid: x = cross-stream (flamelet profile axis), y = streamwise (periodic, extruded).
cross_min, cross_max = -10.0, 10.0
points_per_cross = 8.0 * args.scale
stream_min, stream_max = -5.0, 5.0
num_y = 32  # weno_order=5 needs n+1 >= 25 (num_stcls_min*weno_order)

t_step_stop = 200
t_step_save = 100

sol = ct.Solution(ctfile)
cross_coord, grid = flamelet_ic.compute_grid(vort_thickness, cross_min, cross_max, points_per_cross, stream_min, stream_max, num_y)
fluid = flamelet_ic.reference_fluid_properties(sol, temperature_ox, pressure, mole_fraction_ox)

ic_dir = os.path.join(current_dir, "IC")
if not flamelet_ic.ic_cache_valid(ic_dir, "000000", len(cross_coord)):
    import jax.numpy as jnp
    from pyrometheus.codegen.python import PythonCodeGenerator
    from pyrometheus.flamelets.make_pyro import make_pyro_object

    pyro_cls = PythonCodeGenerator.get_thermochem_class(sol)
    pyro_gas = make_pyro_object(pyro_cls, jnp)

    flamelet_ic.generate_ic_files(
        output_dir=ic_dir,
        sol=sol,
        pyro_gas=pyro_gas,
        cross_coord=cross_coord,
        pressure=pressure,
        temperature_ox=temperature_ox,
        temperature_fu=temperature_fu,
        fuel=fuel,
        mole_fraction_ox=mole_fraction_ox,
        mole_fraction_fu=mole_fraction_fu,
        vort_thickness=vort_thickness,
        mach_c=mach_c,
        num_iter=num_iter,
        cold=not args.hot,
    )

case = {
    "run_time_info": "T",
    "x_domain%beg": grid["x_domain_beg"],
    "x_domain%end": grid["x_domain_end"],
    "y_domain%beg": grid["y_domain_beg"],
    "y_domain%end": grid["y_domain_end"],
    "m": grid["m"],
    "n": grid["n"],
    "p": 0,
    "cyl_coord": "F",
    "dt": 1.0e-9,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mixture_err": "F",
    "mpp_lim": "F",
    "time_stepper": 3,
    "avg_state": 1,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "weno_Re_flux": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "num_patches": 1,
    "num_fluids": 1,
    "viscous": "T",
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "T",
    "files_dir": ic_dir,
    "file_extension": "000000",
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "fluid_pp(1)%gamma": 1.0 / (fluid["gamma"] - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / fluid["viscosity"],
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%hcid": 273,
    "patch_icpp(1)%x_centroid": 0.5 * (grid["x_domain_beg"] + grid["x_domain_end"]),
    "patch_icpp(1)%y_centroid": 0.5 * (grid["y_domain_beg"] + grid["y_domain_end"]),
    "patch_icpp(1)%length_x": grid["x_domain_end"] - grid["x_domain_beg"],
    "patch_icpp(1)%length_y": grid["y_domain_end"] - grid["y_domain_beg"],
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": pressure,
    "patch_icpp(1)%alpha_rho(1)": 1,
    "patch_icpp(1)%alpha(1)": 1,
    "cantera_file": ctfile,
}

print(json.dumps(case))
