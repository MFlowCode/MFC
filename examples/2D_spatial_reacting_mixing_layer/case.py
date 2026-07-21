#!/usr/bin/env python3
"""Spatially-evolving reacting (H2/air) mixing layer: a 1-D flamelet cross-stream
profile provides both the initial condition and (via bc_x%beg=-17's auto-derivation
from the IC at that face) the inflow condition, and the Wei & Freund (JFM 2005) body
force (bf_spatial_support) seeds the downstream instability.

AXIS CONVENTION -- the reverse of the temporal case (2D_reacting_mixing_layer)'s x/y
roles, and deliberately so: `spatial_bf`'s forcing in src/simulation/m_body_forces.fpp
hardcodes MFC's x as the advecting/streamwise direction (theta_x includes a conv_vel*t
term; theta_y does not), so here:
    MFC x  ->  streamwise (Dirichlet inflow at x-beg, outflow at x-end)
    MFC y  ->  cross-stream (free/extrapolation both ends)
This uses hcid=274 (full 2D field, no extrusion), not hcid=273 -- hcid=273's extrusion
mechanically requires cross-stream=x, which conflicts with spatial_bf's fixed axis
convention. See flamelet_ic.py's write_hcid274_ic docstring.
"""

import argparse
import contextlib
import json
import os
import sys

import cantera as ct
import flamelet_ic
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
ctfile = "h2o2.yaml"

parser = argparse.ArgumentParser(prog="2D_spatial_reacting_mixing_layer", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--scale", type=float, default=1.0, help="Scales grid resolution; use <1 for cheap runs.")
parser.add_argument(
    "--hot",
    action="store_true",
    help="Run the full flamelet Newton/BDF solve for a physically-converged reacting profile (slow; skipped by default). See 2D_reacting_mixing_layer/case.py for why the default is cold.",
)
args = parser.parse_args()

# Physical parameters: same representative H2/air mixing layer as the temporal case.
pressure = 101_325.0
temperature_ox = 500.0
temperature_fu = 300.0
fuel = "H2"
mole_fraction_ox = 0.21
mole_fraction_fu = 1.0
vort_thickness = 1.0e-3
# Spatial (lab-frame, co-flowing) Mach numbers -- unlike the temporal case's single
# symmetric mach_c (u_ox = -u_fu, zero mean by construction), a spatial mixing layer
# needs a genuinely nonzero mean/convective velocity, so both streams flow downstream
# at their own absolute speed. H2's sound speed (~1300 m/s at 300K) is ~3x air's
# (~450 m/s at 500K), so equal-ish Mach numbers would give a near-zero or even
# wrong-signed delta_u; mach_fu is kept low to give a substantial, positive
# delta_u = u_ox - u_fu (oxidizer coflow faster than the fuel jet).
mach_ox = 0.9
mach_fu = 0.1
num_iter = 5

# Grid: x = streamwise (inflow at 0), y = cross-stream.
stream_min, stream_max = 0.0, 15.0
cross_min, cross_max = -10.0, 10.0
points_per_stream = 8.0 * args.scale
points_per_cross = 8.0 * args.scale

t_step_stop = 200
t_step_save = 100

sol = ct.Solution(ctfile)
stream_coord, cross_coord, grid = flamelet_ic.compute_grid_spatial(vort_thickness, cross_min, cross_max, points_per_cross, stream_min, stream_max, points_per_stream)
fluid = flamelet_ic.reference_fluid_properties(sol, temperature_ox, pressure, mole_fraction_ox)

# Cheap (no flamelet solve), needed for spatial_bf%conv_vel and the forcing frequency
# ladder regardless of whether IC/ is cached. Redirected: case.py's contract requires
# its entire stdout to be exactly one JSON line, and streams() prints diagnostics.
with contextlib.redirect_stdout(sys.stderr):
    stream_ox, stream_fu, _ = flamelet_ic.streams(sol, fuel, pressure, temperature_ox, temperature_fu, mole_fraction_ox, mole_fraction_fu, vort_thickness, mach_ox, mach_fu)
u_ox, u_fu = stream_ox[4], stream_fu[4]
conv_vel = 0.5 * (u_ox + u_fu)

# Forcing frequency ladder (Ho & Huerre-style Strouhal scaling), matching the technique
# used by the tool this example was adapted from. Seeded for reproducibility -- this
# case is a registered regression test, so the forcing pattern must be deterministic.
strouhal = 0.032
freq_nom = 4 * conv_vel * strouhal / vort_thickness
rng = np.random.default_rng(seed=0)
bf_freq = 2 * np.pi * 0.25 * freq_nom * (np.arange(8) + rng.uniform(-0.5, 0.5, size=8))
bf_phase = rng.uniform(0, 2 * np.pi, size=8)

ic_dir = os.path.join(current_dir, "IC")
# Key the cache on grid size + mode + physics so a cached IC isn't silently reused across
# a --hot/cold switch or a physical-parameter change that leaves the line count unchanged.
cache_key = {
    "cold": not args.hot,
    "lines": len(stream_coord) * len(cross_coord),
    "vort_thickness": vort_thickness,
    "temperature_ox": temperature_ox,
    "temperature_fu": temperature_fu,
    "mach_ox": mach_ox,
    "mach_fu": mach_fu,
    "mole_fraction_ox": mole_fraction_ox,
    "mole_fraction_fu": mole_fraction_fu,
    "num_iter": num_iter,
}
if not flamelet_ic.ic_cache_valid(ic_dir, "000000", len(stream_coord) * len(cross_coord), cache_key):
    import jax.numpy as jnp
    from pyrometheus.codegen.python import PythonCodeGenerator
    from pyrometheus.flamelets.make_pyro import make_pyro_object

    pyro_cls = PythonCodeGenerator.get_thermochem_class(sol)
    pyro_gas = make_pyro_object(pyro_cls, jnp)

    flamelet_ic.generate_ic_files_spatial(
        output_dir=ic_dir,
        sol=sol,
        pyro_gas=pyro_gas,
        stream_coord=stream_coord,
        cross_coord=cross_coord,
        pressure=pressure,
        temperature_ox=temperature_ox,
        temperature_fu=temperature_fu,
        fuel=fuel,
        mole_fraction_ox=mole_fraction_ox,
        mole_fraction_fu=mole_fraction_fu,
        vort_thickness=vort_thickness,
        mach_ox=mach_ox,
        mach_fu=mach_fu,
        num_iter=num_iter,
        cold=not args.hot,
    )
    flamelet_ic.write_cache_key(ic_dir, cache_key)

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
    "bc_x%beg": -17,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
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
    "bf_spatial_support": "T",
    "spatial_bf%amp": 1.0e2,
    "spatial_bf%x_centroid": 7.5 * vort_thickness,
    "spatial_bf%y_centroid": 0.0,
    "spatial_bf%conv_vel": conv_vel,
    "spatial_bf%sigma": 1.0 / (1.5 * vort_thickness) ** 2,
    **{f"spatial_bf%freq({i + 1})": float(bf_freq[i]) for i in range(8)},
    **{f"spatial_bf%phase({i + 1})": float(bf_phase[i]) for i in range(8)},
    "fluid_pp(1)%gamma": 1.0 / (fluid["gamma"] - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / fluid["viscosity"],
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%hcid": 274,
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
