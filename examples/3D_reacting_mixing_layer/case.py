#!/usr/bin/env python3
"""3D temporal reacting (H2/air) mixing layer, initialized from a 1-D flamelet solve
extruded via hcid=371. See flamelet_ic.py's module docstring for the axis convention
(x = streamwise, y = cross-stream/profile, z = spanwise) and for how the 3D velocity
perturbation that seeds the cascade is baked in: in-plane (x,y) via flamelet_ic.py's
perturb_xy, spanwise (z) via a closed-form modulation added directly in Fortran.
"""

import argparse
import json
import os

import cantera as ct
import flamelet_ic

current_dir = os.path.dirname(os.path.abspath(__file__))
ctfile = os.path.join(current_dir, "sandiego.yaml")

parser = argparse.ArgumentParser(prog="3D_reacting_mixing_layer", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--scale", type=float, default=1.0, help="Scales cross-stream grid resolution; use <1 for cheap runs.")
# See examples/2D_reacting_mixing_layer/case.py for why the default is the cold (mollified,
# non-reacting) profile and --hot runs the full flamelet Newton/BDF solve.
parser.add_argument("--hot", action="store_true", help="Run the full flamelet Newton/BDF solve for a physically-converged reacting profile (slow; skipped by default).")
args = parser.parse_args()

# Physical parameters: same H2/air mixing layer as the 2D temporal case, but supersonic
# (mach_c=1.5 instead of 0.3) -- this case exists specifically to probe compressibility
# effects the 2D/subsonic case can't. Uses the San Diego mechanism (sandiego.yaml, shipped
# alongside this case) rather than the 2D case's h2o2.yaml -- all tuning below (dilution,
# resolution, flamelet chi) is calibrated against it; swapping mechanisms invalidates it.
pressure = 101_325.0
temperature_ox = 500.0
temperature_fu = 300.0
fuel = "H2"
mole_fraction_ox = 0.21
# Fuel stream diluted with N2 (was 1.0, pure H2): pure H2 has such a low molecular weight
# that its sound speed (~1318 m/s at 300K) is over 4x c_ox, so hitting mach_c=1.5 via the
# mean-sound-speed convective-Mach definition below required an unphysically extreme
# absolute velocity split (delta_u ~ 2650 m/s, Mach ~5.9 relative to the oxidizer stream) --
# the two 3D production runs both crashed (VCFL=Inf) under this, and doubling y-resolution
# barely changed the crash timing, pointing at the base flow itself rather than under-
# resolution. Diluting to X_H2=0.5 (c_fu ~ 483 m/s, delta_u ~ 1394 m/s, Ma ~ 3.1 in the
# oxidizer frame) also moves Z_st from near 0 (pure-fuel-side) to ~0.30, away from the
# domain edge.
mole_fraction_fu = 0.5
vort_thickness = 1.0e-3
mach_c = 1.5
num_iter = 5

# Grid: x = streamwise (periodic), y = cross-stream (flamelet profile axis), z = spanwise
# (periodic). Quarter-scale of Wang et al. (C&F 2024)'s temporal mixing-layer DNS gold
# standard: their domain is [60,80,40] delta_omega at a uniform 14 points/delta_omega in
# all three directions (840x1120x560 -- we can't afford that budget-wise, so we keep their
# resolution DENSITY and shrink the box to 1/4 in each direction: [15,20,10] delta_omega,
# 210x280x140 at 14 pts/delta_omega uniform. That run crashed (VCFL=Inf, ~t_step=33000)
# after developing Mach-1.5 shocklets in the braids (local velocities ~1500 m/s, 4-5x the
# base stream speed) -- points_per_cross (below) doubles y's density on top of this
# baseline to better resolve those gradients, independent of x/z.
points_per_delta = 14.0 * args.scale
cross_min, cross_max = -10.0, 10.0  # y: 20 delta_omega
# y carries the flame/shear-layer gradients (and now the shocklets from the Mach-1.5
# braids that crashed the previous run) -- doubled relative to x/z instead of touching the
# numerics (Riemann solver, WENO flux type, reaction sub-stepping) or backing off the
# perturbation amplitude again. dt is left at 1e-9 s unchanged (per prior 2D-case experience).
points_per_cross = 2.0 * points_per_delta
stream_min, stream_max = -7.5, 7.5  # x: 15 delta_omega
num_x = round((stream_max - stream_min) * points_per_delta)
span_min, span_max = -5.0, 5.0  # z: 10 delta_omega
num_z = round((span_max - span_min) * points_per_delta)

t_step_stop = 20000000
# Was 10000 -- at the quarter-scale grid's measured 1-GPU rate (3.82 s/step), a job would
# need >10.6 hours just to reach the FIRST checkpoint. The run that just timed out spent 4
# GPU-hours reaching only 3766/10000 steps -- zero recoverable progress. Lowered so real
# progress gets banked well within a single job's walltime even at conservative (non-ideal)
# multi-GPU scaling.
t_step_save = 5000

sol = ct.Solution(ctfile)
cross_coord, x_coord, grid = flamelet_ic.compute_grid_3d(vort_thickness, cross_min, cross_max, points_per_cross, stream_min, stream_max, num_x, span_min, span_max, num_z)
fluid = flamelet_ic.reference_fluid_properties(sol, temperature_ox, pressure, mole_fraction_ox)

ic_dir = os.path.join(current_dir, "IC")
# Key the cache on grid size + mode + physics so a cached IC isn't silently reused across
# a --hot/cold switch or a physical-parameter change that leaves the line count unchanged.
cache_key = {
    "cold": not args.hot,
    "lines": len(x_coord) * len(cross_coord),
    "vort_thickness": vort_thickness,
    "temperature_ox": temperature_ox,
    "temperature_fu": temperature_fu,
    "mach_c": mach_c,
    "mole_fraction_ox": mole_fraction_ox,
    "mole_fraction_fu": mole_fraction_fu,
    "num_iter": num_iter,
}
if not flamelet_ic.ic_cache_valid(ic_dir, "000000", len(x_coord) * len(cross_coord), cache_key):
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
        x_coord=x_coord,
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
    flamelet_ic.write_cache_key(ic_dir, cache_key)

case = {
    "run_time_info": "T",
    "x_domain%beg": grid["x_domain_beg"],
    "x_domain%end": grid["x_domain_end"],
    "y_domain%beg": grid["y_domain_beg"],
    "y_domain%end": grid["y_domain_end"],
    "z_domain%beg": grid["z_domain_beg"],
    "z_domain%end": grid["z_domain_end"],
    "m": grid["m"],
    "n": grid["n"],
    "p": grid["p"],
    "cyl_coord": "F",
    # Scaled down from the previous 1e-9 s in proportion to the finer grid spacing
    # (8 -> 14 pts/delta_omega, dx shrinks by 8/14 ~ 0.571) to hold the same CFL margin.
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
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "bc_z%beg": -1,
    "bc_z%end": -1,
    "num_patches": 1,
    "num_fluids": 1,
    "viscous": "T",
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "T",
    # Unity-Lewis, to match the flamelet solve's own diffusivity assumption (flamelet_ic.py's
    # diffusivity() uses D_k = k/(rho*cp) for every species) -- default is 1 (mixture-averaged),
    # which is what tmp-3D-1 ran with.
    "chem_params%transport_model": 2,
    "files_dir": ic_dir,
    "file_extension": "000000",
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "fluid_pp(1)%gamma": 1.0 / (fluid["gamma"] - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / fluid["viscosity"],
    "patch_icpp(1)%geometry": 9,
    # hcid=371: hcid=370's read (flamelet base state + in-plane (x,y) perturbation from
    # flamelet_ic.py's perturb_xy, ported from ../tmp-mixlyr/perturb.py) plus a closed-form
    # spanwise (z) modulation added directly in Fortran (src/common/include/
    # 3dHardcodedIC.fpp), so the IC has genuine 3D structure from step 0. Not
    # MFC's own mixlayer_perturb: that routine hardcodes an absolute wavenumber range
    # (effectively >6m wavelengths) that's meaningless for this millimeter-scale domain --
    # confirmed by measuring essentially zero real x/z structure despite large point-wise
    # perturbation variance, both with the Fortran default mixlayer_perturb_k0 and after
    # trying to rescale it (that parameter only reweights an already length-scale-fixed
    # wavenumber range, it doesn't rescale it).
    "patch_icpp(1)%hcid": 371,
    "patch_icpp(1)%x_centroid": 0.5 * (grid["x_domain_beg"] + grid["x_domain_end"]),
    "patch_icpp(1)%y_centroid": 0.5 * (grid["y_domain_beg"] + grid["y_domain_end"]),
    "patch_icpp(1)%z_centroid": 0.5 * (grid["z_domain_beg"] + grid["z_domain_end"]),
    "patch_icpp(1)%length_x": grid["x_domain_end"] - grid["x_domain_beg"],
    "patch_icpp(1)%length_y": grid["y_domain_end"] - grid["y_domain_beg"],
    "patch_icpp(1)%length_z": grid["z_domain_end"] - grid["z_domain_beg"],
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": pressure,
    "patch_icpp(1)%alpha_rho(1)": 1,
    "patch_icpp(1)%alpha(1)": 1,
    "cantera_file": ctfile,
}

print(json.dumps(case))
