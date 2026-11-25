#!/usr/bin/env python3
"""Axisymmetric focused transducer (support 6) in water; check_focus.py compares the focal pressure with the on-axis
solution of a spherical cap. --probe records the pressure at the focus every step."""

import argparse
import json

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mfc", type=json.loads, default={})
parser.add_argument("--refine", type=int, default=1, help="Refinement factor of the 100 x 50 grid (dt scales with it)")
parser.add_argument("--probe", action="store_true", help="Write the focal pressure (D/probe1_prim.dat) every step")
args = parser.parse_args()
k = args.refine
# Probe in the cell just past the focus (source x + focal length) on the axis; off the cell faces
probe = {"probe_wrt": "T", "fd_order": 2, "num_probes": 1, "probe(1)%x": 0.00006 + 0.00054 + 2.5e-6 / k, "probe(1)%y": 2.5e-6 / k} if args.probe else {}

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 0.001,
            "y_domain%beg": 0.0,
            "y_domain%end": 0.0005,
            "cyl_coord": "T",
            "m": 100 * k - 1,
            "n": 50 * k - 1,
            "p": 0,
            "dt": 2e-9 / k,
            "t_step_start": 0,
            "t_step_stop": 400 * k,
            "t_step_save": 20 * k,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": "5eq",
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": "rk3",
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "teno": "T",
            "teno_CT": 1e-8,
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": "hllc",
            "wave_speeds": "direct",
            "avg_state": "arithmetic",
            "bc_x%beg": -6,
            "bc_x%end": -6,
            "bc_y%beg": -2,
            "bc_y%end": -6,
            # Formatted Database Files Structure Parameters
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1 Liquid
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0005,
            "patch_icpp(1)%y_centroid": 0.00025,
            "patch_icpp(1)%length_x": 0.001,
            "patch_icpp(1)%length_y": 0.0005,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1e05,
            "patch_icpp(1)%alpha_rho(1)": 1100,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Acoustic source
            "acoustic_source": "T",
            "num_source": 1,
            "acoustic(1)%support": 6,
            "acoustic(1)%loc(1)": 0.00006,
            "acoustic(1)%loc(2)": 0.0,
            "acoustic(1)%pulse": 2,
            "acoustic(1)%npulse": 1,
            "acoustic(1)%mag": 1.0,
            "acoustic(1)%gauss_sigma_time": 2e-8,
            "acoustic(1)%foc_length": 0.00054,
            "acoustic(1)%aperture": 0.0008,
            "acoustic(1)%delay": 1e-7,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%pi_inf": 4.4e00 * 5.57e08 / (4.4e00 - 1.0e00),
            **probe,
        }
    )
)
