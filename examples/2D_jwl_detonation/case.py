#!/usr/bin/env python3
# Real-world 2D detonation: a cylindrical TNT charge detonating in free air,
# driven by a Rocflu-style program burn on a moderate uniform grid.
#
# The charge is a low-density products reservoir at ambient pressure (with the
# products JWL EOS a full-density 1630 kg/m^3 charge already carries ~6 GPa of
# cold-curve pressure, so the unreacted charge is represented as an expanded
# reservoir at ambient p). A lighting-time front then sweeps outward from the
# centre at the TNT CJ speed D_cj and releases the detonation energy E0/rho0 over
# a reaction band of width pb_width, driving the blast into the surrounding air.
#
# NOTE on grid stretching: this case was first built with a stretched mesh (fine
# core, coarse far field). Grid stretching proved incompatible with the strong
# program-burn energy source here -- with stretch_x/y enabled the run goes NaN
# within ~16 steps regardless of dt (halving dt does not delay it, so it is not a
# CFL limit) or reaction-band width, while the uniform mesh below is robust. The
# strong reactive source on the non-uniform (stretched) reconstruction is the
# trigger; the same NaN appears for a high-pressure products burst on a stretched
# mesh. The moderate 200x200 uniform grid is therefore used for the shipped case.
import json

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": -1.0,
            "x_domain%end": 1.0,
            "y_domain%beg": -1.0,
            "y_domain%end": 1.0,
            "m": 199,
            "n": 199,
            "p": 0,
            "dt": 2.5e-8,
            "t_step_start": 0,
            "t_step_stop": 8000,
            "t_step_save": 1000,
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Rocflu-style program burn, initiated at the charge centre.
            "prog_burn": "T",
            "pb_D_cj": 6930.0,
            "pb_width": 0.03,
            "pb_x_det": 0.0,
            "pb_y_det": 0.0,
            # Patch 1: ambient air.
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 2.0,
            "patch_icpp(1)%length_y": 2.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 101325.0,
            "patch_icpp(1)%alpha_rho(1)": 1.63e-5,
            "patch_icpp(1)%alpha_rho(2)": 1.22499998775,
            "patch_icpp(1)%alpha(1)": 1.0e-8,
            "patch_icpp(1)%alpha(2)": 0.99999999,
            # Patch 2: cylindrical charge (products reservoir), radius 0.1 m.
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%radius": 0.1,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 101325.0,
            "patch_icpp(2)%alpha_rho(1)": 99.999999,
            "patch_icpp(2)%alpha_rho(2)": 1.225e-8,
            "patch_icpp(2)%alpha(1)": 0.99999999,
            "patch_icpp(2)%alpha(2)": 1.0e-8,
            # Fluid 1: TNT JWL products.
            "fluid_pp(1)%eos": 2,
            "fluid_pp(1)%gamma": 2.5,
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%cv": 613.5,
            "fluid_pp(1)%jwl_A": 3.712e11,
            "fluid_pp(1)%jwl_B": 3.231e9,
            "fluid_pp(1)%jwl_R1": 4.15,
            "fluid_pp(1)%jwl_R2": 0.95,
            "fluid_pp(1)%jwl_omega": 0.30,
            "fluid_pp(1)%jwl_rho0": 1630.0,
            "fluid_pp(1)%jwl_E0": 1.0089e10,
            "fluid_pp(1)%jwl_air_e0": 2.5575e5,
            "fluid_pp(1)%jwl_air_rho0": 1.225,
            # Fluid 2: air.
            "fluid_pp(2)%eos": 1,
            "fluid_pp(2)%gamma": 2.5,
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%cv": 717.5,
        }
    )
)
