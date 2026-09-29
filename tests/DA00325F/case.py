#!/usr/bin/env python3
#
# /fastscratch/bwilfong3/software/MFC-Wilfong/tests/DA00325F/case.py:
# 3D -> Example -> TaylorGreenVortex_projection

import json
import argparse

parser = argparse.ArgumentParser(
    prog="/fastscratch/bwilfong3/software/MFC-Wilfong/tests/DA00325F/case.py",
    description="/fastscratch/bwilfong3/software/MFC-Wilfong/tests/DA00325F/case.py: 3D -> Example -> TaylorGreenVortex_projection",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default='{}', metavar="DICT",
                    help="MFC's toolchain's internal state.")

ARGS = vars(parser.parse_args())

case = {
    "run_time_info": "F",
    "m": 25,
    "n": 25,
    "p": 25,
    "dt": 0.0012643239977267494,
    "t_step_start": 0,
    "t_step_stop": 21,
    "t_step_save": 21,
    "num_patches": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "recon_type": 1,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "format": 1,
    "precision": 2,
    "patch_icpp(1)%pres": 101325.0,
    "patch_icpp(1)%alpha_rho(1)": 1.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(2)%pres": 0.5,
    "patch_icpp(2)%alpha_rho(1)": 0.5,
    "patch_icpp(2)%alpha(1)": 1.0,
    "patch_icpp(3)%pres": 0.1,
    "patch_icpp(3)%alpha_rho(1)": 0.125,
    "patch_icpp(3)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 2.5000000000000004,
    "fluid_pp(1)%eos": 1,
    "fluid_pp(1)%cv": 0.0,
    "fluid_pp(1)%qv": 0.0,
    "fluid_pp(1)%qvp": 0.0,
    "bubbles_euler": "F",
    "bubble_model": 3,
    "polytropic": "T",
    "polydisperse": "F",
    "thermal": 3,
    "patch_icpp(1)%r0": 1,
    "patch_icpp(1)%v0": 0,
    "patch_icpp(2)%r0": 1,
    "patch_icpp(2)%v0": 0,
    "patch_icpp(3)%r0": 1,
    "patch_icpp(3)%v0": 0,
    "qbmm": "F",
    "dist_type": 2,
    "poly_sigma": 0.3,
    "sigR": 0.1,
    "sigV": 0.1,
    "rhoRV": 0.0,
    "acoustic_source": "F",
    "num_source": 1,
    "acoustic(1)%loc(1)": 0.5,
    "acoustic(1)%mag": 0.2,
    "acoustic(1)%length": 0.25,
    "acoustic(1)%dir": 1.0,
    "acoustic(1)%npulse": 1,
    "acoustic(1)%pulse": 1,
    "rdma_mpi": "F",
    "bubbles_lagrange": "F",
    "lag_params%nBubs_glb": 1,
    "lag_params%solver_approach": 0,
    "lag_params%cluster_type": 2,
    "lag_params%pressure_corrector": "F",
    "lag_params%smooth_type": 1,
    "lag_params%epsilonb": 1.0,
    "lag_params%heatTransfer_model": "F",
    "lag_params%massTransfer_model": "F",
    "lag_params%valmaxvoid": 0.9,
    "x_domain%beg": -3.141592653589793,
    "x_domain%end": 3.141592653589793,
    "y_domain%beg": -3.141592653589793,
    "y_domain%end": 3.141592653589793,
    "z_domain%beg": -3.141592653589793,
    "z_domain%end": 3.141592653589793,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "bc_z%beg": -1,
    "bc_z%end": -1,
    "viscous": "T",
    "proj_method": "T",
    "low_Mach": 0,
    "parallel_io": "F",
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%hcid": 380,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": 6.283185307179586,
    "patch_icpp(1)%length_y": 6.283185307179586,
    "patch_icpp(1)%length_z": 6.283185307179586,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "fluid_pp(1)%pi_inf": 35109112.500000015,
    "fluid_pp(1)%Re(1)": 42.48128632361878,
    "file_per_process": "F"
}
mods = {}

if "post_process" in ARGS["mfc"]["targets"]:
    mods = {"parallel_io": "T", "cons_vars_wrt": "T", "prim_vars_wrt": "T", "alpha_rho_wrt(1)": "T", "rho_wrt": "T", "mom_wrt(1)": "T", "vel_wrt(1)": "T", "E_wrt": "T", "pres_wrt": "T", "alpha_wrt(1)": "T", "gamma_wrt": "T", "heat_ratio_wrt": "T", "pi_inf_wrt": "T", "pres_inf_wrt": "T", "c_wrt": "T"}
    if case['p'] != 0:
        mods.update({"fd_order": 1, "omega_wrt(1)": "T", "omega_wrt(2)": "T", "omega_wrt(3)": "T"})
else:
    mods = {"parallel_io": "F", "prim_vars_wrt": "F"}

print(json.dumps({**case, **mods}))
