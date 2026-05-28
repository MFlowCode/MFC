"""
Namelist target mapping for Fortran codegen.

NAMELIST_VARS maps each Fortran namelist variable (struct root or simple scalar)
to the set of MFC executables whose namelist it appears in.

CASE_OPT_EXCLUDE is the set of simulation namelist variables excluded under
MFC_CASE_OPTIMIZATION (they become compile-time constants instead).

When adding a new parameter:
  1. Add to definitions.py (type, constraints, etc.)
  2. Add the namelist root variable to NAMELIST_VARS with its target set
  3. Run ./mfc.sh generate to regenerate the .fpp files
"""

from typing import Dict, Set

# All three targets
_ALL = {"pre", "sim", "post"}
_PRE_SIM = {"pre", "sim"}
_SIM_POST = {"sim", "post"}

NAMELIST_VARS: Dict[str, Set[str]] = {
    # --- Grid (all targets) ---
    "m": _ALL,
    "n": _ALL,
    "p": _ALL,
    "cyl_coord": _ALL,
    "x_domain": {"pre", "sim"},
    "y_domain": {"pre", "sim"},
    "z_domain": {"pre", "sim"},
    "x_output": {"post"},
    "y_output": {"post"},
    "z_output": {"post"},
    # --- Grid stretching (pre only) ---
    "stretch_x": {"pre"},
    "stretch_y": {"pre"},
    "stretch_z": {"pre"},
    "a_x": {"pre"},
    "a_y": {"pre"},
    "a_z": {"pre"},
    "x_a": {"pre", "sim"},
    "y_a": {"pre", "sim"},
    "z_a": {"pre", "sim"},
    "x_b": {"pre", "sim"},
    "y_b": {"pre", "sim"},
    "z_b": {"pre", "sim"},
    "loops_x": {"pre"},
    "loops_y": {"pre"},
    "loops_z": {"pre"},
    # --- Time ---
    "dt": {"sim"},
    "t_step_start": _ALL,
    "t_step_stop": {"sim", "post"},
    "t_step_save": {"sim", "post"},
    "t_step_print": {"sim"},
    "t_step_old": {"pre", "sim"},
    "time_stepper": {"sim"},
    "t_stop": {"sim", "post"},
    "t_save": {"sim", "post"},
    "cfl_target": {"sim", "post"},
    "cfl_adap_dt": _ALL,
    "cfl_const_dt": _ALL,
    "n_start": _ALL,
    "n_start_old": {"pre"},
    "adap_dt": {"sim"},
    "adap_dt_tol": {"sim"},
    "adap_dt_max_iters": {"sim"},
    # --- Physics model ---
    "model_eqns": _ALL,
    "num_fluids": {"pre", "post"},
    "mpp_lim": _ALL,
    "relax": _ALL,
    "relax_model": _ALL,
    "palpha_eps": _PRE_SIM,
    "ptgalpha_eps": _PRE_SIM,
    # --- WENO / reconstruction ---
    "weno_order": {"pre", "post"},
    "weno_eps": {"sim"},
    "teno_CT": {"sim"},
    "wenoz_q": {"sim"},
    "mp_weno": {"sim"},
    "weno_avg": {"sim"},
    "weno_Re_flux": {"sim"},
    "null_weights": {"sim"},
    "muscl_eps": {"sim"},
    "recon_type": {"pre", "post"},
    "muscl_order": {"pre", "post"},
    "muscl_lim": set(),
    "int_comp": {"sim"},
    "ic_eps": {"sim"},
    "ic_beta": {"sim"},
    # --- Riemann solver ---
    "riemann_solver": {"sim"},
    "wave_speeds": {"sim"},
    "avg_state": {"sim", "post"},
    "low_Mach": {"sim"},
    # --- MHD ---
    "mhd": {"pre", "post"},
    "hyper_cleaning": _ALL,
    "hyper_cleaning_speed": {"sim"},
    "hyper_cleaning_tau": {"sim"},
    "Bx0": _ALL,
    # --- BCs ---
    "bc_x": _ALL,
    "bc_y": _ALL,
    "bc_z": _ALL,
    "num_bc_patches": _ALL,
    "patch_bc": {"pre"},
    # --- ICs (pre only) ---
    "num_patches": {"pre"},
    "patch_icpp": {"pre"},
    # --- Fluid properties ---
    "fluid_pp": _ALL,
    "bub_pp": _ALL,
    "rhoref": _ALL,
    "pref": _ALL,
    # --- Bubbles ---
    "bubbles_euler": _ALL,
    "bubbles_lagrange": _ALL,
    "R0ref": _ALL,
    "nb": {"pre", "post"},
    "polytropic": _ALL,
    "thermal": _ALL,
    "Ca": _ALL,
    "Web": _ALL,
    "Re_inv": _ALL,
    "polydisperse": _ALL,
    "poly_sigma": _ALL,
    "qbmm": _ALL,
    "sigma": _ALL,
    "adv_n": _ALL,
    "bubble_model": {"sim"},
    "sigR": {"pre", "post"},
    "sigV": {"pre"},
    "dist_type": {"pre"},
    "rhoRV": {"pre"},
    "lag_params": {"sim"},
    # --- Lagrangian output (post) ---
    "lag_header": {"post"},
    "lag_txt_wrt": {"post"},
    "lag_db_wrt": {"post"},
    "lag_id_wrt": {"post"},
    "lag_pos_wrt": {"post"},
    "lag_pos_prev_wrt": {"post"},
    "lag_vel_wrt": {"post"},
    "lag_rad_wrt": {"post"},
    "lag_rvel_wrt": {"post"},
    "lag_r0_wrt": {"post"},
    "lag_rmax_wrt": {"post"},
    "lag_rmin_wrt": {"post"},
    "lag_dphidt_wrt": {"post"},
    "lag_pres_wrt": {"post"},
    "lag_mv_wrt": {"post"},
    "lag_mg_wrt": {"post"},
    "lag_betaT_wrt": {"post"},
    "lag_betaC_wrt": {"post"},
    # --- Elasticity ---
    "hypoelasticity": _ALL,
    "hyperelasticity": _ALL,
    # --- Surface tension ---
    "surface_tension": _ALL,
    # --- Relativity ---
    "relativity": _ALL,
    # --- Immersed boundaries ---
    "ib": _ALL,
    "num_ibs": _ALL,
    "patch_ib": {"pre", "sim"},
    "collision_model": {"sim"},
    "coefficient_of_restitution": {"sim"},
    "collision_time": {"sim"},
    "ib_coefficient_of_friction": {"sim"},
    "ib_state_wrt": {"sim", "post"},
    # --- Continuum damage ---
    "cont_damage": _ALL,
    "tau_star": {"sim"},
    "cont_damage_s": {"sim"},
    "alpha_bar": {"sim"},
    # --- IGR ---
    "igr": {"pre", "post"},
    "igr_order": {"pre", "post"},
    "down_sample": _ALL,
    # --- Probes (sim) ---
    "probe_wrt": {"sim"},
    "num_probes": {"sim"},
    "probe": {"sim"},
    "integral_wrt": {"sim"},
    "num_integrals": {"sim"},
    "integral": {"sim"},
    "fd_order": {"sim", "post"},
    # --- Acoustic sources (sim) ---
    "acoustic_source": {"sim"},
    "num_source": {"sim"},
    "acoustic": {"sim"},
    # --- Chemistry ---
    "chem_params": {"sim"},
    # --- Body forces (sim) ---
    "bf_x": {"sim"},
    "bf_y": {"sim"},
    "bf_z": {"sim"},
    "k_x": {"sim"},
    "k_y": {"sim"},
    "k_z": {"sim"},
    "w_x": {"sim"},
    "w_y": {"sim"},
    "w_z": {"sim"},
    "p_x": {"sim"},
    "p_y": {"sim"},
    "p_z": {"sim"},
    "g_x": {"sim"},
    "g_y": {"sim"},
    "g_z": {"sim"},
    # --- Viscous (pre) ---
    "viscous": {"pre"},
    # --- Output ---
    "precision": _ALL,
    "parallel_io": _ALL,
    "file_per_process": _ALL,
    "prim_vars_wrt": {"sim", "post"},
    "cons_vars_wrt": {"post"},
    "run_time_info": {"sim"},
    "fft_wrt": _ALL,
    "pi_fac": {"pre", "sim"},
    # --- Post-process output ---
    "format": {"post"},
    "output_partial_domain": {"post"},
    "rho_wrt": {"post"},
    "E_wrt": {"post"},
    "pres_wrt": {"post"},
    "c_wrt": {"post"},
    "omega_wrt": {"post"},
    "qm_wrt": {"post"},
    "liutex_wrt": {"post"},
    "schlieren_wrt": {"post"},
    "schlieren_alpha": {"post"},
    "gamma_wrt": {"post"},
    "heat_ratio_wrt": {"post"},
    "pi_inf_wrt": {"post"},
    "pres_inf_wrt": {"post"},
    "alpha_rho_wrt": {"post"},
    "mom_wrt": {"post"},
    "vel_wrt": {"post"},
    "flux_wrt": {"post"},
    "alpha_wrt": {"post"},
    "cf_wrt": {"post"},
    "chem_wrt_T": {"post"},
    "chem_wrt_Y": {"post"},
    "alt_soundspeed": {"sim", "post"},
    "mixture_err": {"sim", "post"},
    "flux_lim": {"post"},
    "sim_data": {"post"},
    "alpha_rho_e_wrt": {"post"},
    "G": {"post"},
    # --- Pre-process IC perturbations ---
    "perturb_flow": {"pre"},
    "perturb_flow_fluid": {"pre"},
    "perturb_flow_mag": {"pre"},
    "perturb_sph": {"pre"},
    "perturb_sph_fluid": {"pre"},
    "fluid_rho": {"pre"},
    "mixlayer_vel_profile": {"pre"},
    "mixlayer_vel_coef": {"pre"},
    "mixlayer_perturb": {"pre"},
    "mixlayer_perturb_nk": {"pre"},
    "mixlayer_perturb_k0": {"pre"},
    "pre_stress": {"pre"},
    "elliptic_smoothing": {"pre"},
    "elliptic_smoothing_iters": {"pre"},
    "simplex_perturb": {"pre"},
    "simplex_params": {"pre"},
    # --- Pre-process restart ---
    "old_grid": {"pre"},
    "old_ic": {"pre"},
    # --- Sim-specific physics ---
    "rdma_mpi": {"sim"},
    "alf_factor": {"sim"},
    "num_igr_iters": {"sim"},
    "num_igr_warm_start_iters": {"sim"},
    "igr_iter_solver": {"sim"},
    "igr_pres_lim": {"sim"},
    "nv_uvm_out_of_core": {"sim"},
    "nv_uvm_igr_temps_on_gpu": {"sim"},
    "nv_uvm_pref_gpu": {"sim"},
    # --- Logistics ---
    "case_dir": _ALL,
}

# Variables excluded from the sim namelist when MFC_CASE_OPTIMIZATION is active
# (they become compile-time integer/logical parameters instead).
CASE_OPT_EXCLUDE: Set[str] = {
    "nb",
    "mapped_weno",
    "wenoz",
    "teno",
    "wenoz_q",
    "weno_order",
    "num_fluids",
    "mhd",
    "relativity",
    "igr_order",
    "viscous",
    "igr_iter_solver",
    "igr",
    "igr_pres_lim",
    "recon_type",
    "muscl_order",
    "muscl_lim",
}

# Add CASE_OPT_EXCLUDE vars to NAMELIST_VARS for sim target
# (they appear in the namelist when NOT using case optimization)
for _v in CASE_OPT_EXCLUDE:
    if _v not in NAMELIST_VARS:
        NAMELIST_VARS[_v] = {"sim"}
    else:
        NAMELIST_VARS[_v].add("sim")
