"""
MFC Case Parameter Type Definitions.

This module provides exports from the central parameter registry (mfc.params).
All parameter definitions are sourced from the registry.

Exports:
    ALL: Dict of all parameters {name: ParamType}
    IGNORE: Parameters to skip during certain operations
    CASE_OPTIMIZATION: Parameters that can be hard-coded for GPU builds
    SCHEMA: JSON schema for fastjsonschema validation
    get_validator(): Returns compiled JSON schema validator
    get_input_dict_keys(): Get parameter keys for a target
"""
# pylint: disable=import-outside-toplevel

from ..state import ARG


def _load_all_params():
    """Load all parameters as {name: ParamType} dict."""
    from ..params import REGISTRY
    return {name: param.param_type for name, param in REGISTRY.all_params.items()}


def _load_case_optimization_params():
    """Get params that can be hard-coded for GPU optimization."""
    from ..params import REGISTRY
    return [name for name, param in REGISTRY.all_params.items() if param.case_optimization]


def _build_schema():
    """Build JSON schema from registry."""
    from ..params import REGISTRY
    return REGISTRY.get_json_schema()


def _get_validator_func():
    """Get the cached validator from registry."""
    from ..params import REGISTRY
    return REGISTRY.get_validator()


# Parameters to ignore during certain operations
IGNORE = ["cantera_file", "chemistry"]

# Parameters that are only valid for specific targets (not in their Fortran namelist)
# These get excluded when generating .inp files for other targets
# Extracted from master branch's separate PRE_PROCESS/SIMULATION/POST_PROCESS dicts

# Params ONLY in simulation (not in pre_process OR post_process)
_SIMULATION_ONLY = {
    # Core simulation params
    "run_time_info", "dt", "t_step_print", "time_stepper",
    # WENO/reconstruction (not in post_process)
    "weno_eps", "teno_CT", "wenoz_q", "mapped_weno", "wenoz", "teno",
    "mp_weno", "weno_avg", "weno_Re_flux", "muscl_lim",
    # Riemann solver
    "riemann_solver", "wave_speeds", "low_Mach",
    # Probes/integrals
    "probe_wrt", "num_probes", "integral_wrt", "num_integrals",
    # Physics
    "null_weights", "bubble_model", "acoustic_source", "num_source",
    # Adaptive dt
    "adap_dt", "adap_dt_tol", "adap_dt_max_iters",
    # MPI
    "rdma_mpi",
    # Damage
    "tau_star", "cont_damage_s", "alpha_bar",
    # MHD
    "powell",
    # IGR
    "num_igr_iters", "num_igr_warm_start_iters", "alf_factor",
    "igr_iter_solver", "igr_pres_lim",
    # Interface compression
    "int_comp", "ic_eps", "ic_beta",
    # GPU memory
    "nv_uvm_out_of_core", "nv_uvm_igr_temps_on_gpu", "nv_uvm_pref_gpu",
}

# Params in simulation AND post_process but NOT pre_process
_SIM_AND_POST = {
    "t_step_stop", "t_step_save", "t_stop", "t_save", "cfl_target", "n_start",
    "alt_soundspeed", "mixture_err", "fd_order", "avg_state", "prim_vars_wrt",
}

_POST_PROCESS_ONLY = {
    # Output format
    "format", "coarsen_silo", "fourier_modes",
    # Field output flags
    "alpha_rho_wrt", "rho_wrt", "mom_wrt", "vel_wrt", "E_wrt", "pres_wrt",
    "alpha_wrt", "gamma_wrt", "heat_ratio_wrt", "pi_inf_wrt", "pres_inf_wrt",
    "cons_vars_wrt", "c_wrt", "omega_wrt", "qm_wrt", "liutex_wrt",
    "schlieren_wrt", "schlieren_alpha", "kappa_wrt",
    "flux_lim", "flux_wrt", "cf_wrt",
    # Lagrange bubble output
    "lag_header", "lag_txt_wrt", "lag_id_wrt", "lag_pos_wrt", "lag_pos_prev_wrt",
    "lag_vel_wrt", "lag_rad_wrt", "lag_rvel_wrt", "lag_rmin_wrt", "lag_rmax_wrt",
    "lag_pres_wrt", "lag_db_wrt", "lag_dphidt_wrt", "lag_mv_wrt", "lag_mg_wrt",
    # Other
    "output_partial_domain", "sim_data",
}

# Combined dict of all parameters
ALL = _load_all_params()

# Parameters that can be hard-coded for GPU case optimization
CASE_OPTIMIZATION = _load_case_optimization_params()

# JSON schema for validation
SCHEMA = _build_schema()


def get_input_dict_keys(target_name: str) -> list:
    """
    Get parameter keys for a given target.

    Args:
        target_name: One of 'pre_process', 'simulation', 'post_process'

    Returns:
        List of parameter names valid for that target
    """
    keys = list(ALL.keys())

    # Filter out params that don't belong to this target's Fortran namelist
    if target_name == "pre_process":
        # pre_process excludes simulation-only, sim+post shared, and post-only params
        exclude = _SIMULATION_ONLY | _SIM_AND_POST | _POST_PROCESS_ONLY
        keys = [k for k in keys if k not in exclude]
    elif target_name == "post_process":
        # post_process excludes only truly simulation-only params
        keys = [k for k in keys if k not in _SIMULATION_ONLY]
    # simulation accepts all params

    # Case optimization filtering for simulation
    if ARG("case_optimization", dflt=False) and target_name == "simulation":
        keys = [k for k in keys if k not in CASE_OPTIMIZATION]

    return keys


def get_validator():
    """Get the cached JSON schema validator."""
    return _get_validator_func()
