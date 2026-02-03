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
_SIMULATION_ONLY = {
    "run_time_info", "dt", "t_step_stop", "t_step_save", "t_step_print",
    "time_stepper", "weno_eps", "teno_CT", "wenoz_q", "mapped_weno", "wenoz",
    "teno", "mp_weno", "weno_avg", "weno_Re_flux", "riemann_solver", "wave_speeds",
    "avg_state", "prim_vars_wrt", "alt_soundspeed", "null_weights", "mixture_err",
    "fd_order", "num_probes", "probe_wrt", "bubble_model", "acoustic_source",
    "num_source", "integral_wrt", "num_integrals", "rdma_mpi", "adap_dt",
    "adap_dt_tol", "adap_dt_max_iters", "t_stop", "t_save", "cfl_target",
    "low_Mach", "viscous", "powell", "tau_star", "cont_damage_s", "alpha_bar",
    "num_igr_iters", "num_igr_warm_start_iters", "alf_factor", "igr_iter_solver",
    "igr_pres_lim", "muscl_lim", "int_comp", "ic_eps", "ic_beta",
    "nv_uvm_out_of_core", "nv_uvm_igr_temps_on_gpu", "nv_uvm_pref_gpu",
}

_POST_PROCESS_ONLY = {
    "format", "precision_post", "coarsen_silo", "fourier_modes", "alpha_rho_wrt",
    "rho_wrt", "mom_wrt", "vel_wrt", "flux_lim", "flux_wrt", "E_wrt", "pres_wrt",
    "alpha_wrt", "gamma_wrt", "heat_ratio_wrt", "pi_inf_wrt", "pres_inf_wrt",
    "cons_vars_wrt", "prim_vars_wrt_post", "c_wrt", "omega_wrt", "schlieren_wrt",
    "schlieren_alpha", "fd_order_post", "mag_wrt", "kappa_wrt", "tensor_wrt",
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

    # Filter out params that don't belong to this target
    if target_name == "pre_process":
        keys = [k for k in keys if k not in _SIMULATION_ONLY and k not in _POST_PROCESS_ONLY]
    elif target_name == "post_process":
        keys = [k for k in keys if k not in _SIMULATION_ONLY]
    # simulation gets all params

    # Case optimization filtering for simulation
    if ARG("case_optimization", dflt=False) and target_name == "simulation":
        keys = [k for k in keys if k not in CASE_OPTIMIZATION]

    return keys


def get_validator():
    """Get the cached JSON schema validator."""
    return _get_validator_func()
