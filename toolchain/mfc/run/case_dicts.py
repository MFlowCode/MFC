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

import re
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


def _get_target_params():
    """Get valid params for each target by parsing Fortran namelists."""
    from ..params.namelist_parser import get_target_params
    return get_target_params()


# Parameters to ignore during certain operations
IGNORE = ["cantera_file", "chemistry"]

# Combined dict of all parameters
ALL = _load_all_params()

# Parameters that can be hard-coded for GPU case optimization
CASE_OPTIMIZATION = _load_case_optimization_params()

# JSON schema for validation
SCHEMA = _build_schema()


def _is_param_valid_for_target(param_name: str, target_name: str) -> bool:
    """
    Check if a parameter is valid for a given target.

    Uses the Fortran namelist definitions as the source of truth.
    Handles indexed params like "patch_icpp(1)%geometry" by checking base name.

    Args:
        param_name: The parameter name (may include indices)
        target_name: One of 'pre_process', 'simulation', 'post_process'

    Returns:
        True if the parameter is valid for the target
    """
    target_params = _get_target_params().get(target_name, set())

    # Extract base parameter name (before any index or attribute)
    # e.g., "patch_icpp(1)%geometry" -> "patch_icpp"
    # e.g., "fluid_pp(2)%gamma" -> "fluid_pp"
    # e.g., "acoustic(1)%loc(1)" -> "acoustic"
    match = re.match(r'^([a-zA-Z_][a-zA-Z0-9_]*)', param_name)
    if match:
        base_name = match.group(1)
        return base_name in target_params

    return param_name in target_params


def get_input_dict_keys(target_name: str) -> list:
    """
    Get parameter keys for a given target.

    Uses the Fortran namelist definitions as the source of truth.
    Only returns params whose base name is in the target's namelist.

    Args:
        target_name: One of 'pre_process', 'simulation', 'post_process'

    Returns:
        List of parameter names valid for that target
    """
    keys = [k for k in ALL.keys() if _is_param_valid_for_target(k, target_name)]

    # Case optimization filtering for simulation
    if ARG("case_optimization", dflt=False) and target_name == "simulation":
        keys = [k for k in keys if k not in CASE_OPTIMIZATION]

    return keys


def get_validator():
    """Get the cached JSON schema validator."""
    return _get_validator_func()
