"""
MFC Case Parameter Type Definitions.

This module provides backward-compatible exports from the central parameter
registry (mfc.params). All parameter definitions are now sourced from the
registry, eliminating the previous dual type system.

Exports:
    COMMON, PRE_PROCESS, SIMULATION, POST_PROCESS: Stage-specific param dicts
    ALL: Combined dict of all parameters
    IGNORE: Parameters to skip during certain operations
    CASE_OPTIMIZATION: Parameters that can be hard-coded for GPU builds
    SCHEMA: JSON schema for fastjsonschema validation
    get_validator(): Returns compiled JSON schema validator
    get_input_dict_keys(): Get parameter keys for a target
"""
# pylint: disable=import-outside-toplevel

from ..state import ARG


def _load_stage_dicts():
    """
    Load parameter definitions from the central registry.

    Returns dicts mapping parameter names to their ParamType.
    Uses caching to avoid repeated iteration over ~3300 parameters.
    """
    from functools import lru_cache
    from ..params import REGISTRY
    from ..params.schema import Stage

    @lru_cache(maxsize=8)
    def params_for_stage(stage, include_common=True):
        """Get params for a stage as {name: ParamType} dict (cached)."""
        result = {}
        for name, param in REGISTRY.all_params.items():
            if stage in param.stages:
                result[name] = param.param_type
            elif include_common and Stage.COMMON in param.stages:
                result[name] = param.param_type
        return result

    common = params_for_stage(Stage.COMMON, include_common=False)
    pre = params_for_stage(Stage.PRE_PROCESS, include_common=True)
    sim = params_for_stage(Stage.SIMULATION, include_common=True)
    post = params_for_stage(Stage.POST_PROCESS, include_common=True)

    return common, pre, sim, post


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

# Load parameter definitions from registry
COMMON, PRE_PROCESS, SIMULATION, POST_PROCESS = _load_stage_dicts()

# Parameters to ignore during certain operations
IGNORE = ["cantera_file", "chemistry"]

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
    result = {
        "pre_process": PRE_PROCESS,
        "simulation": SIMULATION,
        "post_process": POST_PROCESS
    }.get(target_name, {}).keys()

    if not ARG("case_optimization") or target_name != "simulation":
        return list(result)

    return [x for x in result if x not in CASE_OPTIMIZATION]


def get_validator():
    """Get the cached JSON schema validator."""
    return _get_validator_func()
