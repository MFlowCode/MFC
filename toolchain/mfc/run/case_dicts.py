"""
MFC Case Parameter Type Definitions.

This module provides exports from the central parameter registry (mfc.params).
All parameter definitions are sourced from the registry.

Exports:
    ALL: Family-aware mapping of all parameters {name: ParamType}
    IGNORE: Parameters to skip during certain operations
    CASE_OPTIMIZATION: Parameters that can be hard-coded for GPU builds
    SCHEMA: JSON schema for fastjsonschema validation
    get_validator(): Returns compiled JSON schema validator
    is_valid_for_target(): Check if a parameter is valid for a target
"""

import re

from ..state import ARG


class _ParamTypeMapping:
    """
    Mapping-like object that wraps REGISTRY for {name: ParamType} access.

    Supports O(1) containment checks for indexed family params (e.g.,
    ``'patch_ib(500000)%geometry' in ALL`` is True) without enumerating
    all possible indices.

    For iteration, yields scalar params plus one example per family attr.
    """

    def __init__(self):
        from ..params import REGISTRY

        self._registry = REGISTRY

    def __contains__(self, key):
        return self._registry.is_known_param(key)

    def __getitem__(self, key):
        param_def = self._registry.get_param_def(key)
        if param_def is None:
            raise KeyError(key)
        return param_def.param_type

    def get(self, key, default=None):
        param_def = self._registry.get_param_def(key)
        return param_def.param_type if param_def is not None else default

    def keys(self):
        return self._registry.all_params.keys()

    def values(self):
        return (p.param_type for p in self._registry.all_params.values())

    def items(self):
        return ((n, p.param_type) for n, p in self._registry.all_params.items())

    def __iter__(self):
        return iter(self._registry.all_params)

    def __len__(self):
        return len(self._registry.all_params)


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

# Family-aware mapping of all parameters — supports O(1) lookup for indexed families
ALL = _ParamTypeMapping()

# Parameters that can be hard-coded for GPU case optimization
CASE_OPTIMIZATION = _load_case_optimization_params()

# JSON schema for validation
SCHEMA = _build_schema()

# Regex to extract the base name from indexed params
_BASE_NAME_RE = re.compile(r"^([a-zA-Z_][a-zA-Z0-9_]*)")


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
    match = _BASE_NAME_RE.match(param_name)
    if match:
        base_name = match.group(1)
        return base_name in target_params

    return param_name in target_params


class _TargetKeySet:
    """
    Set-like object for checking if a param is valid for a specific target.

    Supports ``key in target_key_set`` via base-name matching against the
    Fortran namelist, plus optionally filtering out case-optimization params.
    Does not enumerate all possible indexed family members.
    """

    def __init__(self, target_name: str, filter_case_opt: bool = False):
        self._target_name = target_name
        self._filter_case_opt = filter_case_opt
        self._case_opt = set(CASE_OPTIMIZATION) if filter_case_opt else set()

    def __contains__(self, key):
        if self._filter_case_opt and key in self._case_opt:
            return False
        return _is_param_valid_for_target(key, self._target_name)


def get_input_dict_keys(target_name: str):
    """
    Get a set-like object for checking parameter validity for a target.

    Returns an object that supports ``key in result`` for O(1) checks.
    For indexed families, this does NOT enumerate all possible indices —
    it checks the base name against the Fortran namelist.

    Args:
        target_name: One of 'pre_process', 'simulation', 'post_process'

    Returns:
        Set-like object supporting ``in`` operator
    """
    filter_case_opt = ARG("case_optimization", dflt=False) and target_name == "simulation"
    return _TargetKeySet(target_name, filter_case_opt)


def get_validator():
    """Get the cached JSON schema validator."""
    return _get_validator_func()
