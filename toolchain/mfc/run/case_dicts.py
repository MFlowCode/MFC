"""
MFC Case Parameter Type Definitions.

Exports from the central parameter registry (mfc.params).

Exports:
    ALL: Family-aware mapping of all parameters {name: ParamType}
    IGNORE: Parameters to skip during certain operations
    CASE_OPTIMIZATION: Parameters that can be hard-coded for GPU builds
    SCHEMA: JSON schema for fastjsonschema validation
    get_validator(): Returns compiled JSON schema validator
    get_input_dict_keys(): Get set-like object for target parameter checking
"""

import re
from collections.abc import Mapping

from ..state import ARG


class _ParamTypeMapping(Mapping):
    """Read-only {name: ParamType} view over REGISTRY.all_params.

    Delegates containment and lookup to the registry's family-aware mapping,
    so indexed families like ``patch_ib(500000)%geometry`` resolve in O(1).
    """

    def __init__(self):
        from ..params import REGISTRY

        self._view = REGISTRY.all_params

    def __contains__(self, key):
        return key in self._view

    def __getitem__(self, key):
        return self._view[key].param_type

    def __iter__(self):
        return iter(self._view)

    def __len__(self):
        return len(self._view)


def _registry():
    from ..params import REGISTRY

    return REGISTRY


IGNORE = ["cantera_file", "chemistry"]
ALL = _ParamTypeMapping()
CASE_OPTIMIZATION = [n for n, p in _registry().all_params.items() if p.case_optimization]
SCHEMA = _registry().get_json_schema()

_BASE_NAME_RE = re.compile(r"^([a-zA-Z_][a-zA-Z0-9_]*)")


def _is_param_valid_for_target(param_name: str, target_name: str) -> bool:
    from ..params.namelist_parser import get_target_params

    target_params = get_target_params().get(target_name, set())
    match = _BASE_NAME_RE.match(param_name)
    return match.group(1) in target_params if match else param_name in target_params


class _TargetKeySet:
    """Set-like object for checking param validity for a specific target.

    Supports ``key in obj`` via base-name matching against the Fortran namelist,
    optionally filtering out case-optimization params.
    """

    def __init__(self, target_name: str, filter_case_opt: bool = False):
        self._target_name = target_name
        self._case_opt = set(CASE_OPTIMIZATION) if filter_case_opt else set()

    def __contains__(self, key):
        if key in self._case_opt:
            return False
        return _is_param_valid_for_target(key, self._target_name)


def get_input_dict_keys(target_name: str):
    """Return a set-like object for checking parameter validity for a target.

    Supports ``key in result`` for O(1) checks. Does NOT enumerate indexed family
    members — checks the base name against the Fortran namelist.
    """
    filter_case_opt = ARG("case_optimization", dflt=False) and target_name == "simulation"
    return _TargetKeySet(target_name, filter_case_opt)


def get_validator():
    """Return the cached JSON schema validator."""
    return _registry().get_validator()
