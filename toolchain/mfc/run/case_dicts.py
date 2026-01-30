"""
MFC Case Parameter Type Definitions.

This module defines parameter types for MFC case files. The parameter definitions
are now sourced from the central registry (mfc.params), providing a single source
of truth for all ~3,300 parameters.

The ParamType enum defines JSON schema types for validation.
"""

import fastjsonschema

from enum import Enum
from ..state import ARG
from functools import cache


class ParamType(Enum):
    """
    Parameter types for JSON schema validation.

    Each type maps to a JSON schema fragment used by fastjsonschema
    to validate case file parameters.
    """
    INT = {"type": "integer"}
    REAL = {"type": "number"}
    LOG = {"enum": ["T", "F"]}
    STR = {"type": "string"}

    _ANALYTIC_INT = {"type": ["integer", "string"]}
    _ANALYTIC_REAL = {"type": ["number", "string"]}

    def analytic(self):
        """Return the analytic version of this type (allows string expressions)."""
        if self == self.INT:
            return self._ANALYTIC_INT
        if self == self.REAL:
            return self._ANALYTIC_REAL
        return self.STR


def _load_from_registry():
    """
    Load parameter definitions from the central registry.

    Returns:
        Tuple of (COMMON, PRE_PROCESS, SIMULATION, POST_PROCESS) dicts
    """
    from ..params.generators.case_dicts_gen import get_stage_dict
    from ..params.schema import Stage

    return (
        get_stage_dict(Stage.COMMON, include_common=False),
        get_stage_dict(Stage.PRE_PROCESS, include_common=True),
        get_stage_dict(Stage.SIMULATION, include_common=True),
        get_stage_dict(Stage.POST_PROCESS, include_common=True),
    )


# Load parameter definitions from registry
COMMON, PRE_PROCESS, SIMULATION, POST_PROCESS = _load_from_registry()

# Parameters to ignore during certain operations
IGNORE = ["cantera_file", "chemistry"]

# Build combined ALL dict
ALL = COMMON.copy()
ALL.update(PRE_PROCESS)
ALL.update(SIMULATION)
ALL.update(POST_PROCESS)


def _get_case_optimization_params():
    """Get params that can be hard-coded for GPU optimization from registry."""
    from ..params import REGISTRY
    return [name for name, param in REGISTRY.all_params.items() if param.case_optimization]


# Parameters that can be hard-coded for GPU case optimization (from registry)
CASE_OPTIMIZATION = _get_case_optimization_params()

# Build JSON schema for validation
_properties = {k: v.value for k, v in ALL.items()}

SCHEMA = {
    "type": "object",
    "properties": _properties,
    "additionalProperties": False
}


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
        return result

    return [x for x in result if x not in CASE_OPTIMIZATION]


@cache
def get_validator():
    """Get the cached JSON schema validator."""
    return fastjsonschema.compile(SCHEMA)
