"""
Case Dictionary Schema Generator (Minimal).

Generates case_dicts.py compatible type schemas from the parameter registry.
"""

from typing import Dict, TYPE_CHECKING
from ..schema import Stage, ParamType as RegistryParamType

if TYPE_CHECKING:
    from mfc.run.case_dicts import ParamType as CaseDictsParamType


def get_stage_dict(stage: Stage, include_common: bool = True) -> Dict[str, "CaseDictsParamType"]:
    """
    Generate a case_dicts.py compatible dict for a stage.

    Args:
        stage: The stage to generate dict for
        include_common: If True and stage != COMMON, include COMMON params too

    Returns:
        Dict mapping parameter names to ParamType enum values
    """
    from mfc.run.case_dicts import ParamType as CDParamType
    from .. import definitions  # noqa: F401  pylint: disable=unused-import
    from ..registry import REGISTRY

    type_map = {
        RegistryParamType.INT: CDParamType.INT,
        RegistryParamType.REAL: CDParamType.REAL,
        RegistryParamType.LOG: CDParamType.LOG,
        RegistryParamType.STR: CDParamType.STR,
        RegistryParamType.ANALYTIC_INT: CDParamType._ANALYTIC_INT,
        RegistryParamType.ANALYTIC_REAL: CDParamType._ANALYTIC_REAL,
    }

    result = {}
    for name, param in REGISTRY.all_params.items():
        if stage in param.stages:
            result[name] = type_map[param.param_type]
        elif include_common and stage != Stage.COMMON and Stage.COMMON in param.stages:
            result[name] = type_map[param.param_type]

    return result
