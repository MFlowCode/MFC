"""
Case Dictionary Schema Generator (Minimal).

Generates case_dicts.py compatible type schemas from the parameter registry.
Includes verification to ensure type systems stay synchronized.
"""
# pylint: disable=import-outside-toplevel

from typing import Dict, List, TYPE_CHECKING
from ..schema import Stage, ParamType as RegistryParamType

if TYPE_CHECKING:
    from ...run.case_dicts import ParamType as CaseDictsParamType


def _get_type_map() -> Dict[RegistryParamType, "CaseDictsParamType"]:
    """
    Get the mapping from registry ParamType to case_dicts ParamType.

    Returns:
        Dict mapping each RegistryParamType to corresponding CaseDictsParamType.

    Raises:
        RuntimeError: If any registry type is not mapped.
    """
    from ...run.case_dicts import ParamType as CDParamType

    type_map = {
        RegistryParamType.INT: CDParamType.INT,
        RegistryParamType.REAL: CDParamType.REAL,
        RegistryParamType.LOG: CDParamType.LOG,
        RegistryParamType.STR: CDParamType.STR,
        RegistryParamType.ANALYTIC_INT: CDParamType._ANALYTIC_INT,  # pylint: disable=protected-access
        RegistryParamType.ANALYTIC_REAL: CDParamType._ANALYTIC_REAL,  # pylint: disable=protected-access
    }

    # Verify all registry types are mapped
    missing = set(RegistryParamType) - set(type_map.keys())
    if missing:
        raise RuntimeError(
            f"Type system sync error: Registry ParamTypes not mapped to case_dicts: {missing}. "
            "Update _get_type_map() in case_dicts_gen.py to include these types."
        )

    return type_map


def verify_type_systems() -> List[str]:
    """
    Verify that the two type systems (params.ParamType and case_dicts.ParamType) are synchronized.

    This check ensures:
    1. All registry ParamTypes have a mapping to case_dicts ParamTypes
    2. The mapping produces valid case_dicts types

    Returns:
        List of error messages (empty if all checks pass).
    """
    errors = []

    try:
        type_map = _get_type_map()
    except RuntimeError as e:
        errors.append(str(e))
        return errors

    # Verify each mapping produces a valid case_dicts ParamType
    from ...run.case_dicts import ParamType as CDParamType

    valid_cd_types = set(CDParamType)
    for reg_type, cd_type in type_map.items():
        if cd_type not in valid_cd_types:
            errors.append(
                f"Type mapping error: {reg_type} maps to invalid case_dicts type {cd_type}"
            )

    return errors


def get_stage_dict(stage: Stage, include_common: bool = True) -> Dict[str, "CaseDictsParamType"]:
    """
    Generate a case_dicts.py compatible dict for a stage.

    Args:
        stage: The stage to generate dict for
        include_common: If True and stage != COMMON, include COMMON params too

    Returns:
        Dict mapping parameter names to ParamType enum values

    Raises:
        RuntimeError: If type systems are out of sync.
    """
    from .. import definitions  # noqa: F401  pylint: disable=unused-import
    from ..registry import REGISTRY

    type_map = _get_type_map()

    result = {}
    for name, param in REGISTRY.all_params.items():
        if stage in param.stages:
            result[name] = type_map[param.param_type]
        elif include_common and stage != Stage.COMMON and Stage.COMMON in param.stages:
            result[name] = type_map[param.param_type]

    return result
