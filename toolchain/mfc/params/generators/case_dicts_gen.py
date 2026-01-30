"""
Case Dictionary Schema Generator.

Generates case_dicts.py compatible type schemas from the parameter registry.
This allows maintaining parameter types in a single place while still producing
the format expected by MFC's case loading system.
"""

from typing import Dict, List, TYPE_CHECKING
from ..schema import Stage
from ..schema import ParamType as RegistryParamType
from ..registry import ParamRegistry

if TYPE_CHECKING:
    from mfc.run.case_dicts import ParamType as CaseDictsParamType


def generate_case_dicts_schema(registry: ParamRegistry, stage: Stage) -> Dict[str, str]:
    """
    Generate a case_dicts.py compatible schema dict for a stage.

    Args:
        registry: The parameter registry
        stage: The stage to generate schema for

    Returns:
        Dict mapping parameter names to type tags (e.g., {"m": "int", "dt": "real"})
    """
    params = registry.get_params_by_stage(stage)
    return {name: param.type_tag for name, param in sorted(params.items())}


def get_stage_dict(stage: Stage, include_common: bool = True) -> Dict[str, "CaseDictsParamType"]:
    """
    Generate a case_dicts.py compatible dict for a stage.

    This returns a dict mapping parameter names to case_dicts.ParamType enum values,
    suitable for direct use in case_dicts.py.

    Args:
        stage: The stage to generate dict for (COMMON, PRE_PROCESS, etc.)
        include_common: If True and stage != COMMON, include COMMON params too
                        (matches case_dicts.py's PRE_PROCESS = COMMON.copy() pattern)

    Returns:
        Dict mapping parameter names to ParamType enum values
    """
    # Import here to avoid circular imports
    from mfc.run.case_dicts import ParamType as CDParamType
    from ..definitions import load_all_definitions

    # Ensure definitions are loaded
    load_all_definitions()

    from ..registry import REGISTRY

    # Map registry ParamType to case_dicts ParamType
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
        # Check if this param applies to the requested stage
        # Include COMMON params for non-COMMON stages if include_common=True
        if stage in param.stages:
            result[name] = type_map[param.param_type]
        elif include_common and stage != Stage.COMMON and Stage.COMMON in param.stages:
            result[name] = type_map[param.param_type]

    return result


def get_all_stage_dicts() -> Dict[str, Dict[str, "CaseDictsParamType"]]:
    """
    Generate all stage dicts for case_dicts.py.

    Returns:
        Dict with keys 'COMMON', 'PRE_PROCESS', 'SIMULATION', 'POST_PROCESS'
    """
    return {
        'COMMON': get_stage_dict(Stage.COMMON),
        'PRE_PROCESS': get_stage_dict(Stage.PRE_PROCESS),
        'SIMULATION': get_stage_dict(Stage.SIMULATION),
        'POST_PROCESS': get_stage_dict(Stage.POST_PROCESS),
    }


def generate_case_dicts_code(registry: ParamRegistry) -> str:
    """
    Generate the full case_dicts.py file content.

    This produces Python code that defines COMMON, PRE_PROCESS, SIMULATION,
    POST_PROCESS, and ALL dictionaries matching the existing case_dicts.py format.

    Args:
        registry: The parameter registry with all parameters registered

    Returns:
        Python source code as a string
    """
    lines = [
        '"""',
        'MFC Case Parameter Type Definitions.',
        '',
        'AUTO-GENERATED from mfc.params registry - Do not edit manually.',
        'To regenerate: python -m mfc.params.generators.case_dicts_gen',
        '"""',
        '',
    ]

    # Generate each stage dict
    for stage in [Stage.COMMON, Stage.PRE_PROCESS, Stage.SIMULATION, Stage.POST_PROCESS]:
        stage_name = stage.name
        params = registry.get_params_by_stage(stage)

        # Only include params that are specifically for this stage (not inherited)
        if stage != Stage.COMMON:
            common_params = registry.get_params_by_stage(Stage.COMMON)
            params = {k: v for k, v in params.items() if k not in common_params}

        lines.append(f'{stage_name} = {{')

        # Group by category if available
        by_category: Dict[str, List[tuple]] = {}
        for name, param in sorted(params.items()):
            cat = param.category or "other"
            if cat not in by_category:
                by_category[cat] = []
            by_category[cat].append((name, param.type_tag))

        for category, items in sorted(by_category.items()):
            if len(by_category) > 1:
                lines.append(f'    # {category}')
            for name, type_tag in sorted(items):
                lines.append(f'    "{name}": "{type_tag}",')

        lines.append('}')
        lines.append('')

    # Generate ALL dict
    lines.append('# Combined dictionary of all parameters')
    lines.append('ALL = {')
    lines.append('    **COMMON,')
    lines.append('    **PRE_PROCESS,')
    lines.append('    **SIMULATION,')
    lines.append('    **POST_PROCESS,')
    lines.append('}')

    return '\n'.join(lines)


def compare_with_existing(registry: ParamRegistry, existing_path: str) -> Dict[str, List[str]]:
    """
    Compare registry parameters with existing case_dicts.py.

    Args:
        registry: The parameter registry
        existing_path: Path to existing case_dicts.py

    Returns:
        Dict with keys 'missing_in_registry', 'missing_in_existing', 'type_mismatches'
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location("case_dicts", existing_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    existing_all = getattr(module, 'ALL', {})
    registry_all = {}

    for stage in Stage:
        for name, param in registry.get_params_by_stage(stage).items():
            registry_all[name] = param.type_tag

    missing_in_registry = [k for k in existing_all if k not in registry_all]
    missing_in_existing = [k for k in registry_all if k not in existing_all]
    type_mismatches = [
        f"{k}: registry={registry_all[k]}, existing={existing_all[k]}"
        for k in registry_all
        if k in existing_all and registry_all[k] != existing_all[k]
    ]

    return {
        'missing_in_registry': missing_in_registry,
        'missing_in_existing': missing_in_existing,
        'type_mismatches': type_mismatches,
    }


if __name__ == "__main__":
    from ..registry import REGISTRY
    # Import definitions to populate registry
    from ..definitions import load_all_definitions
    load_all_definitions()

    print(generate_case_dicts_code(REGISTRY))
