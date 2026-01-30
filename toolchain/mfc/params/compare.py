"""
Comparison tool for registry vs existing case_dicts.py.

Compares parameter names and types between the new registry and
the existing case_dicts.py to identify gaps and mismatches.
"""

from typing import Dict, List, Set, Tuple
from pathlib import Path

from .schema import ParamType, Stage
from .registry import REGISTRY


def get_case_dicts_params() -> Dict[str, Dict[str, str]]:
    """Load parameters from existing case_dicts.py."""
    from mfc.run.case_dicts import COMMON, PRE_PROCESS, SIMULATION, POST_PROCESS, ParamType as CDParamType

    def type_to_str(pt) -> str:
        """Convert case_dicts ParamType to string."""
        if pt == CDParamType.INT:
            return "int"
        elif pt == CDParamType.REAL:
            return "real"
        elif pt == CDParamType.LOG:
            return "log"
        elif pt == CDParamType.STR:
            return "str"
        elif pt == CDParamType._ANALYTIC_INT:
            return "analytic:int"
        elif pt == CDParamType._ANALYTIC_REAL:
            return "analytic:real"
        else:
            # Handle .analytic() results
            val = pt.get("type") if isinstance(pt, dict) else None
            if val == ["integer", "string"]:
                return "analytic:int"
            elif val == ["number", "string"]:
                return "analytic:real"
            return str(pt)

    return {
        "COMMON": {k: type_to_str(v) for k, v in COMMON.items()},
        "PRE_PROCESS": {k: type_to_str(v) for k, v in PRE_PROCESS.items()},
        "SIMULATION": {k: type_to_str(v) for k, v in SIMULATION.items()},
        "POST_PROCESS": {k: type_to_str(v) for k, v in POST_PROCESS.items()},
    }


def get_registry_params() -> Dict[str, Dict[str, str]]:
    """Get parameters from registry organized by stage."""
    result = {
        "COMMON": {},
        "PRE_PROCESS": {},
        "SIMULATION": {},
        "POST_PROCESS": {},
    }

    for name, param in REGISTRY.all_params.items():
        type_str = param.type_tag
        for stage in param.stages:
            stage_name = stage.name
            if stage_name in result:
                result[stage_name][name] = type_str

    return result


def compare_params() -> Dict[str, any]:
    """
    Compare registry parameters with case_dicts.py.

    Returns dict with:
    - missing_in_registry: params in case_dicts but not registry
    - extra_in_registry: params in registry but not case_dicts
    - type_mismatches: params with different types
    - matching: params that match exactly
    """
    cd_params = get_case_dicts_params()
    reg_params = get_registry_params()

    # Get all unique params from case_dicts
    cd_all = set()
    for stage_params in cd_params.values():
        cd_all.update(stage_params.keys())

    # Get all unique params from registry
    reg_all = set(REGISTRY.all_params.keys())

    missing = cd_all - reg_all
    extra = reg_all - cd_all
    common = cd_all & reg_all

    # Check type mismatches for common params
    type_mismatches = []
    matching = []

    for name in common:
        # Get type from case_dicts (use PRE_PROCESS as primary)
        cd_type = None
        for stage in ["PRE_PROCESS", "SIMULATION", "POST_PROCESS", "COMMON"]:
            if name in cd_params[stage]:
                cd_type = cd_params[stage][name]
                break

        # Get type from registry
        reg_type = REGISTRY.get(name).type_tag

        if cd_type != reg_type:
            type_mismatches.append((name, cd_type, reg_type))
        else:
            matching.append(name)

    return {
        "missing_in_registry": sorted(missing),
        "extra_in_registry": sorted(extra),
        "type_mismatches": type_mismatches,
        "matching": sorted(matching),
        "summary": {
            "case_dicts_total": len(cd_all),
            "registry_total": len(reg_all),
            "missing": len(missing),
            "extra": len(extra),
            "mismatches": len(type_mismatches),
            "matching": len(matching),
        }
    }


def print_comparison_report():
    """Print a detailed comparison report."""
    from .definitions import load_all_definitions
    load_all_definitions()

    result = compare_params()
    summary = result["summary"]

    print("=" * 70)
    print("Registry vs case_dicts.py Comparison")
    print("=" * 70)

    print(f"\ncase_dicts.py: {summary['case_dicts_total']} parameters")
    print(f"Registry:      {summary['registry_total']} parameters")

    print(f"\n  Matching:              {summary['matching']}")
    print(f"  Missing in registry:   {summary['missing']}")
    print(f"  Extra in registry:     {summary['extra']}")
    print(f"  Type mismatches:       {summary['mismatches']}")

    coverage = summary['matching'] / summary['case_dicts_total'] * 100
    print(f"\n  Coverage: {coverage:.1f}%")

    if result["type_mismatches"]:
        print("\n" + "-" * 70)
        print("TYPE MISMATCHES:")
        print("-" * 70)
        for name, cd_type, reg_type in result["type_mismatches"][:20]:
            print(f"  {name}")
            print(f"    case_dicts: {cd_type}")
            print(f"    registry:   {reg_type}")

    if result["missing_in_registry"]:
        print("\n" + "-" * 70)
        print("MISSING IN REGISTRY (sample of 30):")
        print("-" * 70)

        # Group by prefix
        by_prefix = {}
        for name in result["missing_in_registry"]:
            if "%" in name:
                prefix = name.split("%")[0]
            elif "(" in name:
                prefix = name.split("(")[0]
            else:
                prefix = name
            if prefix not in by_prefix:
                by_prefix[prefix] = []
            by_prefix[prefix].append(name)

        for prefix, names in sorted(by_prefix.items(), key=lambda x: -len(x[1]))[:15]:
            print(f"  {prefix}: {len(names)} params")
            for name in names[:3]:
                print(f"    - {name}")
            if len(names) > 3:
                print(f"    ... and {len(names) - 3} more")


def categorize_missing() -> Dict[str, List[str]]:
    """Categorize missing parameters by family."""
    from .definitions import load_all_definitions
    load_all_definitions()

    result = compare_params()
    missing = result["missing_in_registry"]

    categories = {
        "patch_ib": [],
        "patch_bc": [],
        "acoustic": [],
        "probe": [],
        "integral": [],
        "bc_*": [],
        "simplex_params": [],
        "lag_params": [],
        "chem_params": [],
        "fluid_rho": [],
        "post_process_wrt": [],
        "other": [],
    }

    for name in missing:
        categorized = False
        for cat in categories:
            if cat.endswith("*"):
                prefix = cat[:-1]
                if name.startswith(prefix):
                    categories[cat].append(name)
                    categorized = True
                    break
            elif name.startswith(cat) or cat in name:
                categories[cat].append(name)
                categorized = True
                break
        if not categorized:
            categories["other"].append(name)

    return {k: v for k, v in categories.items() if v}


if __name__ == "__main__":
    print_comparison_report()
