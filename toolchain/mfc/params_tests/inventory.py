"""
Parameter Inventory Export Tool.

Exports all MFC parameters with their types and stages to JSON for analysis.
"""

import re
import json
from pathlib import Path
from typing import Dict, Any

from ..run.case_dicts import (
    COMMON, PRE_PROCESS, SIMULATION, POST_PROCESS, ALL, ParamType
)


def get_param_type_name(param_type) -> str:
    """Convert ParamType to string name."""
    if isinstance(param_type, ParamType):
        return param_type.name
    # Handle analytic types (dict)
    if isinstance(param_type, dict):
        if param_type.get("type") == ["integer", "string"]:
            return "ANALYTIC_INT"
        if param_type.get("type") == ["number", "string"]:
            return "ANALYTIC_REAL"
        if param_type.get("type") == "string":
            return "STR"
    return "UNKNOWN"


def export_parameter_inventory() -> Dict[str, Any]:
    """Export complete parameter inventory with metadata."""
    inventory = {
        "metadata": {
            "total_parameters": len(ALL),
            "common_count": len(COMMON),
            "pre_process_count": len(PRE_PROCESS),
            "simulation_count": len(SIMULATION),
            "post_process_count": len(POST_PROCESS),
        },
        "parameters": {},
        "by_stage": {
            "common": [],
            "pre_process_only": [],
            "simulation_only": [],
            "post_process_only": [],
        },
        "by_type": {
            "INT": [],
            "REAL": [],
            "LOG": [],
            "STR": [],
            "ANALYTIC_INT": [],
            "ANALYTIC_REAL": [],
        }
    }

    # Categorize parameters
    common_keys = set(COMMON.keys())
    pre_only = set(PRE_PROCESS.keys()) - common_keys
    sim_only = set(SIMULATION.keys()) - common_keys
    post_only = set(POST_PROCESS.keys()) - common_keys

    for param_name, param_type in sorted(ALL.items()):
        type_name = get_param_type_name(param_type)

        # Determine which stages this parameter is valid for
        stages = []
        if param_name in COMMON:
            stages = ["common", "pre_process", "simulation", "post_process"]
        else:
            if param_name in PRE_PROCESS:
                stages.append("pre_process")
            if param_name in SIMULATION:
                stages.append("simulation")
            if param_name in POST_PROCESS:
                stages.append("post_process")

        param_info = {
            "type": type_name,
            "stages": stages,
        }

        # Detect pattern-based parameters
        if "(" in param_name:
            # Extract pattern (e.g., "patch_icpp(1)%x_centroid" -> "patch_icpp({id})%x_centroid")
            param_pattern = re.sub(r'\((\d+)\)', r'({id})', param_name)
            param_pattern = re.sub(r'\((\d+),\s*(\d+)\)', r'({id1}, {id2})', param_pattern)
            param_info["pattern"] = param_pattern

        inventory["parameters"][param_name] = param_info

        # Categorize by type
        if type_name in inventory["by_type"]:
            inventory["by_type"][type_name].append(param_name)

    # Categorize by stage
    inventory["by_stage"]["common"] = sorted(common_keys)
    inventory["by_stage"]["pre_process_only"] = sorted(pre_only)
    inventory["by_stage"]["simulation_only"] = sorted(sim_only)
    inventory["by_stage"]["post_process_only"] = sorted(post_only)

    return inventory


def export_parameter_patterns() -> Dict[str, Any]:
    """Extract unique parameter patterns (for dynamic parameters)."""
    patterns = {}
    for param_name, param_type in ALL.items():
        if "(" not in param_name:
            continue

        # Normalize the pattern
        normalized = re.sub(r'\((\d+)\)', r'({N})', param_name)
        normalized = re.sub(r'\((\d+),\s*(\d+)\)', r'({N}, {M})', normalized)

        if normalized not in patterns:
            patterns[normalized] = {
                "examples": [],
                "type": get_param_type_name(param_type),
                "count": 0
            }
        patterns[normalized]["examples"].append(param_name)
        patterns[normalized]["count"] += 1

    # Trim examples to max 3
    for pattern_data in patterns.values():
        pattern_data["examples"] = pattern_data["examples"][:3]

    return patterns


def save_inventory(output_path: Path = None):
    """Save parameter inventory to JSON file."""
    if output_path is None:
        output_path = Path(__file__).parent / "param_inventory.json"

    inventory = export_parameter_inventory()
    inventory["patterns"] = export_parameter_patterns()

    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(inventory, f, indent=2)

    return output_path


def print_inventory_summary():
    """Print a summary of the parameter inventory."""
    inventory = export_parameter_inventory()
    patterns = export_parameter_patterns()

    print("=" * 60)
    print("MFC Parameter Inventory Summary")
    print("=" * 60)
    print(f"Total parameters: {inventory['metadata']['total_parameters']}")
    print(f"  - Common:       {inventory['metadata']['common_count']}")
    print(f"  - Pre-process:  {inventory['metadata']['pre_process_count']}")
    print(f"  - Simulation:   {inventory['metadata']['simulation_count']}")
    print(f"  - Post-process: {inventory['metadata']['post_process_count']}")
    print()
    print("By type:")
    for type_name, params in inventory["by_type"].items():
        print(f"  - {type_name}: {len(params)}")
    print()
    print(f"Dynamic parameter patterns: {len(patterns)}")
    print("Top patterns:")
    sorted_patterns = sorted(patterns.items(), key=lambda x: -x[1]["count"])[:10]
    for pattern, info in sorted_patterns:
        print(f"  - {pattern}: {info['count']} instances ({info['type']})")


if __name__ == "__main__":
    print_inventory_summary()
    path = save_inventory()
    print(f"\nInventory saved to: {path}")
