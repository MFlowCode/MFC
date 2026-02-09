"""
Parameter Inventory Export Tool.

Exports all MFC parameters with their types and tags to JSON for analysis.
"""

import re
import json
from pathlib import Path
from typing import Dict, Any

from ..run.case_dicts import ALL
from ..params import REGISTRY
from ..params.schema import ParamType


def get_param_type_name(param_type) -> str:
    """Convert ParamType to string name."""
    if isinstance(param_type, ParamType):
        return param_type.name
    return "UNKNOWN"


def export_parameter_inventory() -> Dict[str, Any]:
    """Export complete parameter inventory with metadata."""
    # Count by type
    by_type = {
        "INT": [],
        "REAL": [],
        "LOG": [],
        "STR": [],
        "ANALYTIC_INT": [],
        "ANALYTIC_REAL": [],
    }

    # Count by tag
    by_tag = {}
    for tag in REGISTRY.get_all_tags():
        by_tag[tag] = []

    inventory = {
        "metadata": {
            "total_parameters": len(ALL),
        },
        "parameters": {},
        "by_type": by_type,
        "by_tag": by_tag,
    }

    for param_name, param_type in sorted(ALL.items()):
        type_name = get_param_type_name(param_type)
        param = REGISTRY.all_params.get(param_name)

        param_info = {
            "type": type_name,
            "tags": sorted(param.tags) if param else [],
        }

        # Detect pattern-based parameters
        if "(" in param_name:
            # Extract pattern (e.g., "patch_icpp(1)%x_centroid" -> "patch_icpp({id})%x_centroid")
            param_pattern = re.sub(r'\((\d+)\)', r'({id})', param_name)
            param_pattern = re.sub(r'\((\d+),\s*(\d+)\)', r'({id1}, {id2})', param_pattern)
            param_info["pattern"] = param_pattern

        inventory["parameters"][param_name] = param_info

        # Categorize by type
        if type_name in by_type:
            by_type[type_name].append(param_name)

        # Categorize by tag
        if param:
            for tag in param.tags:
                if tag in by_tag:
                    by_tag[tag].append(param_name)

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
    print()
    print("By type:")
    for type_name, params in inventory["by_type"].items():
        print(f"  - {type_name}: {len(params)}")
    print()
    print("By feature tag:")
    for tag, params in sorted(inventory["by_tag"].items()):
        if params:
            print(f"  - {tag}: {len(params)}")
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
