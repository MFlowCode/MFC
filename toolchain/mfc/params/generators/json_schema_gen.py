"""
JSON Schema Generator for MFC Case Files.

Generates VS Code / PyCharm compatible JSON Schema for case file auto-completion.
"""
# pylint: disable=import-outside-toplevel

import json
from typing import Dict, Any
from ..schema import ParamType
from ..registry import REGISTRY
from .. import definitions  # noqa: F401  pylint: disable=unused-import


def _param_type_to_json_schema(param_type: ParamType, constraints: Dict = None) -> Dict[str, Any]:
    """Convert ParamType to JSON Schema type definition."""
    base_schemas = {
        ParamType.INT: {"type": "integer"},
        ParamType.REAL: {"type": "number"},
        ParamType.LOG: {"type": "string", "enum": ["T", "F"]},
        ParamType.STR: {"type": "string"},
        # Analytic types allow strings (expressions) or their base type
        ParamType.ANALYTIC_INT: {"oneOf": [{"type": "integer"}, {"type": "string"}]},
        ParamType.ANALYTIC_REAL: {"oneOf": [{"type": "number"}, {"type": "string"}]},
    }

    schema = base_schemas.get(param_type, {"type": "string"}).copy()

    # Add constraints
    if constraints:
        if "choices" in constraints and param_type in (ParamType.INT, ParamType.REAL):
            schema["enum"] = constraints["choices"]
        if "min" in constraints:
            schema["minimum"] = constraints["min"]
        if "max" in constraints:
            schema["maximum"] = constraints["max"]

    return schema


def generate_json_schema(include_descriptions: bool = True) -> Dict[str, Any]:
    """
    Generate JSON Schema for MFC case file parameters.

    Args:
        include_descriptions: Include parameter descriptions in schema

    Returns:
        JSON Schema dict
    """
    from ..descriptions import get_description

    properties = {}
    all_params = []

    for name, param in sorted(REGISTRY.all_params.items()):
        prop_schema = _param_type_to_json_schema(param.param_type, param.constraints)

        if include_descriptions:
            # Get description from descriptions module
            desc = get_description(name)
            if desc:
                prop_schema["description"] = desc

        # Add deprecation notice if applicable
        if param.dependencies and "deprecated" in param.dependencies:
            prop_schema["deprecated"] = True
            prop_schema["deprecationMessage"] = param.dependencies["deprecated"]

        properties[name] = prop_schema
        all_params.append(name)

    schema = {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "$id": "https://mflowcode.github.io/schemas/mfc-case.json",
        "title": "MFC Case File Schema",
        "description": "Schema for MFC (Multi-component Flow Code) simulation case parameters",
        "type": "object",
        "properties": properties,
        "additionalProperties": False,
    }

    return schema


def generate_vscode_settings() -> Dict[str, Any]:
    """Generate VS Code settings snippet for JSON Schema association."""
    return {
        "json.schemas": [
            {
                "fileMatch": ["case.py", "**/case.py"],
                "url": "./mfc-case-schema.json"
            }
        ],
        "yaml.schemas": {
            "./mfc-case-schema.json": ["case.yaml", "**/case.yaml"]
        }
    }


def write_json_schema(output_path: str, include_descriptions: bool = True) -> None:
    """
    Write JSON Schema to file.

    Args:
        output_path: Path to write schema file
        include_descriptions: Include parameter descriptions
    """
    schema = generate_json_schema(include_descriptions)

    with open(output_path, 'w') as f:
        json.dump(schema, f, indent=2)


def get_schema_stats() -> Dict[str, int]:
    """Get statistics about the generated schema."""
    from ..descriptions import get_description

    schema = generate_json_schema(include_descriptions=False)
    props = schema.get("properties", {})

    stats = {
        "total_params": len(props),
        "with_constraints": sum(1 for p in props.values() if "enum" in p or "minimum" in p or "maximum" in p),
        "with_descriptions": sum(1 for name in REGISTRY.all_params if get_description(name)),
    }

    return stats
