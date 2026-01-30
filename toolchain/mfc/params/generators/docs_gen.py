"""
Parameter Documentation Generator.

Generates markdown documentation for all MFC case parameters,
organized by family with descriptions, types, and constraints.
"""

from typing import Dict, List, Tuple
from collections import defaultdict
import re

from ..schema import ParamType
from ..registry import REGISTRY
from ..descriptions import get_description
from .. import definitions  # noqa: F401  pylint: disable=unused-import


def _get_family(name: str) -> str:
    """Extract family name from parameter (e.g., 'patch_icpp' from 'patch_icpp(1)%vel(1)')."""
    # Handle indexed parameters
    match = re.match(r'^([a-z_]+)', name)
    if match:
        base = match.group(1)
        # Check if it's a known family pattern
        if any(name.startswith(f"{base}(") or name.startswith(f"{base}%") for _ in [1]):
            return base
    return "general"


def _type_to_str(param_type: ParamType) -> str:
    """Convert ParamType to readable string."""
    return {
        ParamType.INT: "Integer",
        ParamType.REAL: "Real",
        ParamType.LOG: "Logical (T/F)",
        ParamType.STR: "String",
        ParamType.ANALYTIC_INT: "Integer or Expression",
        ParamType.ANALYTIC_REAL: "Real or Expression",
    }.get(param_type, str(param_type))


def _format_constraints(param) -> str:
    """Format constraints as readable string."""
    if not param.constraints:
        return ""

    parts = []
    c = param.constraints
    if "choices" in c:
        parts.append(f"Values: {c['choices']}")
    if "min" in c:
        parts.append(f"Min: {c['min']}")
    if "max" in c:
        parts.append(f"Max: {c['max']}")

    return ", ".join(parts)


def generate_parameter_docs() -> str:
    """Generate markdown documentation for all parameters."""
    lines = [
        "# MFC Case Parameters Reference",
        "",
        "> **Auto-generated** from parameter registry",
        "> ",
        "> Regenerate with: `./mfc.sh generate --json-schema`",
        "",
        "## Overview",
        "",
        f"MFC supports **{len(REGISTRY.all_params):,}** case parameters organized into families.",
        "",
        "**Quick search:** Use `./mfc.sh params <query>` to search parameters from the command line.",
        "",
        "## Parameter Families",
        "",
    ]

    # Group parameters by family
    families: Dict[str, List[Tuple[str, any]]] = defaultdict(list)
    for name, param in sorted(REGISTRY.all_params.items()):
        family = _get_family(name)
        families[family].append((name, param))

    # Sort families by size (largest first), but put "general" last
    sorted_families = sorted(
        families.items(),
        key=lambda x: (x[0] == "general", -len(x[1]), x[0])
    )

    # Table of contents
    lines.append("| Family | Count | Description |")
    lines.append("|--------|-------|-------------|")

    family_descriptions = {
        "general": "Core simulation parameters (grid, time, model, etc.)",
        "patch_icpp": "Initial condition patch parameters",
        "patch_ib": "Immersed boundary patch parameters",
        "patch_bc": "Boundary condition patch parameters",
        "fluid_pp": "Fluid material properties",
        "acoustic": "Acoustic source parameters",
        "bc_x": "X-direction boundary conditions",
        "bc_y": "Y-direction boundary conditions",
        "bc_z": "Z-direction boundary conditions",
        "probe": "Probe/monitoring point parameters",
        "integral": "Integral region parameters",
        "simplex_params": "Simplex noise perturbation parameters",
        "chem_wrt_Y": "Chemistry species output parameters",
        "bub_pp": "Bubble property parameters",
        "lag_params": "Lagrangian particle parameters",
    }

    for family, params in sorted_families:
        desc = family_descriptions.get(family, "")
        anchor = family.replace("_", "-")
        lines.append(f"| [{family}](#{anchor}) | {len(params)} | {desc} |")

    lines.append("")
    lines.append("---")
    lines.append("")

    # Document each family
    for family, params in sorted_families:
        lines.append(f"## {family}")
        lines.append("")

        desc = family_descriptions.get(family, "")
        if desc:
            lines.append(f"*{desc}*")
            lines.append("")

        lines.append(f"**{len(params)} parameters**")
        lines.append("")

        # For large families, show a collapsed pattern
        if len(params) > 20:
            # Group by pattern
            patterns: Dict[str, List[str]] = defaultdict(list)
            for name, _ in params:
                # Extract pattern (replace indices with N)
                pattern = re.sub(r'\(\d+\)', '(N)', name)
                patterns[pattern].append(name)

            lines.append("### Patterns")
            lines.append("")
            lines.append("| Pattern | Example | Description |")
            lines.append("|---------|---------|-------------|")

            for pattern, examples in sorted(patterns.items()):
                example = examples[0]
                desc = get_description(example) or ""
                # Truncate long descriptions
                if len(desc) > 60:
                    desc = desc[:57] + "..."
                lines.append(f"| `{pattern}` | `{example}` | {desc} |")

            lines.append("")
        else:
            # Show full table for small families
            lines.append("| Parameter | Type | Description |")
            lines.append("|-----------|------|-------------|")

            for name, param in params:
                type_str = _type_to_str(param.param_type)
                desc = get_description(name) or ""
                constraints = _format_constraints(param)
                if constraints:
                    desc = f"{desc} ({constraints})" if desc else constraints
                # Truncate long descriptions
                if len(desc) > 80:
                    desc = desc[:77] + "..."
                lines.append(f"| `{name}` | {type_str} | {desc} |")

            lines.append("")

        lines.append("---")
        lines.append("")

    # Add footer
    lines.extend([
        "## Command Line Reference",
        "",
        "Search parameters using the CLI:",
        "",
        "```bash",
        "# Search for parameters",
        "./mfc.sh params weno",
        "",
        "# Show parameter descriptions",
        "./mfc.sh params weno -d",
        "",
        "# List all families",
        "./mfc.sh params -f",
        "",
        "# Filter by type",
        "./mfc.sh params -t real weno",
        "```",
        "",
    ])

    return "\n".join(lines)


def write_parameter_docs(output_path: str) -> int:
    """Write parameter documentation to file.

    Returns:
        Number of parameters documented
    """
    content = generate_parameter_docs()
    with open(output_path, 'w') as f:
        f.write(content)
    return len(REGISTRY.all_params)
