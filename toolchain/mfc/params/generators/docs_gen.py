"""
Parameter Documentation Generator.

Generates markdown documentation for all MFC case parameters,
organized by family with descriptions, types, and constraints.
"""

from typing import Any, Dict, List, Tuple
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


def _escape_percent(s: str) -> str:
    """Escape % for Doxygen (% is a special character, use %% to get literal %)."""
    return s.replace('%', '%%')


def _parse_paren_content(name: str, start: int) -> Tuple[str, int]:
    """Parse content within parentheses, return (content, end_index) or ('', -1) if invalid."""
    j = start + 1
    paren_content = []
    while j < len(name) and name[j] != ')':
        paren_content.append(name[j])
        j += 1
    if j < len(name):
        return ''.join(paren_content), j
    return '', -1


def _collapse_indices(name: str) -> str:
    """
    Collapse numeric indices to placeholders for pattern grouping.

    Examples:
        patch_icpp(1)%vel(2) -> patch_icpp(N)%vel(M)
        simplex_params%perturb_dens_offset(1, 2) -> simplex_params%perturb_dens_offset(N, M)
        bc_x%vel_in(1) -> bc_x%vel_in(N)
    """
    placeholders = ['N', 'M', 'K', 'L', 'P', 'Q']
    placeholder_idx = 0
    result = []
    i = 0

    while i < len(name):
        if name[i] != '(':
            result.append(name[i])
            i += 1
            continue

        # Found opening paren, look for indices
        content, end_idx = _parse_paren_content(name, i)
        if end_idx == -1:
            result.append(name[i])
            i += 1
            continue

        # Check if content is numeric indices (possibly comma-separated)
        parts = [p.strip() for p in content.split(',')]
        if not all(p.isdigit() for p in parts):
            result.append(name[i])
            i += 1
            continue

        # Replace each index with a placeholder
        new_parts = []
        for _ in parts:
            ph = placeholders[placeholder_idx] if placeholder_idx < len(placeholders) else '?'
            new_parts.append(ph)
            placeholder_idx += 1
        result.append('(' + ', '.join(new_parts) + ')')
        i = end_idx + 1

    return ''.join(result)


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


def generate_parameter_docs() -> str:  # pylint: disable=too-many-locals,too-many-statements
    """Generate markdown documentation for all parameters."""
    lines = [
        "@page parameters Case Parameters Reference",
        "",
        "# Case Parameters Reference",
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
    families: Dict[str, List[Tuple[str, Any]]] = defaultdict(list)
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
        # Post-processing output flags
        "alpha_rho_wrt": "Partial density output flags",
        "alpha_rho_e_wrt": "Partial density-energy output flags",
        "alpha_wrt": "Volume fraction output flags",
        "kappa_wrt": "Curvature output flags",
        "schlieren_alpha": "Numerical schlieren coefficients",
        "mom_wrt": "Momentum output flags",
        "vel_wrt": "Velocity output flags",
        "flux_wrt": "Flux output flags",
        "omega_wrt": "Vorticity output flags",
        # Domain and output regions
        "x_domain": "X-direction domain parameters",
        "y_domain": "Y-direction domain parameters",
        "z_domain": "Z-direction domain parameters",
        "x_output": "X-direction output region",
        "y_output": "Y-direction output region",
        "z_output": "Z-direction output region",
        # Other
        "fluid_rho": "Fluid reference densities",
        "chem_params": "Chemistry model parameters",
    }

    for family, params in sorted_families:
        desc = family_descriptions.get(family, "")
        # Use family name directly as anchor (GitHub keeps underscores in heading anchors)
        lines.append(f"| [{family}](#{family}) | {len(params)} | {desc} |")

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

        # Group by pattern (collapse indices to N, M, etc.)
        patterns: Dict[str, List[str]] = defaultdict(list)
        for name, _ in params:
            pattern = _collapse_indices(name)
            patterns[pattern].append(name)

        # Use pattern view if it reduces rows, otherwise show full table
        if len(patterns) < len(params):
            # Pattern view - shows collapsed patterns
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
                # Escape % for Doxygen
                pattern_escaped = _escape_percent(pattern)
                example_escaped = _escape_percent(example)
                lines.append(f"| `{pattern_escaped}` | `{example_escaped}` | {desc} |")

            lines.append("")
        else:
            # Full table - no patterns to collapse
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
                # Escape % for Doxygen
                name_escaped = _escape_percent(name)
                lines.append(f"| `{name_escaped}` | {type_str} | {desc} |")

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
