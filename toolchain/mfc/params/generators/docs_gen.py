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
from ..ast_analyzer import analyze_case_validator, classify_message
from .. import definitions  # noqa: F401  pylint: disable=unused-import


def _get_family(name: str) -> str:
    """Extract family name from parameter (e.g., 'patch_icpp' from 'patch_icpp(1)%vel(1)')."""
    # Handle indexed parameters
    match = re.match(r'^([a-zA-Z_]+)', name)
    if match:
        base = match.group(1)
        if name.startswith(f"{base}(") or name.startswith(f"{base}%"):
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
    """Format constraints as readable string, skipping choices already shown via value_labels."""
    if not param.constraints:
        return ""

    parts = []
    c = param.constraints
    if "choices" in c and "value_labels" not in c:
        parts.append(f"Values: {c['choices']}")
    if "min" in c:
        parts.append(f"Min: {c['min']}")
    if "max" in c:
        parts.append(f"Max: {c['max']}")

    return ", ".join(parts)


def _build_param_name_pattern():
    """Build a regex pattern that matches known parameter names at word boundaries.

    Uses longest-match-first to avoid partial matches (e.g., 'model_eqns' before 'model').
    Only matches names that look like identifiers (avoids matching 'm' inside 'must').
    """
    all_names = sorted(REGISTRY.all_params.keys(), key=len, reverse=True)
    # Only include names >= 2 chars to avoid false positives with single-letter params
    # and names that are simple identifiers (no % or parens, which need escaping)
    safe_names = [n for n in all_names if len(n) >= 2 and re.match(r'^[a-zA-Z_]\w*$', n)]
    if not safe_names:
        return None
    pattern = r'\b(' + '|'.join(re.escape(n) for n in safe_names) + r')\b'
    return re.compile(pattern)


# Matches compound param names like bub_pp%mu_g, fluid_pp(1)%Re(1), x_output%beg
_COMPOUND_NAME_RE = re.compile(r'\b\w+(?:\([^)]*\))?(?:%\w+(?:\([^)]*\))?)+')


def _backtick_params(msg: str, pattern) -> str:
    """Wrap parameter names in backticks for markdown rendering.

    Handles three cases in order:
    1. Compound names with % (e.g. bub_pp%mu_g, x_output%beg)
    2. Known registry param names (e.g. model_eqns, weno_order)
    3. Snake_case identifiers not in registry (e.g. cluster_type, smooth_type)
    """
    # 1. Wrap compound names (word%word patterns) — must come first
    msg = _COMPOUND_NAME_RE.sub(lambda m: f'`{m.group(0)}`', msg)

    # 2. Wrap known simple param names, only outside existing backtick spans
    if pattern is not None:
        parts = msg.split('`')
        for i in range(0, len(parts), 2):
            parts[i] = pattern.sub(r'`\1`', parts[i])
        msg = '`'.join(parts)

    # 3. Wrap remaining snake_case identifiers (at least one underscore)
    parts = msg.split('`')
    for i in range(0, len(parts), 2):
        parts[i] = re.sub(r'\b([a-z]\w*_\w+)\b', r'`\1`', parts[i])
    msg = '`'.join(parts)

    return msg


def _escape_pct_outside_backticks(text: str) -> str:
    """Escape % as %% for Doxygen, but not inside backtick code spans."""
    parts = text.split('`')
    for i in range(0, len(parts), 2):
        parts[i] = parts[i].replace('%', '%%')
    return '`'.join(parts)


# Lazily initialized at module level on first use
_PARAM_PATTERN = None


def _get_param_pattern():
    global _PARAM_PATTERN  # noqa: PLW0603  pylint: disable=global-statement
    if _PARAM_PATTERN is None:
        _PARAM_PATTERN = _build_param_name_pattern()
    return _PARAM_PATTERN


def _build_reverse_dep_map() -> Dict[str, List[Tuple[str, str]]]:
    """Build map from target param -> [(relation, source_param), ...] from DEPENDENCIES."""
    from ..definitions import DEPENDENCIES  # pylint: disable=import-outside-toplevel
    reverse: Dict[str, List[Tuple[str, str]]] = {}
    for param, dep in DEPENDENCIES.items():
        if "when_true" in dep:
            wt = dep["when_true"]
            for r in wt.get("requires", []):
                reverse.setdefault(r, []).append(("required by", param))
            for r in wt.get("recommends", []):
                reverse.setdefault(r, []).append(("recommended for", param))
        if "when_value" in dep:
            for val, subspec in dep["when_value"].items():
                for r in subspec.get("requires", []):
                    reverse.setdefault(r, []).append(("required by", f"{param}={val}"))
                for r in subspec.get("recommends", []):
                    reverse.setdefault(r, []).append(("recommended for", f"{param}={val}"))
    return reverse


_REVERSE_DEPS = None


def _get_reverse_deps():
    global _REVERSE_DEPS  # noqa: PLW0603  pylint: disable=global-statement
    if _REVERSE_DEPS is None:
        _REVERSE_DEPS = _build_reverse_dep_map()
    return _REVERSE_DEPS


# Feature tag → human-readable label, in priority order
_TAG_LABELS = [
    ("bubbles", "Bubble model parameter"),
    ("mhd", "MHD parameter"),
    ("chemistry", "Chemistry parameter"),
    ("time", "Time-stepping parameter"),
    ("grid", "Grid parameter"),
    ("weno", "WENO parameter"),
    ("viscosity", "Viscosity parameter"),
    ("elasticity", "Elasticity parameter"),
    ("surface_tension", "Surface tension parameter"),
    ("acoustic", "Acoustic parameter"),
    ("ib", "Immersed boundary parameter"),
    ("probes", "Probe/integral parameter"),
    ("riemann", "Riemann solver parameter"),
    ("relativity", "Relativity parameter"),
    ("output", "Output parameter"),
]

# Prefix → label for untagged params
_PREFIX_LABELS = [
    ("mixlayer_", "Mixing layer parameter"),
    ("nv_uvm_", "GPU memory management"),
    ("ic_", "Initial condition parameter"),
]

# BC sub-parameter attribute → usage annotation
# Keyed by the attribute name after % (without index suffix)
_BC_ATTR_ANNOTATIONS = {
    "grcbc_in": "Enables GRCBC subsonic inflow (bc type -7)",
    "grcbc_out": "Enables GRCBC subsonic outflow (bc type -8)",
    "grcbc_vel_out": "GRCBC velocity outlet (requires `grcbc_out`)",
    "vel_in": "Inlet velocity component (used with `grcbc_in`)",
    "vel_out": "Outlet velocity component (used with `grcbc_vel_out`)",
    "pres_in": "Inlet pressure (used with `grcbc_in`)",
    "pres_out": "Outlet pressure (used with `grcbc_out`)",
    "alpha_rho_in": "Inlet partial density per fluid (used with `grcbc_in`)",
    "alpha_in": "Inlet volume fraction per fluid (used with `grcbc_in`)",
    "vb1": "Boundary velocity component 1 at domain begin",
    "vb2": "Boundary velocity component 2 at domain begin",
    "vb3": "Boundary velocity component 3 at domain begin",
    "ve1": "Boundary velocity component 1 at domain end",
    "ve2": "Boundary velocity component 2 at domain end",
    "ve3": "Boundary velocity component 3 at domain end",
}

# patch_bc attribute → usage annotation
_PATCH_BC_ANNOTATIONS = {
    "geometry": "Patch shape: 1=line, 2=circle, 3=rectangle",
    "type": "BC type applied within patch region",
    "dir": "Patch normal direction (1=x, 2=y, 3=z)",
    "loc": "Domain boundary (-1=begin, 1=end)",
    "centroid": "Patch center coordinate",
    "length": "Patch dimension",
    "radius": "Patch radius (geometry=2)",
}

# simplex_params attribute → usage annotation
_SIMPLEX_ANNOTATIONS = {
    "perturb_dens": "Enable simplex density perturbation",
    "perturb_dens_freq": "Density perturbation frequency",
    "perturb_dens_scale": "Density perturbation amplitude",
    "perturb_dens_offset": "Density perturbation offset seed",
    "perturb_vel": "Enable simplex velocity perturbation",
    "perturb_vel_freq": "Velocity perturbation frequency",
    "perturb_vel_scale": "Velocity perturbation amplitude",
    "perturb_vel_offset": "Velocity perturbation offset seed",
}

# fluid_pp attribute → usage annotation (EOS params only)
_FLUID_PP_ANNOTATIONS = {
    "gamma": "Specific heat ratio (EOS)",
    "pi_inf": "Stiffness pressure (EOS)",
    "cv": "Specific heat at constant volume",
    "qv": "Heat of formation",
    "qvp": "Heat of formation derivative",
}


def _get_attr_annotation(param_name: str) -> str:
    """Look up annotation for compound param names (e.g. bc_x%vel_in(1) → 'Inlet velocity')."""
    if '%' not in param_name:
        return ""
    # Extract attribute after last %: "bc_x%vel_in(1)" → "vel_in(1)"
    attr_full = param_name.split('%')[-1]
    # Strip index suffix: "vel_in(1)" → "vel_in", "centroid(2)" → "centroid"
    attr_base = re.match(r'^([a-zA-Z_0-9]+)', attr_full)
    if not attr_base:
        return ""
    attr_key = attr_base.group(1)
    # Determine family prefix: "bc_x%..." → "bc_", "patch_bc(1)%..." → "patch_bc"
    prefix = param_name.split('%')[0]
    if re.match(r'bc_[xyz]', prefix):
        return _BC_ATTR_ANNOTATIONS.get(attr_key, "")
    if prefix.startswith("patch_bc"):
        return _PATCH_BC_ANNOTATIONS.get(attr_key, "")
    if prefix.startswith("simplex_params"):
        return _SIMPLEX_ANNOTATIONS.get(attr_key, "")
    if prefix.startswith("fluid_pp"):
        return _FLUID_PP_ANNOTATIONS.get(attr_key, "")
    return ""


def _format_tag_annotation(param_name: str, param) -> str:  # pylint: disable=too-many-locals
    """
    Return a short annotation for params with no schema constraints and no AST rules.

    Checks (in order): own DEPENDENCIES, output flag tags, reverse dependencies,
    feature tag labels, prefix-group labels, and compound-name attribute annotations.
    """
    # 1. Own DEPENDENCIES info
    if param.dependencies:
        dep = param.dependencies
        if "when_true" in dep:
            wt = dep["when_true"]
            if "requires" in wt:
                req = ", ".join(f"`{r}`" for r in wt["requires"])
                return f"Requires {req} when enabled"
            if "requires_value" in wt:
                parts = []
                for k, vals in wt["requires_value"].items():
                    parts.append(f"`{k}` in {vals}")
                return "Requires " + ", ".join(parts)
            if "recommends" in wt:
                rec = ", ".join(f"`{r}`" for r in wt["recommends"])
                return f"Recommends {rec}"

    # 2. Tag-based output flag label (specific labels for LOG output params)
    if "output" in param.tags and param.param_type == ParamType.LOG:
        if "bubbles" in param.tags:
            return "Lagrangian output flag"
        if "chemistry" in param.tags:
            return "Chemistry output flag"
        return "Post-processing output flag"

    # 3. Reverse dependencies (params required/recommended by other features)
    reverse = _get_reverse_deps()
    if param_name in reverse:
        entries = reverse[param_name]
        parts = []
        for relation, source in entries[:2]:
            parts.append(f"Required by `{source}`" if relation == "required by"
                         else f"Recommended for `{source}`")
        return "; ".join(parts)

    # 4. Feature tag labels (broader tag-based annotation)
    for tag, label in _TAG_LABELS:
        if tag in param.tags:
            return label

    # 5. Prefix group labels (for untagged params with common prefixes)
    for prefix, label in _PREFIX_LABELS:
        if param_name.startswith(prefix):
            return label

    # 6. Compound-name attribute annotations (bc_x%vel_in, patch_bc%geometry, etc.)
    attr_ann = _get_attr_annotation(param_name)
    if attr_ann:
        return attr_ann

    return ""


def _format_validator_rules(param_name: str, by_trigger: Dict[str, list],  # pylint: disable=too-many-locals
                            by_param: Dict[str, list] | None = None) -> str:
    """Format AST-extracted validator rules for a parameter's Constraints column.

    Gets rules where this param is the trigger.  Falls back to by_param
    (rules that mention this param) when no trigger rules exist.
    """
    rules = by_trigger.get(param_name, [])
    if not rules and by_param:
        rules = by_param.get(param_name, [])
    if not rules:
        return ""

    pattern = _get_param_pattern()

    # Deduplicate messages (same message can appear from multiple loop iterations)
    seen = set()
    unique_rules = []
    for r in rules:
        if r.message not in seen:
            seen.add(r.message)
            unique_rules.append(r)

    # Classify and pick representative messages
    requirements = []
    incompatibilities = []
    ranges = []
    others = []

    for r in unique_rules:
        kind = classify_message(r.message)
        msg = _backtick_params(r.message, pattern)
        if kind == "requirement":
            requirements.append(msg)
        elif kind == "incompatibility":
            incompatibilities.append(msg)
        elif kind == "range":
            ranges.append(msg)
        else:
            others.append(msg)

    # Build concise output - show up to 3 rules total, prioritized
    parts = []
    budget = 3
    for group in [requirements, incompatibilities, ranges, others]:
        for msg in group:
            if budget <= 0:
                break
            parts.append(msg)
            budget -= 1

    return "; ".join(parts)


def generate_parameter_docs() -> str:  # pylint: disable=too-many-locals,too-many-statements
    """Generate markdown documentation for all parameters."""
    # AST-extract rules from case_validator.py
    analysis = analyze_case_validator()
    by_trigger = analysis["by_trigger"]
    by_param = analysis["by_param"]

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
            # Check if any member of a pattern has constraints
            pattern_has_constraints = False
            for _pattern, examples in patterns.items():
                for ex in examples:
                    p = REGISTRY.all_params[ex]
                    if p.constraints or ex in by_trigger or ex in by_param:
                        pattern_has_constraints = True
                        break
                if pattern_has_constraints:
                    break

            lines.append("### Patterns")
            lines.append("")
            if pattern_has_constraints:
                lines.append("| Pattern | Example | Description | Constraints |")
                lines.append("|---------|---------|-------------|-------------|")
            else:
                lines.append("| Pattern | Example | Description |")
                lines.append("|---------|---------|-------------|")

            for pattern, examples in sorted(patterns.items()):
                example = examples[0]
                desc = get_description(example) or ""
                # Truncate long descriptions
                if len(desc) > 60:
                    desc = desc[:57] + "..."
                # Escape % for Doxygen (even inside backtick code spans)
                pattern_escaped = _escape_percent(pattern)
                example_escaped = _escape_percent(example)
                desc = _escape_percent(desc)
                if pattern_has_constraints:
                    p = REGISTRY.all_params[example]
                    constraints = _format_constraints(p)
                    deps = _format_validator_rules(example, by_trigger, by_param)
                    extra = "; ".join(filter(None, [constraints, deps]))
                    if not extra:
                        extra = _format_tag_annotation(example, p)
                    extra = _escape_pct_outside_backticks(extra)
                    lines.append(f"| `{pattern_escaped}` | `{example_escaped}` | {desc} | {extra} |")
                else:
                    lines.append(f"| `{pattern_escaped}` | `{example_escaped}` | {desc} |")

            lines.append("")
        else:
            # Full table - no patterns to collapse
            lines.append("| Parameter | Type | Description | Constraints |")
            lines.append("|-----------|------|-------------|-------------|")

            for name, param in params:
                type_str = _type_to_str(param.param_type)
                desc = get_description(name) or ""
                # Truncate long descriptions
                if len(desc) > 80:
                    desc = desc[:77] + "..."
                constraints = _format_constraints(param)
                deps = _format_validator_rules(name, by_trigger, by_param)
                extra = "; ".join(filter(None, [constraints, deps]))
                if not extra:
                    extra = _format_tag_annotation(name, param)
                extra = _escape_pct_outside_backticks(extra)
                # Escape % for Doxygen (even inside backtick code spans)
                name_escaped = _escape_percent(name)
                desc = _escape_percent(desc)
                lines.append(f"| `{name_escaped}` | {type_str} | {desc} | {extra} |")

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
