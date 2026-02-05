"""
MFC Parameter Search and Discovery Command.

Provides CLI access to search and explore MFC's ~3,300 case parameters.
"""
# pylint: disable=import-outside-toplevel

import re
from .state import ARG
from .printer import cons


def params():
    """Execute the params command based on CLI arguments."""
    from .params import REGISTRY
    from .params import definitions  # noqa: F401  pylint: disable=unused-import

    query = ARG("query")
    type_filter = ARG("param_type")
    show_families = ARG("families")
    show_features = ARG("features")
    feature_name = ARG("feature")
    names_only = ARG("names_only")
    show_count = ARG("count")
    limit = ARG("limit")
    describe = ARG("describe")

    # By default, search both names and descriptions (more user-friendly)
    search_descriptions = not names_only

    if show_count:
        _show_statistics(REGISTRY)
    elif show_features:
        _show_feature_groups(REGISTRY)
    elif feature_name:
        _show_feature_params(REGISTRY, feature_name, type_filter, limit, describe)
    elif show_families:
        _show_families(REGISTRY, limit)
    elif query:
        _search_params(REGISTRY, query, type_filter, limit, describe, search_descriptions)
    else:
        _show_statistics(REGISTRY)
        cons.print()
        cons.print("[yellow]Tip:[/yellow] Use './mfc.sh params <query>' to search for parameters")
        cons.print("     Use './mfc.sh params --feature mhd' to see MHD parameters")
        cons.print("     Use './mfc.sh params -F' to see all feature groups")


def _collapse_indexed_params(matches):  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    """
    Collapse indexed parameters into patterns.

    Handles multiple index patterns:
    - Suffix index: bc_z%alpha_in(1) -> bc_z%alpha_in(N)
    - Prefix index: patch_icpp(1)%geometry -> patch_icpp(N)%geometry
    - Both: patch_icpp(1)%alpha(1) -> patch_icpp(N)%alpha(M)
    """
    # Patterns for different index positions
    # Pattern 1: prefix(N)%suffix or prefix(N)%suffix(M)
    prefix_pattern = re.compile(r'^([^(]+)\((\d+)\)%(.+)$')
    # Pattern 2: name(N) or name(N, M) at end
    suffix_pattern = re.compile(r'^(.+)\((\d+)(?:,\s*(\d+))?\)$')

    # Two-level grouping: first by base pattern (with indices replaced), then collect indices
    groups = {}  # normalized_pattern -> {indices: [...], param_type, stages, pattern_type}

    for name, param in matches:
        # Try prefix pattern first: patch_icpp(1)%geometry
        prefix_match = prefix_pattern.match(name)
        if prefix_match:
            prefix = prefix_match.group(1)
            idx1 = int(prefix_match.group(2))
            suffix = prefix_match.group(3)

            # Check if suffix also has an index
            suffix_match = suffix_pattern.match(suffix)
            if suffix_match:
                suffix_base = suffix_match.group(1)
                idx2 = int(suffix_match.group(2))
                idx3 = int(suffix_match.group(3)) if suffix_match.group(3) else None
                # Pattern: prefix(N)%suffix_base(M) or prefix(N)%suffix_base(M, K)
                if idx3 is not None:
                    base_pattern = f"{prefix}(N)%{suffix_base}(M, K)"
                    indices_key = (idx1, idx2, idx3)
                else:
                    base_pattern = f"{prefix}(N)%{suffix_base}(M)"
                    indices_key = (idx1, idx2, None)
            else:
                # Pattern: prefix(N)%suffix
                base_pattern = f"{prefix}(N)%{suffix}"
                indices_key = (idx1, None, None)

            if base_pattern not in groups:
                groups[base_pattern] = {
                    'indices': [],
                    'param_type': param.param_type,
                                    }
            groups[base_pattern]['indices'].append((indices_key, param))
            continue

        # Try suffix-only pattern: name(N) or name(N, M)
        suffix_match = suffix_pattern.match(name)
        if suffix_match:
            base = suffix_match.group(1)
            idx1 = int(suffix_match.group(2))
            idx2 = int(suffix_match.group(3)) if suffix_match.group(3) else None

            if idx2 is not None:
                base_pattern = f"{base}(N, M)"
                indices_key = (idx1, idx2, None)
            else:
                base_pattern = f"{base}(N)"
                indices_key = (idx1, None, None)

            if base_pattern not in groups:
                groups[base_pattern] = {
                    'indices': [],
                    'param_type': param.param_type,
                                    }
            groups[base_pattern]['indices'].append((indices_key, param))
            continue

        # No index pattern - add as-is
        if name not in groups:
            groups[name] = {
                'indices': [(None, param)],
                'param_type': param.param_type,
                            }
        else:
            groups[name]['indices'].append((None, param))

    # Build collapsed results
    collapsed = []

    for pattern, data in sorted(groups.items()):
        indices = data['indices']
        param = indices[0][1]  # Get param from first entry
        count = len(indices)

        if count == 1 and indices[0][0] is None:
            # Non-indexed parameter
            collapsed.append((pattern, param, 1))
        elif count == 1:
            # Single indexed parameter - show with actual index
            idx_tuple = indices[0][0]
            # Reconstruct the actual name
            actual_name = pattern
            if idx_tuple[0] is not None:
                actual_name = actual_name.replace('(N)', f'({idx_tuple[0]})', 1)
            if idx_tuple[1] is not None:
                actual_name = actual_name.replace('(M)', f'({idx_tuple[1]})', 1)
            if idx_tuple[2] is not None:
                actual_name = actual_name.replace('(K)', f'({idx_tuple[2]})', 1)
            collapsed.append((actual_name, param, 1))
        else:
            # Multiple indices - build range string
            range_parts = []

            # Extract index values
            idx1_vals = sorted(set(idx[0] for idx, _ in indices if idx and idx[0] is not None))
            idx2_vals = sorted(set(idx[1] for idx, _ in indices if idx and idx[1] is not None))
            idx3_vals = sorted(set(idx[2] for idx, _ in indices if idx and idx[2] is not None))

            if idx1_vals:
                if idx1_vals == list(range(min(idx1_vals), max(idx1_vals) + 1)):
                    range_parts.append(f"N={min(idx1_vals)}..{max(idx1_vals)}")
                else:
                    range_parts.append(f"N={min(idx1_vals)}..{max(idx1_vals)}")
            if idx2_vals:
                if idx2_vals == list(range(min(idx2_vals), max(idx2_vals) + 1)):
                    range_parts.append(f"M={min(idx2_vals)}..{max(idx2_vals)}")
                else:
                    range_parts.append(f"M={min(idx2_vals)}..{max(idx2_vals)}")
            if idx3_vals:
                range_parts.append(f"K={min(idx3_vals)}..{max(idx3_vals)}")

            range_str = ", ".join(range_parts) if range_parts else ""
            collapsed.append((pattern, param, count, range_str))

    # Sort by name
    collapsed.sort(key=lambda x: x[0])

    return collapsed


def _show_statistics(registry):
    """Show parameter count statistics."""
    cons.print("[bold]MFC Parameter Statistics[/bold]")
    cons.print()

    cons.print(f"  Total parameters: [cyan]{len(registry.all_params)}[/cyan]")

    # Count by type
    by_type = {}
    for param in registry.all_params.values():
        tname = param.param_type.name
        by_type[tname] = by_type.get(tname, 0) + 1

    cons.print()
    cons.print("  By type:")
    for tname, count in sorted(by_type.items(), key=lambda x: -x[1]):
        cons.print(f"    {tname:15} {count:5}")


def _show_feature_groups(registry):
    """Show available feature groups."""
    from .params.descriptions import FEATURE_DESCRIPTIONS

    cons.print("[bold]Feature Groups[/bold]")
    cons.print()
    cons.print("  Use './mfc.sh params --feature <name>' to see parameters for a feature.")
    cons.print()
    cons.print(f"  {'Feature':<20} {'Description'}")
    cons.print(f"  {'-'*20} {'-'*50}")

    # Get all tags from registry and show with descriptions
    all_tags = registry.get_all_tags()
    for tag in sorted(all_tags):
        desc = FEATURE_DESCRIPTIONS.get(tag, "")
        cons.print(f"  [cyan]{tag:<20}[/cyan] {desc}")

    cons.print()
    cons.print("[yellow]Example:[/yellow] ./mfc.sh params --feature mhd")


def _show_feature_params(registry, feature_name, type_filter, limit, describe):
    """Show all parameters for a feature group."""
    from .params.descriptions import FEATURE_DESCRIPTIONS

    # Check if feature exists in registry
    all_tags = registry.get_all_tags()
    if feature_name not in all_tags:
        cons.print(f"[red]Unknown feature group: '{feature_name}'[/red]")
        cons.print()
        cons.print("Available feature groups:")
        for name in sorted(all_tags):
            cons.print(f"  {name}")
        return

    # Get params by tag from registry (single source of truth)
    tagged_params = registry.get_params_by_tag(feature_name)

    # Build matches list
    matches = []
    for name, param in tagged_params.items():
        # Apply type filter
        if type_filter and param.param_type.value != type_filter:
            continue
        matches.append((name, param))

    if not matches:
        cons.print(f"[yellow]No parameters found for feature '{feature_name}'[/yellow]")
        return

    # Collapse indexed parameters
    collapsed = _collapse_indexed_params(matches)

    desc = FEATURE_DESCRIPTIONS.get(feature_name, feature_name.title() + " parameters")
    cons.print(f"[bold]{desc}[/bold] ({len(matches)} params, {len(collapsed)} unique patterns)")
    cons.print()

    # Show collapsed results
    _show_collapsed_results(collapsed[:limit], describe)

    if len(collapsed) > limit:
        cons.print()
        cons.print(f"  [dim]... {len(collapsed) - limit} more patterns (use -n {len(collapsed)} to show all)[/dim]")


def _show_families(registry, limit):
    """Show parameter families grouped by prefix."""
    families = {}
    for name in registry.all_params.keys():
        if "%" in name:
            # Get prefix before %, then strip any index: patch_icpp(1)%x -> patch_icpp
            prefix = name.split("%")[0]
            if "(" in prefix:
                prefix = prefix.split("(")[0]
        elif "(" in name:
            # Simple indexed: chem_wrt_Y(0) -> chem_wrt_Y
            prefix = name.split("(")[0]
        else:
            continue  # Skip simple params
        families[prefix] = families.get(prefix, 0) + 1

    sorted_families = sorted(families.items(), key=lambda x: -x[1])

    cons.print("[bold]Parameter Families[/bold]")
    cons.print()
    cons.print(f"  {'Family':<40} {'Count':>6}")
    cons.print(f"  {'-'*40} {'-'*6}")

    for prefix, count in sorted_families[:limit]:
        cons.print(f"  {prefix:<40} {count:>6}")

    if len(sorted_families) > limit:
        cons.print(f"  ... and {len(sorted_families) - limit} more families")

    cons.print()
    cons.print("[yellow]Tip:[/yellow] Use './mfc.sh params <family>' to see parameters in a family")


def _search_params(registry, query, type_filter, limit, describe=False, search_descriptions=True):  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    """Search for parameters matching a query."""
    from .params.descriptions import get_description

    query_lower = query.lower()
    matches = []
    desc_matches = set()  # Track which params matched via description

    for name, param in registry.all_params.items():
        name_match = query_lower in name.lower()
        desc_match = False

        if search_descriptions and not name_match:
            # Also search in description
            desc = get_description(name)
            if desc and query_lower in desc.lower():
                desc_match = True
                desc_matches.add(name)

        if not name_match and not desc_match:
            continue

        # Apply type filter
        if type_filter and param.param_type.value != type_filter:
            continue

        matches.append((name, param))

    if not matches:
        cons.print(f"[yellow]No parameters found matching '{query}'[/yellow]")
        _suggest_alternatives(registry, query)
        return

    # Collapse indexed parameters
    collapsed = _collapse_indexed_params(matches)

    cons.print(f"[bold]Parameters matching '{query}'[/bold] ({len(matches)} params, {len(collapsed)} unique patterns)")
    cons.print()

    # Show collapsed results (enable describe mode if we matched via description)
    show_describe = describe or (search_descriptions and len(desc_matches) > 0)
    _show_collapsed_results(collapsed[:limit], show_describe)

    if len(collapsed) > limit:
        cons.print()
        cons.print(f"  [dim]... {len(collapsed) - limit} more patterns (use -n {len(collapsed)} to show all)[/dim]")


def _show_collapsed_results(collapsed, describe=False):  # pylint: disable=too-many-branches
    """Show collapsed search results."""
    from .params.descriptions import get_description, get_pattern_description

    # Check if any items have index ranges to show
    has_ranges = any(len(item) == 4 and item[2] > 1 for item in collapsed)

    if describe:
        # Description mode: one param per block with description
        for item in collapsed:
            name = item[0]
            param = item[1]
            count = item[2]
            range_str = item[3] if len(item) == 4 else ""

            # Get description - use pattern description for indexed params
            if "(N)" in name or "(M)" in name:
                desc = get_pattern_description(name)
            else:
                desc = get_description(name)

            cons.print(f"  [cyan]{name}[/cyan]")
            cons.print(f"    Type: {param.param_type.name}")
            if count > 1:
                cons.print(f"    Count: {count}  ({range_str})")
            if desc:
                cons.print(f"    [dim]{desc}[/dim]")
            cons.print()
    else:
        # Compact table mode
        if has_ranges:
            cons.print(f"  {'Parameter':<40} {'Type':12} {'#':>4}  {'Index Range'}")
            cons.print(f"  {'-'*40} {'-'*12} {'-'*4}  {'-'*15}")
        else:
            cons.print(f"  {'Parameter':<40} {'Type':12}")
            cons.print(f"  {'-'*40} {'-'*12}")

        for item in collapsed:
            if len(item) == 4:
                name, param, count, range_str = item
                if count > 1:
                    cons.print(f"  {name:<40} {param.param_type.name:12} {count:>4}  {range_str}")
                else:
                    if has_ranges:
                        cons.print(f"  {name:<40} {param.param_type.name:12} {count:>4}")
                    else:
                        cons.print(f"  {name:<40} {param.param_type.name:12}")
            else:
                name, param, count = item
                if has_ranges:
                    cons.print(f"  {name:<40} {param.param_type.name:12} {count:>4}")
                else:
                    cons.print(f"  {name:<40} {param.param_type.name:12}")


def _suggest_alternatives(registry, query):
    """Suggest similar parameter names."""
    import difflib
    all_names = list(registry.all_params.keys())
    suggestions = difflib.get_close_matches(query, all_names, n=5, cutoff=0.5)

    if suggestions:
        cons.print()
        cons.print("[yellow]Did you mean:[/yellow]")
        for s in suggestions:
            cons.print(f"  {s}")
