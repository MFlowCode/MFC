"""
MFC Parameter Search and Discovery Command.

Provides CLI access to search and explore MFC's ~3,300 case parameters.
"""

import re
from .state import ARG
from .printer import cons


def params():
    """Execute the params command based on CLI arguments."""
    from .params import REGISTRY
    from .params import definitions  # noqa: F401 - loads definitions
    from .params.schema import Stage

    query = ARG("query")
    stage_filter = ARG("stage")
    type_filter = ARG("param_type")
    show_families = ARG("families")
    show_count = ARG("count")
    limit = ARG("limit")

    # Map stage names to Stage enum
    stage_map = {
        "common": Stage.COMMON,
        "pre_process": Stage.PRE_PROCESS,
        "simulation": Stage.SIMULATION,
        "post_process": Stage.POST_PROCESS,
    }

    if show_count:
        _show_statistics(REGISTRY)
    elif show_families:
        _show_families(REGISTRY, limit)
    elif query:
        _search_params(REGISTRY, query, stage_map.get(stage_filter), type_filter, limit)
    elif stage_filter:
        _list_stage_params(REGISTRY, stage_map[stage_filter], type_filter, limit)
    else:
        _show_statistics(REGISTRY)
        cons.print()
        cons.print("[yellow]Tip:[/yellow] Use './mfc.sh params <query>' to search for parameters")
        cons.print("     Use './mfc.sh params -f' to see parameter families")


def _collapse_indexed_params(matches):
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
                    'stages': param.stages,
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
                    'stages': param.stages,
                }
            groups[base_pattern]['indices'].append((indices_key, param))
            continue

        # No index pattern - add as-is
        if name not in groups:
            groups[name] = {
                'indices': [(None, param)],
                'param_type': param.param_type,
                'stages': param.stages,
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
    from .params.schema import Stage

    cons.print("[bold]MFC Parameter Statistics[/bold]")
    cons.print()

    # Count by stage
    by_stage = {}
    for name, param in registry.all_params.items():
        for stage in param.stages:
            by_stage[stage.name] = by_stage.get(stage.name, 0) + 1

    cons.print(f"  Total parameters: [cyan]{len(registry.all_params)}[/cyan]")
    cons.print()
    cons.print("  By stage:")
    for stage in Stage:
        count = by_stage.get(stage.name, 0)
        cons.print(f"    {stage.name:15} {count:5}")

    # Count by type
    by_type = {}
    for param in registry.all_params.values():
        tname = param.param_type.name
        by_type[tname] = by_type.get(tname, 0) + 1

    cons.print()
    cons.print("  By type:")
    for tname, count in sorted(by_type.items(), key=lambda x: -x[1]):
        cons.print(f"    {tname:15} {count:5}")


def _show_families(registry, limit):
    """Show parameter families grouped by prefix."""
    families = {}
    for name in registry.all_params.keys():
        if "%" in name:
            prefix = name.split("%")[0]
        elif "(" in name:
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


def _search_params(registry, query, stage_filter, type_filter, limit):
    """Search for parameters matching a query."""
    query_lower = query.lower()
    matches = []

    for name, param in registry.all_params.items():
        if query_lower not in name.lower():
            continue

        # Apply stage filter
        if stage_filter and stage_filter not in param.stages:
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

    # Show collapsed results
    _show_collapsed_results(collapsed[:limit])

    if len(collapsed) > limit:
        cons.print()
        cons.print(f"  [dim]... {len(collapsed) - limit} more patterns (use -n {len(collapsed)} to show all)[/dim]")


def _list_stage_params(registry, stage, type_filter, limit):
    """List all parameters for a specific stage."""
    params = registry.get_params_by_stage(stage)

    # Apply type filter
    if type_filter:
        params = {k: v for k, v in params.items() if v.param_type.value == type_filter}

    cons.print(f"[bold]{stage.name} Parameters[/bold] ({len(params)} total)")
    cons.print()

    # Group by prefix for summary
    by_prefix = {}
    simple = []
    for name, param in sorted(params.items()):
        if "%" in name or "(" in name:
            if "%" in name:
                prefix = name.split("%")[0]
            else:
                prefix = name.split("(")[0]
            if prefix not in by_prefix:
                by_prefix[prefix] = []
            by_prefix[prefix].append((name, param))
        else:
            simple.append((name, param))

    # Show simple params first
    if simple:
        cons.print("  [cyan]Simple parameters:[/cyan]")
        for name, param in simple[:limit]:
            cons.print(f"    {name:<30} {param.param_type.name:6}")
        if len(simple) > limit:
            cons.print(f"    [dim]... and {len(simple) - limit} more[/dim]")
        cons.print()

    # Show families summary
    if by_prefix:
        cons.print("  [cyan]Parameter families:[/cyan]")
        sorted_families = sorted(by_prefix.items(), key=lambda x: -len(x[1]))
        for prefix, items in sorted_families[:limit]:
            cons.print(f"    {prefix}: {len(items)} parameters")
        if len(sorted_families) > limit:
            cons.print(f"    [dim]... and {len(sorted_families) - limit} more families[/dim]")

    cons.print()
    cons.print("[yellow]Tip:[/yellow] Use './mfc.sh params <family>' to search within a family")


def _show_collapsed_results(collapsed):
    """Show collapsed search results."""
    cons.print(f"  {'Parameter':<40} {'Type':6} {'Count':>6}  {'Range/Stages'}")
    cons.print(f"  {'-'*40} {'-'*6} {'-'*6}  {'-'*20}")

    for item in collapsed:
        if len(item) == 4:
            name, param, count, range_str = item
            stages = ",".join(s.name[:3] for s in param.stages)
            if count > 1:
                cons.print(f"  {name:<40} {param.param_type.name:6} {count:>6}  {range_str}")
            else:
                cons.print(f"  {name:<40} {param.param_type.name:6} {count:>6}  {stages}")
        else:
            name, param, count = item
            stages = ",".join(s.name[:3] for s in param.stages)
            cons.print(f"  {name:<40} {param.param_type.name:6} {count:>6}  {stages}")


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
