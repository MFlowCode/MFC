"""
Parameter Validation with Constraints and Dependencies.

Provides enhanced validation beyond JSON schema type checking.
"""

from typing import Dict, Any, List, Tuple
from .registry import REGISTRY
from . import definitions  # noqa: F401  pylint: disable=unused-import


def validate_constraints(params: Dict[str, Any]) -> List[str]:
    """
    Validate parameter values against their constraints.

    Args:
        params: Dictionary of parameter name -> value

    Returns:
        List of error messages (empty if all valid)
    """
    errors = []

    for name, value in params.items():
        param_def = REGISTRY.all_params.get(name)
        if param_def is None:
            continue  # Unknown param handled elsewhere

        # Skip analytic expressions (strings for numeric params)
        if isinstance(value, str) and param_def.param_type.value.startswith("analytic"):
            continue
        if isinstance(value, str) and param_def.param_type.value in ("int", "real"):
            # This is an analytic expression, skip constraint check
            continue

        param_errors = param_def.validate_value(value)
        errors.extend(param_errors)

    return errors


def check_dependencies(params: Dict[str, Any]) -> Tuple[List[str], List[str]]:  # pylint: disable=too-many-branches
    """
    Check parameter dependencies.

    Args:
        params: Dictionary of parameter name -> value

    Returns:
        Tuple of (errors, warnings)
        - errors: Missing required params
        - warnings: Missing recommended params
    """
    errors = []
    warnings = []

    for name, value in params.items():
        param_def = REGISTRY.all_params.get(name)
        if param_def is None or param_def.dependencies is None:
            continue

        deps = param_def.dependencies

        # Check "when_true" dependencies (for LOG params set to "T")
        if "when_true" in deps and value == "T":
            when_true = deps["when_true"]

            # Required params
            if "requires" in when_true:
                for req in when_true["requires"]:
                    if req not in params:
                        errors.append(
                            f"{name}=T requires '{req}' to be set"
                        )

            # Recommended params
            if "recommends" in when_true:
                for rec in when_true["recommends"]:
                    if rec not in params:
                        warnings.append(
                            f"{name}=T: consider setting '{rec}'"
                        )

        # Check "when_set" dependencies (for any param that's set)
        if "when_set" in deps:
            when_set = deps["when_set"]

            if "requires" in when_set:
                for req in when_set["requires"]:
                    if req not in params:
                        errors.append(
                            f"'{name}' requires '{req}' to be set"
                        )

            if "recommends" in when_set:
                for rec in when_set["recommends"]:
                    if rec not in params:
                        warnings.append(
                            f"'{name}' is set: consider also setting '{rec}'"
                        )

    return errors, warnings


def validate_case(params: Dict[str, Any], warn: bool = True) -> Tuple[List[str], List[str]]:
    """
    Full validation of case parameters.

    Args:
        params: Dictionary of parameter name -> value
        warn: Whether to check for warnings (recommended params)

    Returns:
        Tuple of (errors, warnings)
    """
    errors = []
    warnings = []

    # Check constraints
    constraint_errors = validate_constraints(params)
    errors.extend(constraint_errors)

    # Check dependencies
    if warn:
        dep_errors, dep_warnings = check_dependencies(params)
        errors.extend(dep_errors)
        warnings.extend(dep_warnings)

    return errors, warnings


def format_validation_results(errors: List[str], warnings: List[str]) -> str:
    """Format validation results for display."""
    lines = []

    if errors:
        lines.append("[red]Validation Errors:[/red]")
        for err in errors:
            lines.append(f"  [red]✗[/red] {err}")

    if warnings:
        if lines:
            lines.append("")
        lines.append("[yellow]Warnings:[/yellow]")
        for warn in warnings:
            lines.append(f"  [yellow]![/yellow] {warn}")

    return "\n".join(lines)
