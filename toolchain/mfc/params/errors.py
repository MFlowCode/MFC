"""
Consistent Error Message Formatting for MFC Parameter Validation.

Provides utility functions for creating consistent, user-friendly error
messages across all validation systems (params/validate.py, case_validator.py,
and JSON schema validation).

Error Message Format
--------------------
All error messages follow this structure:
- Parameter name in single quotes: 'param_name'
- Clear description of the problem
- Current value if relevant: got <value>
- Expected value/range if relevant: expected <constraint>

Examples:
- "'weno_order' must be one of [1, 3, 5, 7], got 4"
- "'bubbles_euler'=T requires 'nb' to be set"
- "'m' must be >= 0, got -1"
"""

from typing import Any, List, Optional


def format_param(name: str) -> str:
    """Format a parameter name for error messages."""
    return f"'{name}'"


def format_value(value: Any) -> str:
    """Format a value for error messages."""
    if isinstance(value, str):
        return f"'{value}'"
    return str(value)


def constraint_error(
    param: str,
    constraint_type: str,
    expected: Any,
    got: Any,
) -> str:
    """
    Create a constraint violation error message.

    Args:
        param: Parameter name
        constraint_type: Type of constraint ('choices', 'min', 'max')
        expected: Expected constraint value
        got: Actual value received

    Returns:
        Formatted error message.
    """
    if constraint_type == "choices":
        return f"{format_param(param)} must be one of {expected}, got {format_value(got)}"
    if constraint_type == "min":
        return f"{format_param(param)} must be >= {expected}, got {format_value(got)}"
    if constraint_type == "max":
        return f"{format_param(param)} must be <= {expected}, got {format_value(got)}"
    return f"{format_param(param)} constraint '{constraint_type}' violated: expected {expected}, got {format_value(got)}"


def type_error(param: str, expected_type: str, got: Any) -> str:
    """
    Create a type mismatch error message.

    Args:
        param: Parameter name
        expected_type: Expected type description
        got: Actual value received

    Returns:
        Formatted error message.
    """
    return f"{format_param(param)} must be {expected_type}, got {format_value(got)}"


def dependency_error(
    param: str,
    required_param: str,
    condition: Optional[str] = None,
) -> str:
    """
    Create a missing dependency error message.

    Args:
        param: Parameter that has the dependency
        required_param: Parameter that is required
        condition: Optional condition (e.g., "=T", "> 0")

    Returns:
        Formatted error message.
    """
    if condition:
        return f"{format_param(param)}{condition} requires {format_param(required_param)} to be set"
    return f"{format_param(param)} requires {format_param(required_param)} to be set"


def dependency_recommendation(
    param: str,
    recommended_param: str,
    condition: Optional[str] = None,
) -> str:
    """
    Create a dependency recommendation message.

    Args:
        param: Parameter that has the recommendation
        recommended_param: Parameter that is recommended
        condition: Optional condition

    Returns:
        Formatted recommendation message.
    """
    if condition:
        return f"{format_param(param)}{condition}: consider setting {format_param(recommended_param)}"
    return f"When {format_param(param)} is set, consider also setting {format_param(recommended_param)}"


def required_error(param: str, context: Optional[str] = None) -> str:
    """
    Create a missing required parameter error message.

    Args:
        param: Required parameter name
        context: Optional context (e.g., "when m > 0")

    Returns:
        Formatted error message.
    """
    if context:
        return f"{format_param(param)} must be set {context}"
    return f"{format_param(param)} must be set"


def mutual_exclusion_error(params: List[str], active: List[str]) -> str:
    """
    Create a mutual exclusion error message.

    Args:
        params: List of mutually exclusive parameters
        active: List of parameters that are currently active (conflicting)

    Returns:
        Formatted error message.
    """
    formatted = [format_param(p) for p in params]
    active_formatted = [format_param(p) for p in active]
    return (
        f"Only one of {', '.join(formatted)} can be enabled, "
        f"but {', '.join(active_formatted)} are all enabled"
    )


def dimension_error(param: str, requirement: str) -> str:
    """
    Create a dimensionality constraint error.

    Args:
        param: Parameter name
        requirement: Description of the dimensional requirement

    Returns:
        Formatted error message.
    """
    return f"{format_param(param)}: {requirement}"


def unknown_param_error(param: str, suggestions: Optional[List[str]] = None) -> str:
    """
    Create an error message for an unknown parameter with suggestions.

    Args:
        param: The unknown parameter name.
        suggestions: Optional list of similar valid parameter names.

    Returns:
        Formatted error message with "Did you mean?" if suggestions available.
    """
    base_msg = f"Unknown parameter {format_param(param)}"
    if suggestions:
        if len(suggestions) == 1:
            return f"{base_msg}. Did you mean {format_param(suggestions[0])}?"
        quoted = [format_param(s) for s in suggestions]
        return f"{base_msg}. Did you mean one of: {', '.join(quoted)}?"
    return base_msg


def format_error_list(
    errors: List[str],
    warnings: Optional[List[str]] = None,
    use_rich: bool = True,
) -> str:
    """
    Format a list of errors and warnings for display.

    Args:
        errors: List of error messages
        warnings: Optional list of warning messages
        use_rich: Whether to use Rich markup for colors

    Returns:
        Formatted string with errors and warnings.
    """
    lines = []

    if errors:
        if use_rich:
            lines.append("[red]Validation Errors:[/red]")
            for err in errors:
                lines.append(f"  [red]✗[/red] {err}")
        else:
            lines.append("Validation Errors:")
            for err in errors:
                lines.append(f"  ✗ {err}")

    if warnings:
        if lines:
            lines.append("")
        if use_rich:
            lines.append("[yellow]Warnings:[/yellow]")
            for warn in warnings:
                lines.append(f"  [yellow]![/yellow] {warn}")
        else:
            lines.append("Warnings:")
            for warn in warnings:
                lines.append(f"  ! {warn}")

    return "\n".join(lines)
