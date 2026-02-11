"""
Parameter Validation with Constraints and Dependencies.

Provides enhanced validation beyond JSON schema type checking.

Relationship to case_validator.py
---------------------------------
This module (params/validate.py) provides **generic parameter validation**:
- Type checking (int, real, string, etc.)
- Range constraints (min/max values)
- Choice validation (enum-like constraints)
- Parameter dependencies (requires/recommends)

The case_validator.py module provides **domain-specific physics validation**:
- Cross-parameter consistency checks (e.g., bubble model + polytropic settings)
- Model-specific requirements (e.g., WENO order constraints)
- Stage-specific validation (pre_process, simulation, post_process)
- 50+ physics-aware constraint checks

These modules are complementary:
- params/validate.py: Fast, generic checks that apply to all parameters
- case_validator.py: Comprehensive physics validation for MFC simulations

Typical usage:
    1. JSON schema validation (via mfc-case-schema.json)
    2. Generic constraint validation (via this module)
    3. Physics validation (via case_validator.py)
"""

from typing import Dict, Any, List, Optional, Tuple
from .registry import REGISTRY
from .errors import (
    dependency_error,
    dependency_recommendation,
    dependency_value_error,
    format_error_list,
    unknown_param_error,
)
from .suggest import suggest_parameter
# Note: definitions is imported by params/__init__.py to populate REGISTRY.
# This redundant import ensures REGISTRY is populated even if this module
# is imported directly (e.g., during testing).
from . import definitions  # noqa: F401  pylint: disable=unused-import


def check_unknown_params(params: Dict[str, Any]) -> List[str]:
    """
    Check for unknown parameters and suggest corrections.

    Uses fuzzy matching via rapidfuzz to provide "Did you mean?" suggestions
    for parameter names that don't exist in the registry.

    Args:
        params: Dictionary of parameter name -> value

    Returns:
        List of error messages for unknown parameters with suggestions.
    """
    errors = []

    for name in params.keys():
        if name not in REGISTRY.all_params:
            suggestions = suggest_parameter(name)
            errors.append(unknown_param_error(name, suggestions))

    return errors


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
            continue  # Unknown params handled by check_unknown_params

        # Skip analytic expressions (strings for numeric params)
        if isinstance(value, str) and param_def.param_type.value.startswith("analytic"):
            continue
        if isinstance(value, str) and param_def.param_type.value in ("int", "real"):
            # This is an analytic expression, skip constraint check
            continue

        param_errors = param_def.validate_value(value)
        errors.extend(param_errors)

    return errors


def _check_condition(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    name: str,
    condition: Dict[str, Any],
    condition_label: Optional[str],
    params: Dict[str, Any],
    errors: List[str],
    warnings: List[str],
) -> None:
    """Check a single condition dict (requires, recommends, requires_value)."""
    if "requires" in condition:
        for req in condition["requires"]:
            if req not in params:
                errors.append(dependency_error(name, req, condition_label))

    if "recommends" in condition:
        for rec in condition["recommends"]:
            if rec not in params:
                warnings.append(dependency_recommendation(name, rec, condition_label))

    if "requires_value" in condition:
        for req_param, expected_vals in condition["requires_value"].items():
            if req_param not in params:
                errors.append(dependency_error(name, req_param, condition_label))
            else:
                got = params[req_param]
                if got not in expected_vals:
                    errors.append(dependency_value_error(
                        name, condition_label, req_param, expected_vals, got,
                    ))


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
            _check_condition(name, deps["when_true"], "=T", params, errors, warnings)

        # Check "when_set" dependencies (for any param that's set)
        if "when_set" in deps:
            _check_condition(name, deps["when_set"], None, params, errors, warnings)

        # Check "when_value" dependencies (value-specific conditions)
        if "when_value" in deps:
            for trigger_val, condition in deps["when_value"].items():
                if value == trigger_val:
                    _check_condition(
                        name, condition, f"={trigger_val}", params, errors, warnings,
                    )

    return errors, warnings


def validate_case(
    params: Dict[str, Any],
    warn: bool = True,
    check_unknown: bool = True,
) -> Tuple[List[str], List[str]]:
    """
    Full validation of case parameters.

    Args:
        params: Dictionary of parameter name -> value
        warn: Whether to check for warnings (recommended params)
        check_unknown: Whether to check for unknown parameters

    Returns:
        Tuple of (errors, warnings)
    """
    errors = []
    warnings = []

    # Check for unknown parameters with "did you mean?" suggestions
    if check_unknown:
        unknown_errors = check_unknown_params(params)
        errors.extend(unknown_errors)

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
    """
    Format validation results for display.

    Uses the centralized formatting from errors.py for consistency
    across all validation systems.
    """
    return format_error_list(errors, warnings, use_rich=True)
