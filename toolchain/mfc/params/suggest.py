"""
Fuzzy Matching for Parameter Suggestions.

Provides "did you mean?" functionality for typo detection in parameter names,
constraint keys, and dependency keys using rapidfuzz for fast string matching.

This replaces manual typo detection with automatic fuzzy matching against
the known valid options (parameter names from REGISTRY, valid constraint keys, etc.)
"""

from typing import List, Iterable
from functools import lru_cache

# Import rapidfuzz - falls back gracefully if not installed
try:
    from rapidfuzz import process, fuzz
    RAPIDFUZZ_AVAILABLE = True
except ImportError:
    RAPIDFUZZ_AVAILABLE = False

# Minimum similarity score (0-100) to consider a match
MIN_SIMILARITY_SCORE = 60

# Maximum number of suggestions to return
MAX_SUGGESTIONS = 3


def suggest_similar(
    unknown: str,
    valid_options: Iterable[str],
    min_score: int = MIN_SIMILARITY_SCORE,
    max_suggestions: int = MAX_SUGGESTIONS,
) -> List[str]:
    """
    Find similar strings from valid_options that match the unknown string.

    Uses rapidfuzz for fast fuzzy string matching. Falls back to empty list
    if rapidfuzz is not available.

    Args:
        unknown: The unknown/misspelled string to match.
        valid_options: Iterable of valid strings to match against.
        min_score: Minimum similarity score (0-100) to include a match.
        max_suggestions: Maximum number of suggestions to return.

    Returns:
        List of similar valid options, sorted by similarity (best first).
        Empty list if no good matches found or rapidfuzz not available.
    """
    if not RAPIDFUZZ_AVAILABLE:
        return []

    if not unknown or not valid_options:
        return []

    # Convert to list if needed (rapidfuzz needs indexable sequence)
    options_list = list(valid_options)
    if not options_list:
        return []

    # Use rapidfuzz to find best matches
    # process.extract returns list of (match, score, index) tuples
    matches = process.extract(
        unknown,
        options_list,
        scorer=fuzz.WRatio,  # Weighted ratio handles partial matches well
        limit=max_suggestions,
        score_cutoff=min_score,
    )

    return [match[0] for match in matches]


def format_suggestion(suggestions: List[str]) -> str:
    """
    Format a "did you mean?" suggestion message.

    Args:
        suggestions: List of suggested alternatives.

    Returns:
        Formatted suggestion string, or empty string if no suggestions.
    """
    if not suggestions:
        return ""

    if len(suggestions) == 1:
        return f"Did you mean '{suggestions[0]}'?"
    quoted = [f"'{s}'" for s in suggestions]
    return f"Did you mean one of: {', '.join(quoted)}?"


def suggest_parameter(unknown_param: str) -> List[str]:
    """
    Suggest similar parameter names from the registry.

    Args:
        unknown_param: Unknown parameter name.

    Returns:
        List of similar valid parameter names.
    """
    # Import here to avoid circular import (registry imports definitions which may use suggest)
    from .registry import REGISTRY  # pylint: disable=import-outside-toplevel

    return suggest_similar(unknown_param, REGISTRY.all_params.keys())


def suggest_constraint_key(unknown_key: str) -> List[str]:
    """
    Suggest similar constraint keys.

    Valid constraint keys are: choices, min, max

    Args:
        unknown_key: Unknown constraint key.

    Returns:
        List of similar valid constraint keys.
    """
    valid_keys = {"choices", "min", "max"}
    return suggest_similar(unknown_key, valid_keys)


def suggest_dependency_key(unknown_key: str) -> List[str]:
    """
    Suggest similar top-level dependency keys.

    Valid keys are: when_true, when_set

    Args:
        unknown_key: Unknown dependency key.

    Returns:
        List of similar valid dependency keys.
    """
    valid_keys = {"when_true", "when_set"}
    return suggest_similar(unknown_key, valid_keys)


def suggest_condition_key(unknown_key: str) -> List[str]:
    """
    Suggest similar condition keys within dependencies.

    Valid keys are: requires, recommends

    Args:
        unknown_key: Unknown condition key.

    Returns:
        List of similar valid condition keys.
    """
    valid_keys = {"requires", "recommends"}
    return suggest_similar(unknown_key, valid_keys)


@lru_cache(maxsize=128)
def get_param_suggestions_cached(unknown_param: str) -> tuple:
    """
    Cached version of suggest_parameter for repeated lookups.

    Returns tuple for hashability in cache.
    """
    return tuple(suggest_parameter(unknown_param))


def unknown_param_error(param_name: str) -> str:
    """
    Create an error message for an unknown parameter with suggestions.

    Args:
        param_name: The unknown parameter name.

    Returns:
        Error message with "did you mean?" suggestion if available.
    """
    suggestions = suggest_parameter(param_name)
    suggestion_text = format_suggestion(suggestions)

    if suggestion_text:
        return f"Unknown parameter '{param_name}'. {suggestion_text}"
    return f"Unknown parameter '{param_name}'."


def invalid_key_error(
    context: str,
    invalid_key: str,
    valid_keys: Iterable[str],
) -> str:
    """
    Create an error message for an invalid key with suggestions.

    Args:
        context: Description of what the key is for (e.g., "constraint", "dependency").
        invalid_key: The invalid key that was used.
        valid_keys: The set of valid keys.

    Returns:
        Error message with valid keys listed and "did you mean?" if applicable.
    """
    valid_set = set(valid_keys)
    suggestions = suggest_similar(invalid_key, valid_set)
    suggestion_text = format_suggestion(suggestions)

    base_msg = f"Invalid {context} key '{invalid_key}'. Valid keys are: {sorted(valid_set)}"

    if suggestion_text:
        return f"{base_msg}. {suggestion_text}"
    return base_msg
