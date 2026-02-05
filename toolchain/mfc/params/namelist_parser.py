"""
Parse Fortran namelist definitions to extract valid parameters for each target.

This module reads the Fortran source files and extracts the parameter names
from each target's namelist definition. This ensures the Python toolchain
stays in sync with what the Fortran code actually accepts.
"""

import re
from pathlib import Path
from typing import Dict, Set


def parse_namelist_from_file(filepath: Path) -> Set[str]:
    """
    Parse a Fortran file and extract parameter names from the namelist definition.

    Args:
        filepath: Path to the Fortran source file (m_start_up.fpp)

    Returns:
        Set of parameter names found in the namelist
    """
    content = filepath.read_text()

    # Find the namelist block - starts with "namelist /user_inputs/"
    # and continues until a line without continuation (&) or a blank line
    namelist_match = re.search(
        r'namelist\s+/user_inputs/\s*(.+?)(?=\n\s*\n|\n\s*!(?!\s*&)|\n\s*[a-zA-Z_]+\s*=)',
        content,
        re.DOTALL | re.IGNORECASE
    )

    if not namelist_match:
        raise ValueError(f"Could not find namelist /user_inputs/ in {filepath}")

    namelist_text = namelist_match.group(1)

    # Remove Fortran line continuations (&) and join lines
    namelist_text = re.sub(r'&\s*\n\s*', ' ', namelist_text)

    # Remove preprocessor directives (#:if, #:endif, etc.)
    namelist_text = re.sub(r'#:.*', '', namelist_text)

    # Remove comments (! to end of line, but not inside strings)
    namelist_text = re.sub(r'!.*', '', namelist_text)

    # Extract parameter names - they're comma-separated identifiers
    # Parameter names are alphanumeric with underscores
    found_params = set()
    for match in re.finditer(r'\b([a-zA-Z_][a-zA-Z0-9_]*)\b', namelist_text):
        name = match.group(1)
        # Skip Fortran keywords that might appear
        if name.lower() not in {'namelist', 'user_inputs', 'if', 'endif', 'not'}:
            found_params.add(name)

    return found_params


def parse_all_namelists(mfc_root: Path) -> Dict[str, Set[str]]:
    """
    Parse namelist definitions from all MFC targets.

    Args:
        mfc_root: Path to MFC root directory

    Returns:
        Dict mapping target name to set of valid parameter names
    """
    targets = {
        'pre_process': mfc_root / 'src' / 'pre_process' / 'm_start_up.fpp',
        'simulation': mfc_root / 'src' / 'simulation' / 'm_start_up.fpp',
        'post_process': mfc_root / 'src' / 'post_process' / 'm_start_up.fpp',
    }

    result = {}
    for target_name, filepath in targets.items():
        if not filepath.exists():
            raise FileNotFoundError(f"Fortran source not found: {filepath}")
        result[target_name] = parse_namelist_from_file(filepath)

    return result


def get_mfc_root() -> Path:
    """Get the MFC root directory from this file's location."""
    # This file is at toolchain/mfc/params/namelist_parser.py
    # MFC root is 4 levels up
    return Path(__file__).resolve().parent.parent.parent.parent


# Module-level cache for parsed target params
_TARGET_PARAMS_CACHE: Dict[str, Set[str]] = {}


def get_target_params() -> Dict[str, Set[str]]:
    """
    Get the valid parameters for each target, parsing Fortran if needed.

    Returns:
        Dict mapping target name to set of valid parameter names
    """
    if not _TARGET_PARAMS_CACHE:
        _TARGET_PARAMS_CACHE.update(parse_all_namelists(get_mfc_root()))
    return _TARGET_PARAMS_CACHE


def is_param_valid_for_target(param_name: str, target_name: str) -> bool:
    """
    Check if a parameter is valid for a given target.

    This handles both scalar params (like "m") and indexed params
    (like "patch_icpp(1)%geometry") by checking the base name.

    Args:
        param_name: The parameter name (may include indices like "(1)%attr")
        target_name: One of 'pre_process', 'simulation', 'post_process'

    Returns:
        True if the parameter is valid for the target
    """
    valid_params = get_target_params().get(target_name, set())

    # Extract base parameter name (before any index or attribute)
    # e.g., "patch_icpp(1)%geometry" -> "patch_icpp"
    # e.g., "fluid_pp(2)%gamma" -> "fluid_pp"
    base_match = re.match(r'^([a-zA-Z_][a-zA-Z0-9_]*)', param_name)
    if base_match:
        return base_match.group(1) in valid_params

    return param_name in valid_params


if __name__ == '__main__':
    # Test the parser
    import sys

    try:
        parsed_targets = parse_all_namelists(get_mfc_root())

        print("Parsed namelist parameters:\n")
        for tgt, tgt_params in sorted(parsed_targets.items()):
            print(f"{tgt}: {len(tgt_params)} parameters")
            # Print first 10 as sample
            sorted_list = sorted(tgt_params)
            for param in sorted_list[:10]:
                print(f"  - {param}")
            if len(tgt_params) > 10:
                print(f"  ... and {len(tgt_params) - 10} more")
            print()

        # Show params unique to each target
        print("Parameters unique to each target:\n")
        all_param_names = set.union(*parsed_targets.values())
        for tgt, tgt_params in sorted(parsed_targets.items()):
            other = set.union(*[p for t, p in parsed_targets.items() if t != tgt])
            unique = tgt_params - other
            print(f"{tgt} only ({len(unique)}): {sorted(unique)[:15]}...")
            print()

        # Show params in all targets
        common = set.intersection(*parsed_targets.values())
        print(f"Parameters in ALL targets ({len(common)}): {sorted(common)[:20]}...")

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
