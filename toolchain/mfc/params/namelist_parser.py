"""
Parse Fortran namelist definitions to extract valid parameters for each target.

This module reads the Fortran source files and extracts the parameter names
from each target's namelist definition. This ensures the Python toolchain
stays in sync with what the Fortran code actually accepts.

When Fortran sources are unavailable (e.g. Homebrew installs), the fallback
is computed from NAMELIST_VARS in definitions.py.
"""

import re
from pathlib import Path
from typing import Dict, Optional, Set


def _fallback_params() -> Dict[str, Set[str]]:
    # Lazy import avoids circular dependency (definitions -> namelist_parser).
    from .definitions import NAMELIST_VARS
    from .generators.fortran_gen import TARGETS

    return {full: {v for v, ts in NAMELIST_VARS.items() if short in ts} for short, full in TARGETS}


def parse_namelist_from_file(filepath: Path) -> Set[str]:
    """Parse parameter names from the namelist /user_inputs/ block in an fpp file."""
    from .generators.fortran_gen import resolve_namelist_content

    content = resolve_namelist_content(filepath)

    namelist_match = re.search(
        r"namelist\s+/user_inputs/\s*(.+?)(?=\n\s*\n|\n\s*!(?!\s*&)|\n\s*[a-zA-Z_]+\s*=|$)",
        content,
        re.DOTALL | re.IGNORECASE,
    )
    if not namelist_match:
        raise ValueError(f"Could not find namelist /user_inputs/ in {filepath}")

    namelist_text = namelist_match.group(1)
    namelist_text = re.sub(r"&\s*\n\s*", " ", namelist_text)
    namelist_text = re.sub(r"#:.*", "", namelist_text)
    namelist_text = re.sub(r"!.*", "", namelist_text)

    found_params = set()
    for match in re.finditer(r"\b([a-zA-Z_][a-zA-Z0-9_]*)\b", namelist_text):
        name = match.group(1)
        if name.lower() not in {"namelist", "user_inputs", "if", "endif", "not"}:
            found_params.add(name)
    return found_params


def parse_all_namelists(mfc_root: Path) -> Dict[str, Set[str]]:
    """Parse namelist definitions from all MFC targets.

    Falls back to NAMELIST_VARS when Fortran sources are unavailable.
    """
    src = mfc_root / "src"
    target_files = {
        "pre_process": src / "pre_process" / "m_start_up.fpp",
        "simulation": src / "simulation" / "m_start_up.fpp",
        "post_process": src / "post_process" / "m_start_up.fpp",
    }

    for filepath in target_files.values():
        if not filepath.exists():
            return _fallback_params()

    return {name: parse_namelist_from_file(path) for name, path in target_files.items()}


def parse_fortran_constants(filepath: Path) -> Dict[str, int]:
    """Parse integer parameter constants from a Fortran source file."""
    constants: Dict[str, int] = {}
    pattern = re.compile(r"integer\s*,\s*parameter\s*::\s*(\w+)\s*=\s*(\d+)", re.IGNORECASE)
    try:
        text = filepath.read_text()
    except FileNotFoundError:
        return constants
    for m in pattern.finditer(text):
        constants[m.group(1)] = int(m.group(2))
    return constants


_FORTRAN_CONSTANTS_CACHE: Optional[Dict[str, int]] = None


def get_fortran_constants() -> Dict[str, int]:
    """Get Fortran compile-time constants from m_constants.fpp.

    Cached after first call. Returns an empty dict when src/ is unavailable;
    callers supply inline defaults via _fc(name, default) in definitions.py.
    """
    global _FORTRAN_CONSTANTS_CACHE  # noqa: PLW0603
    if _FORTRAN_CONSTANTS_CACHE is None:
        path = get_mfc_root() / "src" / "common" / "m_constants.fpp"
        _FORTRAN_CONSTANTS_CACHE = parse_fortran_constants(path)
    return _FORTRAN_CONSTANTS_CACHE


def get_mfc_root() -> Path:
    """Return the MFC root directory (4 levels above this file)."""
    return Path(__file__).resolve().parent.parent.parent.parent


_TARGET_PARAMS_CACHE: Dict[str, Set[str]] = {}


def get_target_params() -> Dict[str, Set[str]]:
    """Return valid parameters per target, parsing Fortran sources if needed."""
    if not _TARGET_PARAMS_CACHE:
        _TARGET_PARAMS_CACHE.update(parse_all_namelists(get_mfc_root()))
    return _TARGET_PARAMS_CACHE


def is_param_valid_for_target(param_name: str, target_name: str) -> bool:
    """Return True if param_name is valid for target_name.

    Handles indexed params (patch_icpp(1)%geometry) by checking the base name.
    """
    valid_params = get_target_params().get(target_name, set())
    base_match = re.match(r"^([a-zA-Z_][a-zA-Z0-9_]*)", param_name)
    return base_match.group(1) in valid_params if base_match else param_name in valid_params


if __name__ == "__main__":
    import sys

    try:
        parsed_targets = parse_all_namelists(get_mfc_root())

        for tgt, tgt_params in sorted(parsed_targets.items()):
            print(f"{tgt}: {len(tgt_params)} parameters")
            for param in sorted(tgt_params)[:10]:
                print(f"  - {param}")
            if len(tgt_params) > 10:
                print(f"  ... and {len(tgt_params) - 10} more")
            print()

        common = set.intersection(*parsed_targets.values())
        print(f"Parameters in ALL targets ({len(common)}): {sorted(common)[:20]}...")

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
