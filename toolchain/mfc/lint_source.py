"""Static analysis for Fortran/Fypp source code.

Checks for patterns that indicate copy-paste bugs, non-standard constructs,
and hardcoded assumptions that break under different build configurations.
"""

import re
import sys
from pathlib import Path

# Source directory to scan (relative to repo root)
SRC_DIR = "src"

# Minimum stripped line length to consider for duplicate detection.
# Lines shorter than this (e.g. "end if", "end do") are ignored.
MIN_DUP_LINE_LEN = 40


def _is_comment_or_blank(stripped: str) -> bool:
    """True if stripped line is blank, a Fortran comment, or a Fypp directive."""
    return not stripped or stripped.startswith("!") or stripped.startswith("#:")


def _fortran_fpp_files(src_dir: Path):
    """Yield all .f90 and .fpp files under src/."""
    yield from sorted(src_dir.rglob("*.f90"))
    yield from sorted(src_dir.rglob("*.fpp"))


def _check_single_fypp_list(full_line: str, rel: Path, start_line: int) -> list[str]:
    """Parse one Fypp ``#:for ... in [...]`` line and return errors for duplicates."""
    errors: list[str] = []

    bracket_start = full_line.find("[")
    bracket_end = full_line.rfind("]")
    if not 0 <= bracket_start < bracket_end:
        return errors

    list_content = full_line[bracket_start + 1:bracket_end]
    list_content = list_content.replace("&", "")

    # Extract single- or double-quoted entries
    entries = re.findall(r"['\"]([^'\"]*)['\"]", list_content)

    seen: dict[str, int] = {}
    for pos, entry in enumerate(entries, 1):
        if entry in seen:
            errors.append(
                f"  {rel}:{start_line} Fypp list has duplicate"
                f" entry '{entry}' (positions {seen[entry]}"
                f" and {pos})."
                " Fix: one copy is likely a typo for a"
                " different variable"
            )
        else:
            seen[entry] = pos

    return errors


def check_fypp_list_duplicates(repo_root: Path) -> list[str]:
    """Check for duplicate entries in Fypp ``#:for VAR in [...]`` lists.

    Copy-paste errors in broadcast lists or loop variable lists can silently
    skip a variable while broadcasting another one twice.
    """
    errors: list[str] = []
    src_dir = repo_root / SRC_DIR

    for fpp in sorted(src_dir.rglob("*.fpp")):
        lines = fpp.read_text(encoding="utf-8").splitlines()
        rel = fpp.relative_to(repo_root)

        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith("#:for") and " in " in line and "[" in line:
                start_line = i + 1  # 1-indexed for display

                # Accumulate across Fortran-style '&' continuation lines
                full = line
                while full.rstrip().endswith("&") and i + 1 < len(lines):
                    i += 1
                    full += " " + lines[i].strip()

                errors.extend(_check_single_fypp_list(full, rel, start_line))
            i += 1

    return errors


def check_duplicate_lines(repo_root: Path) -> list[str]:
    """Flag identical adjacent non-trivial source lines.

    Exact duplicate consecutive lines are almost always copy-paste errors:
    a duplicated accumulation, a repeated subroutine argument, etc.
    """
    errors: list[str] = []
    src_dir = repo_root / SRC_DIR

    for src in _fortran_fpp_files(src_dir):
        lines = src.read_text(encoding="utf-8").splitlines()
        rel = src.relative_to(repo_root)

        prev_stripped = ""
        for i, line in enumerate(lines):
            stripped = line.strip()
            if (
                stripped == prev_stripped
                and len(stripped) >= MIN_DUP_LINE_LEN
                and not _is_comment_or_blank(stripped)
            ):
                display = stripped[:72]
                if len(stripped) > 72:
                    display += "..."
                errors.append(
                    f"  {rel}:{i + 1} identical to previous line:"
                    f" '{display}'."
                    " Fix: check for accidental copy-paste"
                )
            prev_stripped = stripped

    return errors


def check_hardcoded_byte_size(repo_root: Path) -> list[str]:
    """Flag ``int(8._wp, ...)`` patterns that assume 8-byte reals.

    When MFC is built in single precision (``wp = real32``), reals are
    4 bytes. Hard-coding 8 makes MPI I/O read/write the wrong amount.
    Use ``storage_size(0._stp)/8`` instead.
    """
    errors: list[str] = []
    src_dir = repo_root / SRC_DIR
    byte_re = re.compile(r"\bint\s*\(\s*8\._wp\b", re.IGNORECASE)

    for src in _fortran_fpp_files(src_dir):
        lines = src.read_text(encoding="utf-8").splitlines()
        rel = src.relative_to(repo_root)

        for i, line in enumerate(lines):
            stripped = line.strip()
            if _is_comment_or_blank(stripped):
                continue
            if byte_re.search(stripped.split("!")[0]):
                errors.append(
                    f"  {rel}:{i + 1} hard-codes 8-byte real size."
                    " Fix: use 'storage_size(0._stp)/8' instead of"
                    " '8._wp'"
                )

    return errors


def main():
    repo_root = Path(__file__).resolve().parents[2]

    all_errors: list[str] = []
    all_errors.extend(check_fypp_list_duplicates(repo_root))
    all_errors.extend(check_duplicate_lines(repo_root))
    all_errors.extend(check_hardcoded_byte_size(repo_root))

    if all_errors:
        print("Fortran/Fypp source analysis failed:")
        for e in all_errors:
            print(e)
        sys.exit(1)


if __name__ == "__main__":
    main()
