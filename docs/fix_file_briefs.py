#!/usr/bin/env python3
"""Ensure @file/@brief headers match actual module/program declarations.

Usage: python3 fix_file_briefs.py [source_dir]
  source_dir defaults to current directory.

For each .fpp/.f90 in src/{pre_process,simulation,post_process,common}:
  1. Parses the first `module <name>` or `program <name>` declaration.
  2. If no @file directive exists in the first 15 lines, prepends a header.
  3. If a "Contains module/program ..." @brief exists, rewrites the name
     to match the source, using @ref for mixed-case Fortran identifiers
     (Doxygen lowercases Fortran namespaces).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

src_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")

DIRS = [
    src_dir / "src" / "pre_process",
    src_dir / "src" / "simulation",
    src_dir / "src" / "post_process",
    src_dir / "src" / "common",
]

# First `module X` or `program X` that isn't `end module/program`.
DECL_RE = re.compile(
    r"^\s*(module(?!\s+procedure)\b|program)\s+(\w+)\s*$",
    re.MULTILINE | re.IGNORECASE,
)

# Any "Contains module/program <name>" in a Doxygen comment line.
# Matches both `!! @brief Contains ...` and `!> @brief Contains ...`.
BRIEF_CONTAINS_RE = re.compile(
    r"^(!!|!>)\s*@brief\s+(Contains (?:module|program) )(.*)",
    re.MULTILINE,
)


def find_entity(text: str) -> tuple[str, str] | None:
    """Return (kind, name) of the first module/program declaration."""
    for m in DECL_RE.finditer(text):
        line_start = text.rfind("\n", 0, m.start()) + 1
        line = text[line_start : m.end()].strip()
        if line.lower().startswith("end"):
            continue
        return m.group(1).lower(), m.group(2)
    return None


def make_ref(name: str) -> str:
    """Return @ref for mixed-case names, plain name for lowercase."""
    lower = name.lower()
    return f'@ref {lower} "{name}"' if lower != name else name


def has_file_directive(text: str) -> bool:
    """Check if @file appears in the first 15 lines."""
    head = "\n".join(text.splitlines()[:15])
    return bool(re.search(r"@file", head))


fixed = 0
for d in DIRS:
    if not d.exists():
        continue
    for f in sorted(list(d.glob("*.fpp")) + list(d.glob("*.f90"))):
        text = f.read_text(encoding="utf-8")
        entity = find_entity(text)
        if entity is None:
            continue
        kind, name = entity
        ref = make_ref(name)

        if not has_file_directive(text):
            # No @file at all — prepend a complete header.
            header = f"!>\n!! @file\n!! @brief Contains {kind} {ref}\n\n"
            text = header + text
            f.write_text(text, encoding="utf-8")
            fixed += 1
            print(f"Added  {f.relative_to(src_dir)}")
            continue

        # Has @file — check if there's a "Contains module/program" brief to fix.
        m = BRIEF_CONTAINS_RE.search(text)
        if m is None:
            continue  # Has @file but no "Contains ..." brief — leave it alone

        current_name = m.group(3).strip()
        if current_name == ref:
            continue  # Already correct

        # Replace the name portion of the brief line.
        new_line = f"{m.group(1)} @brief {m.group(2)}{ref}"
        new_text = text[: m.start()] + new_line + text[m.end() :]

        if new_text != text:
            f.write_text(new_text, encoding="utf-8")
            fixed += 1
            print(f"Fixed  {f.relative_to(src_dir)}: {current_name} -> {ref}")

print(f"Done — {fixed} file(s) updated.")
