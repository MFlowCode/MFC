#!/usr/bin/env python3
"""Generate architecture.md from template + auto-generated module map.

Usage: python3 gen_architecture.py [source_dir]
  source_dir defaults to current directory.

Reads docs/documentation/architecture.md.in as a template, generates the
module map section from docs/module_categories.json and source file briefs,
and writes the final docs/documentation/architecture.md.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

src_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")

_BRIEF_RE = re.compile(r"^!>\s*@brief\s+(.+)", re.IGNORECASE)


def _extract_brief(path: Path) -> str:
    """Extract the module-level @brief from a Fortran source file."""
    lines = path.read_text(encoding="utf-8").splitlines()
    parts: list[str] = []
    in_file_block = False
    collecting = False
    for line in lines[:40]:
        stripped = line.strip()
        if stripped.startswith("!>") and not stripped.startswith("!> @brief"):
            in_file_block = True
            continue
        if in_file_block and stripped.startswith("!!"):
            continue
        if in_file_block:
            in_file_block = False
        m = _BRIEF_RE.match(stripped)
        if m and not collecting:
            brief_text = m.group(1).strip()
            if re.match(r"Contains\s+(module|program)\s+\w+", brief_text):
                continue
            parts.append(brief_text)
            collecting = True
            continue
        if collecting:
            if stripped.startswith("!!"):
                parts.append(stripped.lstrip("! ").strip())
            else:
                break
    if not parts:
        return ""
    text = " ".join(parts)
    dot = text.find(". ")
    if dot != -1:
        text = text[:dot]
    return text.rstrip(". ").strip()


def _find_module_file(name: str) -> Path | None:
    """Find the source file for a module name across all src/ directories."""
    for subdir in ("common", "pre_process", "simulation", "post_process"):
        for ext in (".fpp", ".f90", ".F90"):
            path = src_dir / "src" / subdir / (name + ext)
            if path.exists():
                return path
    return None


def generate_module_map() -> str:
    """Generate the module map markdown from categories + source briefs."""
    categories_file = src_dir / "docs" / "module_categories.json"
    categories = json.loads(categories_file.read_text(encoding="utf-8"))

    lines: list[str] = []
    for entry in categories:
        cat = entry["category"]
        modules = entry["modules"]

        lines.append(f"### {cat}")
        lines.append("| Module | Role |")
        lines.append("|--------|------|")

        for mod in modules:
            path = _find_module_file(mod)
            brief = _extract_brief(path) if path else ""
            lines.append(f"| `{mod}` | {brief} |")

        lines.append("")

    return "\n".join(lines)


def main() -> None:
    template = src_dir / "docs" / "documentation" / "architecture.md.in"
    output = src_dir / "docs" / "documentation" / "architecture.md"

    text = template.read_text(encoding="utf-8")
    module_map = generate_module_map()
    text = text.replace("<!-- MODULE_MAP -->", module_map)

    output.write_text(text, encoding="utf-8")
    print(f"Generated {output}")


if __name__ == "__main__":
    main()
