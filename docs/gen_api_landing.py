#!/usr/bin/env python3
"""Generate API landing pages for MFC documentation.

Usage: python3 gen_api_landing.py [source_dir]
  source_dir defaults to current directory.

Scans src/{target}/*.fpp,*.f90 and src/common/*.fpp,*.f90 to produce module
tables in docs/{target}/readme.md. Also generates a unified API landing page
at docs/api/readme.md that links to all three sub-project APIs.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

src_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")

_BRIEF_RE = re.compile(r"^!>\s*@brief\s+(.+)", re.IGNORECASE)

TARGETS = {
    "pre_process": {
        "title": "MFC Pre-Process",
        "intro": (
            "The pre-process component generates initial conditions and "
            "computational meshes for MFC simulations. It supports patch-based "
            "geometry construction, multi-component material initialization, "
            "and immersed boundary geometry."
        ),
        "siblings": [("Simulation", "simulation"), ("Post-Process", "post_process")],
    },
    "simulation": {
        "title": "MFC Simulation",
        "intro": (
            "The simulation component is the core flow solver. It advances the "
            "governing equations in time using high-order finite-volume methods "
            "on structured grids with GPU acceleration via OpenACC/OpenMP offloading."
        ),
        "siblings": [("Pre-Process", "pre_process"), ("Post-Process", "post_process")],
    },
    "post_process": {
        "title": "MFC Post-Process",
        "intro": (
            "The post-process component reads raw simulation output and computes "
            "derived quantities for visualization. It produces silo/HDF5 files "
            "compatible with VisIt, ParaView, and other visualization tools."
        ),
        "siblings": [("Pre-Process", "pre_process"), ("Simulation", "simulation")],
    },
}


def _extract_brief(path: Path) -> str:
    """Extract the module-level @brief from a Fortran source file.

    Reads the first !> @brief line (outside the @file block) and any
    !! continuation lines, joins them, and returns the first sentence.
    """
    lines = path.read_text(encoding="utf-8").splitlines()
    parts: list[str] = []
    in_file_block = False
    collecting = False
    for line in lines[:40]:
        stripped = line.strip()
        # Track whether we are inside the @file doc-comment block
        if stripped.startswith("!>") and not stripped.startswith("!> @brief"):
            in_file_block = True
            continue
        if in_file_block and stripped.startswith("!!"):
            continue
        if in_file_block:
            in_file_block = False
        # Match module-level @brief (outside @file block)
        m = _BRIEF_RE.match(stripped)
        if m and not collecting:
            brief_text = m.group(1).strip()
            # Skip file-level "Contains module/program X" briefs
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
    # Take first sentence only
    dot = text.find(". ")
    if dot != -1:
        text = text[:dot]
    return text.rstrip(". ").strip()


def get_modules(directory: Path) -> list[tuple[str, str]]:
    """Return sorted list of (module_name, brief) from .fpp and .f90 files."""
    modules: dict[str, str] = {}
    for pattern in ("m_*.fpp", "m_*.f90"):
        for f in directory.glob(pattern):
            if f.stem in modules:
                continue
            modules[f.stem] = _extract_brief(f)
    return sorted(modules.items())


for target, info in TARGETS.items():
    target_modules = get_modules(src_dir / "src" / target)
    common_modules = get_modules(src_dir / "src" / "common")

    out = src_dir / "docs" / target / "readme.md"

    lines = [
        f"@mainpage {info['title']}",
        "",
        info["intro"],
        "",
        "## Modules",
        "",
    ]

    if target_modules:
        lines.append(f"### {info['title'].split()[-1]}")
        lines.append("")
        lines.append("| Module | Description |")
        lines.append("|--------|-------------|")
        for name, brief in target_modules:
            lines.append(f"| @ref {name.lower()} \"{name}\" | {brief} |")
        lines.append("")

    if common_modules:
        lines.append("### Common (shared)")
        lines.append("")
        lines.append("| Module | Description |")
        lines.append("|--------|-------------|")
        for name, brief in common_modules:
            lines.append(f"| @ref {name.lower()} \"{name}\" | {brief} |")
        lines.append("")

    lines.append("## See Also")
    lines.append("")
    lines.append("- [Home & User Guide](../documentation/index.html)")
    for label, sib in info["siblings"]:
        lines.append(f"- [{label} API](../{sib}/index.html)")
    lines.append("")

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"Generated {out}")

# --- Unified API landing page ---
api_lines = [
    "@mainpage API Documentation",
    "",
    "MFC's source code is organized into three components that form a complete "
    "simulation pipeline. Each component has full module-level API documentation.",
    "",
]

for target, info in TARGETS.items():
    label = info["title"].replace("MFC ", "")
    api_lines.append(f"### [{label}](../{target}/index.html)")
    api_lines.append("")
    api_lines.append(info["intro"])
    api_lines.append("")

api_lines.append(
    "All three components share a set of **common modules** for MPI communication, "
    "variable conversion, derived types, and utility functions. "
    "These are documented within each component's API reference."
)
api_lines.append("")

api_out = src_dir / "docs" / "api" / "readme.md"
api_out.parent.mkdir(parents=True, exist_ok=True)
api_out.write_text("\n".join(api_lines), encoding="utf-8")
print(f"Generated {api_out}")
