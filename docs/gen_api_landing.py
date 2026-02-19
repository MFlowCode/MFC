#!/usr/bin/env python3
"""Generate API landing pages for pre_process, simulation, and post_process.

Usage: python3 gen_api_landing.py [source_dir]
  source_dir defaults to current directory.

Scans src/{target}/*.fpp and src/common/*.fpp to produce module tables
in docs/{target}/readme.md. Intro text is defined below per target.
"""

import sys
from pathlib import Path

src_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")

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


def get_modules(directory: Path) -> list[str]:
    """Return sorted list of module names (m_*) from .fpp files."""
    return sorted(
        f.stem for f in directory.glob("m_*.fpp")
    )


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
        for m in target_modules:
            lines.append(f"- @ref {m.lower()} \"{m}\"")
        lines.append("")

    if common_modules:
        lines.append("### Common (shared)")
        lines.append("")
        for m in common_modules:
            lines.append(f"- @ref {m.lower()} \"{m}\"")
        lines.append("")

    lines.append("## See Also")
    lines.append("")
    lines.append("- [Home & User Guide](../documentation/index.html)")
    for label, sib in info["siblings"]:
        lines.append(f"- [{label} API](../{sib}/index.html)")
    lines.append("")

    out.write_text("\n".join(lines))
    print(f"Generated {out}")
