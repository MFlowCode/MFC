#!/usr/bin/env python3
"""
Generate physics constraints documentation from PHYSICS_DOCS metadata
and AST-extracted validation rules.

Produces docs/documentation/physics_constraints.md (Doxygen-compatible).
"""

import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set

HERE = Path(__file__).resolve().parent
CASE_VALIDATOR_PATH = HERE / "case_validator.py"

# Make the toolchain package importable
_toolchain_dir = str(HERE.parent)
if _toolchain_dir not in sys.path:
    sys.path.insert(0, _toolchain_dir)

from mfc.case_validator import PHYSICS_DOCS  # noqa: E402  pylint: disable=wrong-import-position
from mfc.params.ast_analyzer import (  # noqa: E402  pylint: disable=wrong-import-position
    Rule,
    analyze_case_validator,
)

# Canonical category ordering
CATEGORY_ORDER = [
    "Thermodynamic Constraints",
    "Mixture Constraints",
    "Domain and Geometry",
    "Velocity and Dimensional Consistency",
    "Model Equations",
    "Boundary Conditions",
    "Bubble Physics",
    "Feature Compatibility",
    "Numerical Schemes",
    "Acoustic Sources",
    "Post-Processing",
]

_SEVERITY_ICON = {
    "error": "- ",
    "warning": "- \\f$\\triangle\\f$ ",
}

# Regex to detect code-like tokens in validation messages.
# Matches (in order of priority):
#   1. Fortran-style accessors:  fluid_pp(i)%gamma, bc_y%beg
#   2. Known short param names with optional array index: alpha(j), vel(2), Re(1), nb, dt
#   3. Snake_case identifiers with optional array index: model_eqns, weno_order, alpha_rho(j)
#   4. Single-letter grid dimension params: m, n, p
_CODE_RE = re.compile(
    r"(?<!`)\b("
    r"\w+(?:\([^)]*\))?%\w+(?:\([^)]*\))?"
    r"|(?:alpha|vel|Re|dt|nb|sigma|mhd|igr|Bx0|viscous|thermal|polytropic|relativity"
    r"|rhoref|pref|var|wavelength|npulse|hypoelasticity)(?:\([^)]*\))?"
    r"|[a-z]\w*_\w+(?:\([^)]*\))?"
    r"|[mnp]"
    r")(?!\w)"
)


def _format_message(msg: str) -> str:
    """Format a validation message for Doxygen-compatible markdown.

    Handles MFC/Fortran naming conventions:
    - Cleans AST artifacts: istr→i, jstr→j (f-string index variables)
    - Strips runtime value placeholders: (got varname)
    - Wraps code-like parameter references in backticks
    - Escapes Fortran ``%`` accessor for Doxygen (``\\%`` → literal %)
    """
    # Clean AST extraction artifacts: istr/jstr are f-string index vars
    msg = re.sub(r"\bistr\b", "i", msg)
    msg = re.sub(r"\bjstr\b", "j", msg)

    # Remove "(got <varname>)" — Python variable names, not useful in docs
    msg = re.sub(r"\s*\(got \w+\)", "", msg)

    # Clean f-string value placeholders that appear as bare variable names.
    # e.g. "= gamma implies" → "implies", "for pulse = pulse" → "for the given pulse"
    msg = re.sub(r"= (\w+) implies physical (\w+) = (\w+)", r"implies physical \2", msg)
    msg = re.sub(r"= (vel2|vel3) but", "is nonzero but", msg)
    msg = re.sub(r"for (support|pulse) = \1\b", r"for the given \1", msg)
    msg = re.sub(r"for (support|pulse) = (\d+)", r"for \1 = \2", msg)

    # Wrap code-like tokens in code formatting.
    # Doxygen consumes %<word> as "suppress auto-link" even inside backtick
    # and <code> spans. For Fortran % accessors, use split <code> tags with
    # CSS classes that remove internal borders, making them visually seamless.
    def _wrap_code(match: re.Match) -> str:
        token = match.group(1)
        if "%" in token:
            parts = token.split("%")
            result = parts[0]
            for part in parts[1:]:
                result += '%</code><code class="f90r">' + part
            return '<code class="f90l">' + result + "</code>"
        return f"`{token}`"

    msg = _CODE_RE.sub(_wrap_code, msg)

    return msg


def _stages_str(stages: Set[str]) -> str:
    order = ["common", "pre_process", "simulation", "post_process"]
    ordered = [s for s in order if s in stages]
    if ordered:
        return ", ".join(ordered)
    if stages:
        return ", ".join(sorted(stages))
    return "all"


def _severity_badge(rules: List[Rule]) -> str:
    severities = {r.severity for r in rules}
    if "error" in severities and "warning" in severities:
        return "error + warning"
    if "warning" in severities:
        return "warning"
    return "error"


def _collect_stages(rules: List[Rule]) -> Set[str]:
    stages: Set[str] = set()
    for r in rules:
        stages |= r.stages
    return stages


def _render_method(doc: dict, method_rules: List[Rule], lines: List[str]) -> None:
    """Render one PHYSICS_DOCS entry + its AST-extracted rules."""
    lines.append(f"### {doc['title']}\n")

    if "math" in doc:
        lines.append(f"\\f[{doc['math']}\\f]\n")

    lines.append(f"{doc['explanation']}\n")

    if method_rules:
        stages = _collect_stages(method_rules)
        lines.append(f"**Stage:** {_stages_str(stages)} | **Severity:** {_severity_badge(method_rules)}\n")

        seen: Set[str] = set()
        msgs = []
        for r in method_rules:
            if r.message not in seen:
                seen.add(r.message)
                msgs.append(r)
        if msgs:
            lines.append("**Enforced checks:**\n")
            for m in msgs:
                lines.append(f"{_SEVERITY_ICON.get(m.severity, '- ')}{_format_message(m.message)}")
            lines.append("")

    if "exceptions" in doc:
        lines.append("**Exceptions** (constraint does not apply):\n")
        for exc in doc["exceptions"]:
            lines.append(f"- {exc}")
        lines.append("")

    if "references" in doc:
        cites = ", ".join(f"\\cite {r}" for r in doc["references"])
        lines.append(f"**References:** {cites}\n")


    # Undocumented checks are omitted — they are discoverable via
    # @ref case_constraints "Case Creator Guide".


def render(rules: List[Rule]) -> str:
    """Render physics constraints page from PHYSICS_DOCS + AST rules."""
    by_method: Dict[str, List[Rule]] = defaultdict(list)
    for r in rules:
        by_method[r.method].append(r)

    by_category: Dict[str, List[str]] = defaultdict(list)
    for method, doc in PHYSICS_DOCS.items():
        by_category[doc["category"]].append(method)

    lines: List[str] = []
    lines.append("@page physics_constraints Physics Constraints\n")
    lines.append("# Physics Constraints Reference\n")
    lines.append(
        "> Auto-generated from `PHYSICS_DOCS` in `case_validator.py` and "
        "AST-extracted validation rules. Do not edit by hand.\n"
    )
    lines.append(
        "This document catalogs the physics constraints enforced by MFC's case parameter validator. "
        "Constraints are organized by physical category with mathematical justifications.\n"
    )
    lines.append(
        "For parameter syntax and allowed values, see @ref case \"Case Files\" and "
        "the @ref parameters \"Case Parameters\" reference. "
        "For feature compatibility and working examples, see "
        "@ref case_constraints \"Case Creator Guide\".\n"
    )

    extra_categories = [c for c in by_category if c not in CATEGORY_ORDER]
    for category in CATEGORY_ORDER + sorted(extra_categories):
        methods = by_category.get(category)
        if not methods:
            continue
        lines.append("---\n")
        lines.append(f"## {category}\n")
        for method in methods:
            _render_method(PHYSICS_DOCS[method], by_method.get(method, []), lines)

    return "\n".join(lines)


def main() -> None:
    analysis = analyze_case_validator(CASE_VALIDATOR_PATH)
    print(render(analysis["rules"]))


if __name__ == "__main__":
    main()
