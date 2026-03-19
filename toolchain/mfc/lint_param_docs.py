"""Check parameter documentation completeness.

Tier 1 (blocking): Parameters with DESCRIPTIONS entries must appear in case.md.
Tier 2 (warnings): Fortran namelist / REGISTRY sync checks.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


def _import_registry(repo_root: Path):
    """Import REGISTRY from the toolchain."""
    toolchain_dir = str(repo_root / "toolchain")
    if toolchain_dir not in sys.path:
        sys.path.insert(0, toolchain_dir)
    from mfc.params import REGISTRY
    from mfc.params.descriptions import DESCRIPTIONS

    return REGISTRY, DESCRIPTIONS


def _param_to_base_pattern(name: str) -> str:
    """Strip numeric indexes: patch_icpp(1)%vel(2) -> patch_icpp%vel."""
    return re.sub(r"\(\d+(?:,\s*-?\d+)*\)", "", name)


def _build_case_md_index(text: str) -> tuple[set[str], str]:
    """Return (backtick_tokens, normalized_text) with %% -> %."""
    normalized = text.replace("%%", "%")
    tokens = set()
    for m in re.finditer(r"``([^`\n]+)``|`([^`\n]+)`", normalized):
        span = (m.group(1) or m.group(2)).strip()
        tokens.add(span)
    return tokens, normalized


def _param_appears_in_case_md(param_base: str, tokens: set[str], text: str) -> bool:
    """Check if param appears in case.md (lenient: handles bracket shorthand)."""
    if param_base in tokens or f"`{param_base}`" in text:
        return True
    if re.search(rf"`[^`]*\b{re.escape(param_base)}\b[^`]*`", text):
        return True
    # family%attr: check if attr appears standalone (patch table shorthand)
    if "%" in param_base:
        attr = param_base.split("%", 1)[1]
        if f"`{attr}`" in text or f"`{attr}(" in text:
            return True
    # x/y/z axis variants: bc_x%beg matches bc_[x,y,z]%beg[end]
    for axis in ("x", "y", "z"):
        if param_base.startswith(f"{axis}_") or param_base.startswith(f"{axis}%"):
            rest = param_base[len(axis) :]
            for other in ("x", "y", "z"):
                if other != axis:
                    sib = other + rest
                    if re.search(rf"`[^`]*\b{re.escape(sib)}\b[^`]*`", text):
                        return True
            if re.search(rf"`[^`]*\[.*{re.escape(axis)}.*\]{re.escape(rest)}", text):
                return True
        # stretch_x -> stretch_x[y,z] or bc_x%beg -> bc_[x,y,z]%beg
        m = re.match(rf"^(.+_){axis}(.*)$", param_base)
        if m:
            pfx, sfx = m.group(1), m.group(2)
            if re.search(rf"`{re.escape(pfx)}\[.*{re.escape(axis)}.*\]{re.escape(sfx)}", text):
                return True
            for other in ("x", "y", "z"):
                if other != axis:
                    sib = pfx + other + sfx
                    if re.search(rf"`[^`]*\b{re.escape(sib)}\b[^`]*`", text):
                        return True
    # %beg/%end shorthand: param%beg matches if param%end is present
    if param_base.endswith("%beg") or param_base.endswith("%end"):
        stem = param_base.rsplit("%", 1)[0]
        other_sfx = "%end" if param_base.endswith("%beg") else "%beg"
        if re.search(rf"`[^`]*\b{re.escape(stem + other_sfx)}\b[^`]*`", text):
            return True
    return False


def _parse_namelist_params(fpp_path: Path) -> set[str]:
    """Parse parameter names from a namelist /user_inputs/ block in an fpp file."""
    text = fpp_path.read_text(encoding="utf-8")
    params = set()

    in_namelist = False
    accum = ""
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("!") or stripped.startswith("#:") or stripped.startswith("#"):
            continue
        if "namelist /user_inputs/" in stripped.lower():
            in_namelist = True
            accum = stripped.split("/user_inputs/", 1)[1] if "/user_inputs/" in stripped else ""
            continue
        if in_namelist:
            if stripped.startswith("&"):
                stripped = stripped[1:]
            accum += " " + stripped
            if not accum.rstrip().endswith("&"):
                in_namelist = False

    accum = accum.replace("&", " ")
    for raw_token in accum.split(","):
        name = raw_token.strip()
        if name and re.match(r"^[a-zA-Z_]\w*$", name):
            params.add(name)

    return params


def check_descriptions_in_case_md(repo_root: Path) -> list[str]:
    """Tier 1, Check 1+2: Params with DESCRIPTIONS entries should appear in case.md."""
    REGISTRY, DESCRIPTIONS = _import_registry(repo_root)

    KNOWN_UNDOCUMENTED: set[str] = set()

    case_md_path = repo_root / "docs" / "documentation" / "case.md"
    if not case_md_path.exists():
        return ["  docs/documentation/case.md not found"]

    case_md_text = case_md_path.read_text(encoding="utf-8")
    tokens, normalized = _build_case_md_index(case_md_text)

    errors = []
    checked = set()
    for param_name in sorted(DESCRIPTIONS.keys()):
        base = _param_to_base_pattern(param_name)
        if base in checked:
            continue
        checked.add(base)

        if base in KNOWN_UNDOCUMENTED:
            continue

        if not _param_appears_in_case_md(base, tokens, normalized):
            errors.append(f"  docs/documentation/case.md is missing documentation for parameter '{base}' (has description in descriptions.py). Fix: add parameter to case.md")

    return errors


def check_scalar_params_have_descriptions(repo_root: Path) -> list[str]:
    """Tier 3: Scalar (non-indexed, non-family) REGISTRY params must have descriptions."""
    REGISTRY, DESCRIPTIONS = _import_registry(repo_root)

    desc_bases = set()
    for k in DESCRIPTIONS.keys():
        desc_bases.add(_param_to_base_pattern(k))

    errors = []
    for name in sorted(REGISTRY.all_params.keys()):
        if "(" in name or "%" in name:
            continue
        if name not in desc_bases:
            errors.append(f"  REGISTRY param '{name}' has no description in descriptions.py. Fix: add entry to descriptions.py")

    return errors


def _is_derived_type_parent_or_field(name: str, all_params: dict) -> bool:
    """Check if name is a derived type parent (has children) or bare field."""
    # Parent: REGISTRY has params like name%foo or name(1)%foo
    if any(k.startswith(f"{name}%") or k.startswith(f"{name}(") for k in all_params):
        return True
    # Bare field: REGISTRY has params like foo(1)%name or foo%name
    if any(k.endswith(f"%{name}") for k in all_params):
        return True
    return False


def check_namelist_registry_sync(repo_root: Path) -> list[str]:
    """Tier 2: Unknown Fortran namelist params must be in REGISTRY (blocking)."""
    REGISTRY, _ = _import_registry(repo_root)

    errors = []
    startup_files = sorted((repo_root / "src").rglob("m_start_up.fpp"))
    all_params = REGISTRY.all_params

    for fpp in startup_files:
        nl_params = _parse_namelist_params(fpp)
        rel = fpp.relative_to(repo_root)

        for p in sorted(nl_params):
            if REGISTRY.is_known_param(p):
                continue
            if _is_derived_type_parent_or_field(p, all_params):
                continue
            errors.append(f"  {rel} has namelist param '{p}' not in REGISTRY. Fix: add to definitions.py or remove from namelist")

    return errors


def main():
    repo_root = Path(__file__).resolve().parents[2]

    # Tier 1: described params must appear in case.md
    errors = check_descriptions_in_case_md(repo_root)

    # Tier 2: namelist params must be in REGISTRY
    errors.extend(check_namelist_registry_sync(repo_root))

    # Tier 3: scalar REGISTRY params must have descriptions
    errors.extend(check_scalar_params_have_descriptions(repo_root))

    if errors:
        print("Parameter documentation check failed:")
        for e in errors:
            print(e)
        sys.exit(1)


if __name__ == "__main__":
    main()
