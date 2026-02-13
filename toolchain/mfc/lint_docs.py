"""Check that file paths, cite keys, parameters, and @ref targets in docs still exist."""

import re
import sys
from pathlib import Path

# Docs to scan for file path references
DOCS = [
    "docs/documentation/contributing.md",
    "docs/documentation/gpuParallelization.md",
    "docs/documentation/running.md",
    "docs/documentation/case.md",
    "docs/documentation/equations.md",
    ".github/copilot-instructions.md",
]

# Match backtick-wrapped strings that look like repo-relative file paths
PATH_RE = re.compile(r"`((?:src|toolchain|\.github|docs|examples|tests)/[^`]+)`")

# Skip paths with placeholders, globs, or patterns (not real file paths)
SKIP_RE = re.compile(r"[<>*()\[\]{}]|/\.\.\.|%|\$")

# Match \cite Key references in markdown (Doxygen-style)
CITE_RE = re.compile(r"\\cite\s+(\w+)")

# Match backtick-quoted parameter names like `param = value`, `param > 0`, or `param`
PARAM_RE = re.compile(r"`([a-z_][a-z0-9_%]*)(?:\s*[=><][^`]*)?\s*`")

# Parameter-like names to skip (not actual MFC parameters)
PARAM_SKIP = re.compile(
    r"^(src/|toolchain/|\.github|docs/|examples/|tests/)"  # file paths
    r"|^\.(?:true|false)\.$"   # Fortran logicals
    r"|^\d"                     # numeric values
    r"|^[A-Z]"                  # constants/types (uppercase start)
)

# Backtick tokens in case.md that are not real parameters (analytical shorthands,
# stress tensor component names, prose identifiers, hardcoded constants)
CASE_MD_SKIP = {
    # Analytical shorthand variables (stretching formulas, "Analytical Definition" table)
    "eps", "lx", "ly", "lz", "xc", "yc", "zc", "x_cb",
    # Stress tensor component names (descriptive, not params)
    "tau_xx", "tau_xy", "tau_xz", "tau_yy", "tau_yz", "tau_zz",
    # Prose identifiers (example names, math symbols)
    "scaling", "c_h", "thickness",
    # Hardcoded Fortran constants (not case-file params)
    "init_dir", "zeros_default",
}

# Docs to check for parameter references, with per-file skip sets
PARAM_DOCS = {
    "docs/documentation/equations.md": set(),
    "docs/documentation/case.md": CASE_MD_SKIP,
}

# Match @ref page_id patterns
REF_RE = re.compile(r"@ref\s+(\w+)")


def check_docs(repo_root: Path) -> list[str]:
    """Check that file paths referenced in documentation still exist."""
    errors = []
    for doc in DOCS:
        doc_path = repo_root / doc
        if not doc_path.exists():
            continue
        text = doc_path.read_text(encoding="utf-8")
        for match in PATH_RE.finditer(text):
            path_str = match.group(1)
            if SKIP_RE.search(path_str):
                continue
            # Strip trailing punctuation that may have leaked in
            path_str = path_str.rstrip(".,;:!?")
            if not (repo_root / path_str).exists():
                errors.append(f"  {doc} references '{path_str}' but it does not exist")
    return errors


def parse_bib_keys(bib_path: Path) -> set[str]:
    """Parse references.bib and return the set of valid citation keys (lowercase)."""
    text = bib_path.read_text(encoding="utf-8")
    return {m.group(1).lower() for m in re.finditer(r"@\w+\{(\w+),", text)}


def check_cite_keys(repo_root: Path) -> list[str]:
    """Check that \\cite keys in doc files exist in references.bib."""
    bib_path = repo_root / "docs" / "references.bib"
    if not bib_path.exists():
        return []

    valid_keys = parse_bib_keys(bib_path)
    errors = []

    doc_dir = repo_root / "docs" / "documentation"
    if not doc_dir.exists():
        return []

    for md_file in sorted(doc_dir.glob("*.md")):
        text = md_file.read_text(encoding="utf-8")
        rel = md_file.relative_to(repo_root)
        for match in CITE_RE.finditer(text):
            key = match.group(1)
            if key.lower() not in valid_keys:
                errors.append(f"  {rel} uses \\cite {key} but no bib entry found")

    return errors


def _strip_code_blocks(text: str) -> str:
    """Remove fenced code blocks (``` ... ```) from markdown text."""
    lines = text.split("\n")
    result = []
    in_block = False
    for line in lines:
        if line.strip().startswith("```"):
            in_block = not in_block
            continue
        if not in_block:
            result.append(line)
    return "\n".join(result)


def _is_valid_param(param: str, valid_params: set, sub_params: set) -> bool:
    """Check if a param name (possibly with %) is valid against REGISTRY."""
    if "(" in param or ")" in param:
        return True  # Skip indexed refs like patch_icpp(i)%vel(j)

    base = param.split("%")[0] if "%" in param else param

    if base in valid_params or base in sub_params:
        return True

    # Check sub-param part after %
    if "%" in param:
        sub = param.split("%")[-1]
        if sub in sub_params:
            return True

    # Family prefix check (only accept at a boundary to avoid prefix-typo matches)
    if any(p.startswith(base + b) for b in ("(", "%", "_") for p in valid_params):
        return True

    return False


def check_param_refs(repo_root: Path) -> list[str]:  # pylint: disable=too-many-locals
    """Check that parameter names in documentation exist in the MFC registry."""
    # Import REGISTRY from the toolchain
    toolchain_dir = str(repo_root / "toolchain")
    if toolchain_dir not in sys.path:
        sys.path.insert(0, toolchain_dir)
    try:
        from mfc.params import REGISTRY  # pylint: disable=import-outside-toplevel
    except ImportError:
        print("  Warning: could not import REGISTRY, skipping parameter check")
        return []

    valid_params = set(REGISTRY.all_params.keys())
    # Build set of sub-parameter base names (strip trailing (N) indexes)
    _sub_raw = {p.split("%")[-1] for p in valid_params if "%" in p}
    sub_params = set()
    for s in _sub_raw:
        sub_params.add(s)
        base = re.sub(r"\(\d+(?:,\s*\d+)*\)$", "", s)
        if base != s:
            sub_params.add(base)

    errors = []

    for doc_rel, extra_skip in PARAM_DOCS.items():
        doc_path = repo_root / doc_rel
        if not doc_path.exists():
            continue

        text = _strip_code_blocks(doc_path.read_text(encoding="utf-8"))
        seen = set()

        # Check plain params
        for match in PARAM_RE.finditer(text):
            param = match.group(1)
            if param in seen:
                continue
            seen.add(param)

            if PARAM_SKIP.search(param):
                continue
            if len(param) <= 1:
                continue
            if param in extra_skip:
                continue
            if "(" in param or ")" in param:
                continue
            if "[" in param:
                continue  # Bracket shorthand (e.g., x[y,z]_domain%%beg[end])

            # Normalize %% to % for lookup
            normalized = param.replace("%%", "%")
            if not _is_valid_param(normalized, valid_params, sub_params):
                errors.append(f"  {doc_rel} references parameter '{param}' not in REGISTRY")

    return errors


def check_page_refs(repo_root: Path) -> list[str]:
    """Check that @ref targets in docs reference existing page identifiers."""
    doc_dir = repo_root / "docs" / "documentation"
    if not doc_dir.exists():
        return []

    # Collect all @page identifiers
    page_ids = {"citelist"}  # Doxygen built-in
    for md_file in doc_dir.glob("*.md"):
        text = md_file.read_text(encoding="utf-8")
        m = re.search(r"^\s*@page\s+(\w+)", text, flags=re.MULTILINE)
        if m:
            page_ids.add(m.group(1))

    errors = []
    for md_file in sorted(doc_dir.glob("*.md")):
        text = _strip_code_blocks(md_file.read_text(encoding="utf-8"))
        rel = md_file.relative_to(repo_root)
        for match in REF_RE.finditer(text):
            ref_target = match.group(1)
            if ref_target not in page_ids:
                errors.append(f"  {rel} uses @ref {ref_target} but no @page with that ID exists")

    return errors


def main():
    repo_root = Path(__file__).resolve().parents[2]

    all_errors = []
    all_errors.extend(check_docs(repo_root))
    all_errors.extend(check_cite_keys(repo_root))
    all_errors.extend(check_param_refs(repo_root))
    all_errors.extend(check_page_refs(repo_root))

    if all_errors:
        print("Doc reference check failed:")
        for e in all_errors:
            print(e)
        sys.exit(1)


if __name__ == "__main__":
    main()
