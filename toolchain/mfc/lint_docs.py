"""Check that file paths, cite keys, and parameters in docs still exist."""

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


def check_param_refs(repo_root: Path) -> list[str]:
    """Check that parameter names in equations.md exist in the MFC registry."""
    eq_path = repo_root / "docs" / "documentation" / "equations.md"
    if not eq_path.exists():
        return []

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
    # Build set of sub-parameter suffixes (the part after %)
    sub_params = {p.split("%")[-1] for p in valid_params if "%" in p}
    text = eq_path.read_text(encoding="utf-8")
    errors = []
    seen = set()

    for match in PARAM_RE.finditer(text):
        param = match.group(1)
        if param in seen:
            continue
        seen.add(param)

        # Skip non-parameter identifiers
        if PARAM_SKIP.search(param):
            continue
        # Skip single-character names (too ambiguous: m, n, p are grid dims)
        if len(param) <= 1:
            continue
        # Skip names containing % (struct members like x_domain%beg are valid
        # but the base name before % is what matters)
        base = param.split("%")[0] if "%" in param else param
        # Check for indexed parameters: strip trailing _N (e.g., patch_icpp(1)%alpha(1))
        # In docs these appear as e.g. `patch_icpp(i)%vel(j)` — skip indexed refs
        if "(" in param or ")" in param:
            continue

        if base not in valid_params and base not in sub_params:
            # Check if it's a known parameter family prefix
            if not any(p.startswith(base) for p in valid_params):
                errors.append(f"  equations.md references parameter '{param}' not in REGISTRY")

    return errors


def main():
    repo_root = Path(__file__).resolve().parents[2]

    all_errors = []
    all_errors.extend(check_docs(repo_root))
    all_errors.extend(check_cite_keys(repo_root))
    all_errors.extend(check_param_refs(repo_root))

    if all_errors:
        print("Doc reference check failed:")
        for e in all_errors:
            print(e)
        sys.exit(1)


if __name__ == "__main__":
    main()
