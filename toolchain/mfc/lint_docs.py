"""Check that file paths, cite keys, parameters, and @ref targets in docs still exist."""

import re
import sys
from pathlib import Path

# Docs to scan for file path references (all hand-written docs with code refs)
DOCS = [
    "docs/documentation/contributing.md",
    "docs/documentation/gpuParallelization.md",
    "docs/documentation/running.md",
    "docs/documentation/case.md",
    "docs/documentation/equations.md",
    "docs/documentation/testing.md",
    "docs/documentation/getting-started.md",
    "docs/documentation/docker.md",
    "docs/documentation/troubleshooting.md",
    "docs/documentation/visualization.md",
    "docs/documentation/expectedPerformance.md",
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
                errors.append(
                    f"  {doc} references '{path_str}' but it does not exist."
                    " Fix: update the path or remove the reference"
                )
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
                errors.append(
                    f"  {rel} uses \\cite {key} but no bib entry found."
                    " Fix: add entry to docs/references.bib or fix the key"
                )

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

    # Compound params (with %): validate both family prefix and attribute
    if "%" in param:
        sub = param.split("%")[-1]
        family_ok = any(p.startswith(base + "%") for p in valid_params)
        return family_ok and sub in sub_params

    # Simple params: family prefix check at structural boundaries
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
                errors.append(
                    f"  {doc_rel} references parameter '{param}' not in REGISTRY."
                    " Fix: check spelling or add to definitions.py"
                )

    return errors


def check_math_syntax(repo_root: Path) -> list[str]:
    """Check that docs use Doxygen math syntax (\\f$, \\f[) not raw $$/$."""
    doc_dir = repo_root / "docs" / "documentation"
    if not doc_dir.exists():
        return []

    errors = []
    for md_file in sorted(doc_dir.glob("*.md")):
        text = md_file.read_text(encoding="utf-8")
        rel = md_file.relative_to(repo_root)
        in_code = False

        for i, line in enumerate(text.split("\n"), 1):
            if line.strip().startswith("```"):
                in_code = not in_code
                continue
            if in_code:
                continue

            # Strip inline code and Doxygen math before checking
            cleaned = re.sub(r"`[^`\n]+`", "", line)
            cleaned = re.sub(r"\\f\$.*?\\f\$", "", cleaned)
            cleaned = re.sub(r"\\f\[.*?\\f\]", "", cleaned)

            if "$$" in cleaned:
                errors.append(
                    f"  {rel}:{i} uses $$...$$ display math."
                    " Fix: replace $$ with \\f[ and \\f]"
                )
                continue

            for m in re.finditer(r"\$([^$\n]+?)\$", cleaned):
                if re.search(r"\\[a-zA-Z]", m.group(1)):
                    errors.append(
                        f"  {rel}:{i} uses $...$ with LaTeX commands."
                        " Fix: replace $ delimiters with \\f$ and \\f$"
                    )
                    break  # one error per line

    return errors


def _gitignored_docs(repo_root: Path) -> set[str]:
    """Return set of gitignored doc file basenames."""
    import subprocess  # pylint: disable=import-outside-toplevel

    doc_dir = repo_root / "docs" / "documentation"
    try:
        result = subprocess.run(
            ["git", "ls-files", "--ignored", "--exclude-standard", "-o"],
            capture_output=True, text=True, cwd=repo_root, check=False,
        )
        return {
            Path(f).name for f in result.stdout.splitlines()
            if f.startswith(str(doc_dir.relative_to(repo_root)))
        }
    except FileNotFoundError:
        return set()


def check_section_anchors(repo_root: Path) -> list[str]:
    """Check that markdown ](#id) links have matching {#id} definitions."""
    doc_dir = repo_root / "docs" / "documentation"
    if not doc_dir.exists():
        return []

    ignored = _gitignored_docs(repo_root)

    errors = []
    for md_file in sorted(doc_dir.glob("*.md")):
        if md_file.name in ignored:
            continue
        text = md_file.read_text(encoding="utf-8")
        rel = md_file.relative_to(repo_root)

        # Collect all {#id} anchors (outside code blocks)
        anchors = set()
        in_code = False
        for line in text.split("\n"):
            if line.strip().startswith("```"):
                in_code = not in_code
                continue
            if not in_code:
                anchors.update(re.findall(r"\{#([\w-]+)\}", line))

        # Check all ](#id) links
        in_code = False
        for i, line in enumerate(text.split("\n"), 1):
            if line.strip().startswith("```"):
                in_code = not in_code
                continue
            if in_code:
                continue
            for m in re.finditer(r"\]\(#([\w-]+)\)", line):
                if m.group(1) not in anchors:
                    errors.append(
                        f"  {rel}:{i} links to #{m.group(1)}"
                        f" but no {{#{m.group(1)}}} anchor exists."
                        f" Fix: add {{#{m.group(1)}}} to the target"
                        " section header"
                    )

    return errors


def check_doxygen_percent(repo_root: Path) -> list[str]:
    """Check that Fortran % accessors inside backtick code spans use %% (Doxygen escape).

    Doxygen treats %<word> as "suppress auto-link" and silently eats the %
    character, even inside backtick code spans.  Writing %% produces a
    literal %.  Only flag % followed by [a-zA-Z_] (the pattern Doxygen
    consumes); %[, %(, etc. are safe.
    """
    doc_dir = repo_root / "docs" / "documentation"
    if not doc_dir.exists():
        return []

    ignored = _gitignored_docs(repo_root)
    code_span_re = re.compile(r"``([^`\n]+)``|`([^`\n]+)`")
    bad_pct_re = re.compile(r"(?<!%)%(?=[a-zA-Z_])")

    errors = []
    for md_file in sorted(doc_dir.glob("*.md")):
        if md_file.name in ignored:
            continue
        text = md_file.read_text(encoding="utf-8")
        rel = md_file.relative_to(repo_root)
        in_code = False
        for i, line in enumerate(text.split("\n"), 1):
            if line.strip().startswith("```"):
                in_code = not in_code
                continue
            if in_code:
                continue
            for m in code_span_re.finditer(line):
                span = m.group(1) or m.group(2)
                if bad_pct_re.search(span):
                    fixed = bad_pct_re.sub("%%", span)
                    errors.append(
                        f"  {rel}:{i} Doxygen will eat the % in `{span}`."
                        f" Fix: `{fixed}`"
                    )

    return errors


def check_page_refs(repo_root: Path) -> list[str]:
    """Check that @ref targets in docs reference existing page identifiers."""
    doc_dir = repo_root / "docs" / "documentation"
    if not doc_dir.exists():
        return []

    # Collect all @page identifiers
    # Include Doxygen built-ins and auto-generated pages (created by ./mfc.sh generate)
    page_ids = {"citelist", "parameters", "case_constraints", "physics_constraints", "examples", "cli-reference"}
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
                errors.append(
                    f"  {rel} uses @ref {ref_target} but no @page with that ID exists."
                    " Fix: check the page ID or add @page declaration"
                )

    return errors


def check_physics_docs_coverage(repo_root: Path) -> list[str]:
    """Check that all check methods with enforcement calls have PHYSICS_DOCS entries."""
    toolchain_dir = str(repo_root / "toolchain")
    if toolchain_dir not in sys.path:
        sys.path.insert(0, toolchain_dir)
    try:
        from mfc.case_validator import PHYSICS_DOCS  # pylint: disable=import-outside-toplevel
        from mfc.params.ast_analyzer import analyze_case_validator  # pylint: disable=import-outside-toplevel
    except ImportError:
        return []

    # Methods without PHYSICS_DOCS entries. Add a PHYSICS_DOCS entry (with math,
    # references, and explanation) to case_validator.py to remove from this set.
    skip = {
        # Structural/mechanical checks (no physics meaning)
        "check_parameter_types",     # type validation
        "check_output_format",       # output format selection
        "check_restart",             # restart file logistics
        "check_parallel_io_pre_process",  # parallel I/O settings
        "check_misc_pre_process",    # miscellaneous pre-process flags
        "check_bc_patches",          # boundary patch geometry
        "check_grid_stretching",     # grid stretching parameters
        "check_qbmm_pre_process",    # QBMM pre-process settings
        "check_probe_integral_output",  # probe/integral output settings
        "check_finite_difference",   # fd_order value validation
        "check_flux_limiter",        # output dimension requirements
        "check_liutex_post",         # output dimension requirements
        "check_momentum_post",       # output dimension requirements
        "check_velocity_post",       # output dimension requirements
        "check_surface_tension_post",  # output feature dependency
        "check_no_flow_variables",   # output variable selection
        "check_partial_domain",      # output format settings
        "check_perturb_density",     # parameter pairing validation
        "check_qm",                  # output dimension requirements
        "check_chemistry",           # runtime Cantera validation only
        # Awaiting proper physics documentation (math, references, explanation)
        "check_adaptive_time_stepping",
        "check_adv_n",
        "check_body_forces",
        "check_continuum_damage",
        "check_grcbc",
        "check_hyperelasticity",
        "check_ibm",
        "check_igr_simulation",
        "check_mhd_simulation",
        "check_model_eqns_simulation",
        "check_muscl_simulation",
        "check_partial_density",
        "check_qbmm_and_polydisperse",
        "check_riemann_solver",
        "check_schlieren",
        "check_volume_fraction",
        "check_weno_simulation",
    }

    validator_path = repo_root / "toolchain" / "mfc" / "case_validator.py"
    analysis = analyze_case_validator(validator_path)
    rules = analysis["rules"]

    # Find methods that have at least one prohibit/warn call
    methods_with_rules = {r.method for r in rules}

    errors = []
    for method in sorted(methods_with_rules):
        if method in PHYSICS_DOCS:
            continue
        if method in skip:
            continue
        errors.append(
            f"  {method} has validation rules but no PHYSICS_DOCS entry."
            " Fix: add entry to PHYSICS_DOCS in case_validator.py"
            " or add to skip set in lint_docs.py"
        )

    return errors


# Important Python identifiers in contributing.md mapped to files where they must exist.
# If an identifier is renamed or removed, this check catches the stale doc reference.
_CONTRIBUTING_IDENTIFIERS = {
    "PHYSICS_DOCS": "toolchain/mfc/case_validator.py",
    "CONSTRAINTS": "toolchain/mfc/params/definitions.py",
    "DEPENDENCIES": "toolchain/mfc/params/definitions.py",
    "REGISTRY": "toolchain/mfc/params/__init__.py",
}


def check_identifier_refs(repo_root: Path) -> list[str]:
    """Check that important identifiers referenced in contributing.md still exist."""
    doc_path = repo_root / "docs" / "documentation" / "contributing.md"
    if not doc_path.exists():
        return []

    text = _strip_code_blocks(doc_path.read_text(encoding="utf-8"))
    errors = []

    for identifier, source_file in _CONTRIBUTING_IDENTIFIERS.items():
        # Check identifier is actually referenced in the doc
        if f"`{identifier}" not in text:
            continue
        source_path = repo_root / source_file
        if not source_path.exists():
            errors.append(
                f"  contributing.md references `{identifier}` in {source_file}"
                f" but {source_file} does not exist"
            )
            continue
        source_text = source_path.read_text(encoding="utf-8")
        if identifier not in source_text:
            errors.append(
                f"  contributing.md references `{identifier}` but it was not"
                f" found in {source_file}. Fix: update the docs or the identifier"
            )

    return errors


# Match ./mfc.sh <command> patterns (the subcommand name)
_CLI_CMD_RE = re.compile(r"\./mfc\.sh\s+([a-z][a-z_-]*)")


def check_cli_refs(repo_root: Path) -> list[str]:
    """Check that ./mfc.sh commands referenced in docs exist in the CLI schema."""
    toolchain_dir = str(repo_root / "toolchain")
    if toolchain_dir not in sys.path:
        sys.path.insert(0, toolchain_dir)
    try:
        from mfc.cli.commands import MFC_CLI_SCHEMA  # pylint: disable=import-outside-toplevel
    except ImportError:
        return []

    valid_commands = {cmd.name for cmd in MFC_CLI_SCHEMA.commands}
    # Also accept "init" (shell function) and "load" (shell function)
    valid_commands.update({"init", "load"})

    errors = []
    doc_path = repo_root / "docs" / "documentation" / "running.md"
    if not doc_path.exists():
        return []

    text = _strip_code_blocks(doc_path.read_text(encoding="utf-8"))
    seen = set()
    for match in _CLI_CMD_RE.finditer(text):
        cmd = match.group(1)
        if cmd in seen or cmd in valid_commands:
            seen.add(cmd)
            continue
        seen.add(cmd)
        errors.append(
            f"  running.md references './mfc.sh {cmd}' but '{cmd}'"
            " is not a known CLI command."
            " Fix: update the command name or remove the reference"
        )

    return errors


def main():
    repo_root = Path(__file__).resolve().parents[2]

    all_errors = []
    all_errors.extend(check_docs(repo_root))
    all_errors.extend(check_cite_keys(repo_root))
    all_errors.extend(check_param_refs(repo_root))
    all_errors.extend(check_page_refs(repo_root))
    all_errors.extend(check_math_syntax(repo_root))
    all_errors.extend(check_doxygen_percent(repo_root))
    all_errors.extend(check_section_anchors(repo_root))
    all_errors.extend(check_physics_docs_coverage(repo_root))
    all_errors.extend(check_identifier_refs(repo_root))
    all_errors.extend(check_cli_refs(repo_root))

    if all_errors:
        print("Doc reference check failed:")
        for e in all_errors:
            print(e)
        sys.exit(1)


if __name__ == "__main__":
    main()
