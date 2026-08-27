"""Static analysis for Fortran/Fypp source code.

Checks for patterns that indicate copy-paste bugs, non-standard constructs,
and hardcoded assumptions that break under different build configurations.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

# Source directory to scan (relative to repo root)
SRC_DIR = "src"

# Minimum stripped line length to consider for duplicate detection.
# Lines shorter than this (e.g. "end if", "end do") are ignored.
MIN_DUP_LINE_LEN = 40

# Files excluded from the raw OpenACC/OpenMP directive check (these define the macros)
RAW_DIRECTIVE_EXCLUDE = {
    "parallel_macros.fpp",
    "acc_macros.fpp",
    "omp_macros.fpp",
    "shared_parallel_macros.fpp",
    "syscheck.fpp",
}

# Files/dirs excluded from the double-precision intrinsic check
PRECISION_EXCLUDE_DIRS = {"syscheck"}
PRECISION_EXCLUDE_PATTERNS = {"nvtx", "precision_select"}

# MPI proxy source directory -> params-registry target key
MPI_PROXY_TARGETS = {"pre_process": "pre", "simulation": "sim", "post_process": "post"}

# Checker subroutines whose every @:PROHIBIT depends on state the Python
# validator cannot see: the active compiler, the MPI decomposition, or a value
# Cantera fills in at runtime. Constraints between input parameters belong in
# toolchain/mfc/case_validator.py instead -- see check_checker_input_constraints.
#
# Subroutines that mix runtime and input-only checks are deliberately NOT listed.
# s_check_inputs_weno and s_check_inputs_muscl are the case in point: their
# grid-extent checks need per-rank m/n/p, but the muscl_order/int_comp check that
# used to sit alongside them was pure input, and allowlisting the subroutine would
# have let its replacement back in unnoticed. Those lines carry an explicit
# RUNTIME_CHECK_MARKER instead, so anything new added there is still flagged.
RUNTIME_CHECKER_SUBROUTINES = {
    # Compiler conditionals (#ifdef / #if guarded).
    "s_check_amd",
    "s_check_inputs_compilers",
    "s_check_inputs_nvidia_uvm",
    # MPI decomposition: n_global, num_procs_y/z.
    "s_check_total_cells",
    "s_check_inputs_fft",
    # num_species is populated by Cantera at runtime.
    "s_check_inputs_ib_injection",
}

# Opt out of check_checker_input_constraints for a single @:PROHIBIT.
RUNTIME_CHECK_MARKER = "lint: runtime-check"

# Fortran subroutine declaration, allowing any order of the prefixes MFC uses
# (impure/pure/elemental/recursive/module) before the `subroutine` keyword.
_SUBROUTINE_DECL = re.compile(r"(?:(?:impure|pure|elemental|recursive|module|non_recursive)\s+)*subroutine\s+(\w+)")


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

    list_content = full_line[bracket_start + 1 : bracket_end]
    list_content = list_content.replace("&", "")

    # Extract single- or double-quoted entries
    entries = re.findall(r"['\"]([^'\"]*)['\"]", list_content)

    seen: dict[str, int] = {}
    for pos, entry in enumerate(entries, 1):
        if entry in seen:
            errors.append(f"  {rel}:{start_line} Fypp list has duplicate entry '{entry}' (positions {seen[entry]} and {pos}). Fix: one copy is likely a typo for a different variable")
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
            if stripped == prev_stripped and len(stripped) >= MIN_DUP_LINE_LEN and not _is_comment_or_blank(stripped):
                display = stripped[:72]
                if len(stripped) > 72:
                    display += "..."
                errors.append(f"  {rel}:{i + 1} identical to previous line: '{display}'. Fix: check for accidental copy-paste")
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
                errors.append(f"  {rel}:{i + 1} hard-codes 8-byte real size. Fix: use 'storage_size(0._stp)/8' instead of '8._wp'")

    return errors


def check_raw_directives(repo_root: Path) -> list[str]:
    """Flag raw OpenACC/OpenMP directives outside of macro definition files.

    All GPU directives must use the GPU_* Fypp macros from parallel_macros.fpp.
    Note: directives like !$acc start with '!' so they look like Fortran comments.
    We must NOT use _is_comment_or_blank here, and must search the full line.
    """
    errors: list[str] = []
    src_dir = repo_root / SRC_DIR
    directive_re = re.compile(r"!\$acc|!\$omp", re.IGNORECASE)

    for src in _fortran_fpp_files(src_dir):
        if src.name in RAW_DIRECTIVE_EXCLUDE:
            continue
        lines = src.read_text(encoding="utf-8").splitlines()
        rel = src.relative_to(repo_root)

        for i, line in enumerate(lines):
            stripped = line.strip()
            if not stripped or stripped.startswith("#:"):
                continue
            if directive_re.search(stripped):
                errors.append(f"  {rel}:{i + 1} raw OpenACC/OpenMP directive. Fix: use GPU_* Fypp macros instead")

    return errors


def check_double_precision(repo_root: Path) -> list[str]:
    """Flag double-precision-specific intrinsics and type declarations.

    MFC uses generic intrinsics (sqrt, exp, etc.) and kind parameters (wp, stp)
    so that precision can be changed at build time.
    """
    errors: list[str] = []
    src_dir = repo_root / SRC_DIR
    # The d-literal alternative catches double-precision literals with a full
    # mantissa and a signed or multi-digit exponent (5.0d-11, 101325.d0, .5d0,
    # 1013.25d3), not just '[0-9]d0'. The boundaries keep it out of identifiers
    # like cart2d12_coords and out of 'D' edit descriptors like (1D12.4).
    precision_re = re.compile(
        r"\b(?:double_precision|double\s+precision|dsqrt|dexp|dlog|dble|dabs|"
        r"dprod|dmin|dmax|dfloat|dreal|dcos|dsin|dtan|dsign|dtanh|dsinh|dcosh)\b|"
        r"\breal\s*\(\s*[48]\s*\)|"
        r"(?<![A-Za-z0-9_.])(?:[0-9]+\.?[0-9]*|\.[0-9]+)[dD][-+]?[0-9]+(?![A-Za-z0-9_.])",
        re.IGNORECASE,
    )

    for src in _fortran_fpp_files(src_dir):
        if any(p in src.name for p in PRECISION_EXCLUDE_PATTERNS):
            continue
        if any(d in src.parts for d in PRECISION_EXCLUDE_DIRS):
            continue
        lines = src.read_text(encoding="utf-8").splitlines()
        rel = src.relative_to(repo_root)

        for i, line in enumerate(lines):
            stripped = line.strip()
            if _is_comment_or_blank(stripped):
                continue
            code = stripped.split("!")[0]
            match = precision_re.search(code)
            if match:
                errors.append(f"  {rel}:{i + 1} double-precision intrinsic '{match.group()}'. Fix: use generic intrinsics and wp/stp kind parameters")

    return errors


def check_junk_code(repo_root: Path) -> list[str]:
    """Flag junk patterns (..., ---, ===) in Fortran source, including comments.

    Separator comments like ``! ========`` are also forbidden.
    """
    errors: list[str] = []
    src_dir = repo_root / SRC_DIR
    junk_re = re.compile(r"\.\.\.|---|===")

    for src in _fortran_fpp_files(src_dir):
        lines = src.read_text(encoding="utf-8").splitlines()
        rel = src.relative_to(repo_root)

        for i, line in enumerate(lines):
            stripped = line.strip()
            if not stripped or stripped.startswith("#:"):
                continue
            match = junk_re.search(stripped)
            if match:
                errors.append(f"  {rel}:{i + 1} junk code pattern '{match.group()}'. Fix: remove placeholder/separator text")

    return errors


def check_false_integers(repo_root: Path) -> list[str]:
    """Flag bare integer_wp patterns like ``2_wp`` that should be ``2.0_wp``.

    Fortran kind parameters on integers (e.g. 2_wp) produce integers, not reals.
    Almost always the intent is a real literal like 2.0_wp.
    """
    errors: list[str] = []
    src_dir = repo_root / SRC_DIR
    # Match digits followed by _wp, but NOT preceded by '.', digit, 'e', 'E', or '-'
    # (which would indicate a real literal like 1.0_wp or 1e5_wp)
    false_int_re = re.compile(r"(?<![0-9.eE\-])\b[0-9]+_wp\b")

    for src in _fortran_fpp_files(src_dir):
        lines = src.read_text(encoding="utf-8").splitlines()
        rel = src.relative_to(repo_root)

        for i, line in enumerate(lines):
            stripped = line.strip()
            if _is_comment_or_blank(stripped):
                continue
            code = stripped.split("!")[0]
            match = false_int_re.search(code)
            if match:
                errors.append(f"  {rel}:{i + 1} bare integer with _wp kind '{match.group()}'. Fix: use a real literal (e.g. {match.group().replace('_wp', '.0_wp')})")

    return errors


def check_integer_wp(repo_root: Path) -> list[str]:
    """Flag ``integer(wp)`` and ``integer(kind=wp)`` declarations.

    ``wp`` is a floating-point kind parameter; using it as an integer kind is a
    copy-paste error. Integers take the default kind: plain ``integer``.
    """
    errors: list[str] = []
    src_dir = repo_root / SRC_DIR
    integer_wp_re = re.compile(r"\binteger\s*\(\s*(?:kind\s*=\s*)?wp\s*\)", re.IGNORECASE)

    for src in _fortran_fpp_files(src_dir):
        lines = src.read_text(encoding="utf-8").splitlines()
        rel = src.relative_to(repo_root)

        for i, line in enumerate(lines):
            stripped = line.strip()
            if _is_comment_or_blank(stripped):
                continue
            match = integer_wp_re.search(stripped.split("!")[0])
            if match:
                errors.append(f"  {rel}:{i + 1} '{match.group()}' uses a floating-point kind. Fix: use plain 'integer'")

    return errors


def check_junk_comments(repo_root: Path) -> list[str]:
    """Flag junk separator patterns (===, ----+) in Python and shell scripts.

    Three dashes (---) is valid markdown, but four or more is a separator.
    Checks both comment lines and echo/print statements with separator strings.
    """
    errors: list[str] = []
    junk_re = re.compile(r"===|-{4,}")

    # Python files: check comment lines
    for subdir in ["examples", "benchmarks", "toolchain"]:
        d = repo_root / subdir
        if not d.exists():
            continue
        for py in sorted(d.rglob("*.py")):
            lines = py.read_text(encoding="utf-8").splitlines()
            rel = py.relative_to(repo_root)

            for i, line in enumerate(lines):
                stripped = line.strip()
                if not stripped.startswith("#"):
                    continue
                match = junk_re.search(stripped)
                if match:
                    errors.append(f"  {rel}:{i + 1} junk separator pattern '{match.group()}'. Fix: remove separator comment")

    # Shell files: check comments and overly long echo separators
    long_sep_re = re.compile(r"[=]{21,}|-{21,}")
    for subdir in ["toolchain", ".github"]:
        d = repo_root / subdir
        if not d.exists():
            continue
        for sh in sorted(d.rglob("*.sh")):
            lines = sh.read_text(encoding="utf-8").splitlines()
            rel = sh.relative_to(repo_root)

            for i, line in enumerate(lines):
                stripped = line.strip()
                if stripped.startswith("#"):
                    match = junk_re.search(stripped)
                    if match:
                        errors.append(f"  {rel}:{i + 1} junk separator pattern '{match.group()}'. Fix: remove separator comment")
                elif long_sep_re.search(stripped):
                    errors.append(f"  {rel}:{i + 1} echo separator too long (max 20 chars). Fix: shorten to 20 or fewer")

    return errors


_BCAST_CALL_RE = re.compile(r"\bcall\s+MPI_BCAST\s*\(", re.IGNORECASE)
_FYPP_FOR_RE = re.compile(r"#:\s*for\s+(\w+)\s+in\s+\[")
_FYPP_ENDFOR_RE = re.compile(r"#:\s*endfor\b")
_PLACEHOLDER_RE = re.compile(r"\$\{(\w+)\}\$")
_IDENTIFIER_RE = re.compile(r"[A-Za-z_]\w*$")


def _first_bcast_argument(after_paren: str) -> str:
    """Return the first argument of an MPI_BCAST call, given the text after '('."""
    depth = 0
    for pos, ch in enumerate(after_paren):
        if ch == "(":
            depth += 1
        elif ch == ")" and depth > 0:
            depth -= 1
        elif ch in ",)" and depth == 0:
            return after_paren[:pos]
    return after_paren


def _expand_placeholders(arg: str, loops: dict[str, list[str]]) -> list[str]:
    """Substitute Fypp ``${VAR}$`` placeholders with each active loop entry.

    Placeholders whose loop list has no quoted entries (e.g. numeric lists)
    expand to nothing — those arguments cannot name a namelist scalar root.
    """
    match = _PLACEHOLDER_RE.search(arg)
    if not match:
        return [arg]
    expanded: list[str] = []
    for entry in loops.get(match.group(1), []):
        expanded.extend(_expand_placeholders(arg[: match.start()] + entry + arg[match.end() :], loops))
    return expanded


def _normalize_bcast_root(candidate: str) -> str | None:
    """Reduce a broadcast first argument to its root name, or None if not a top-level scalar.

    Struct members (containing '%') are legitimate manual residue; arguments
    that are not plain identifiers after dropping the index part are skipped.
    """
    candidate = candidate.strip()
    if "%" in candidate or "$" in candidate:
        return None
    root = candidate.split("(", 1)[0].strip()
    return root if _IDENTIFIER_RE.fullmatch(root) else None


def _extract_bcast_roots(lines: list[str]) -> list[tuple[int, str]]:
    """Return (line_number, root_name) for every MPI_BCAST first argument.

    Resolves Fypp ``#:for VAR in ['a', 'b', ...]`` loop lists (including '&'
    continuations and nested loops) so list-driven broadcasts are attributed
    to the quoted variable names.
    """
    roots: list[tuple[int, str]] = []
    loop_stack: list[tuple[str, list[str]]] = []

    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        for_match = _FYPP_FOR_RE.match(stripped)
        if for_match:
            full = stripped
            while full.rstrip().endswith("&") and i + 1 < len(lines):
                i += 1
                full += " " + lines[i].strip()
            entries = re.findall(r"['\"]([^'\"]*)['\"]", full[full.find("[") :])
            loop_stack.append((for_match.group(1), entries))
        elif _FYPP_ENDFOR_RE.match(stripped):
            if loop_stack:
                loop_stack.pop()
        elif not stripped.startswith("!"):
            call = _BCAST_CALL_RE.search(stripped)
            if call:
                first_arg = _first_bcast_argument(stripped[call.end() :])
                for candidate in _expand_placeholders(first_arg, dict(loop_stack)):
                    root = _normalize_bcast_root(candidate)
                    if root:
                        roots.append((i + 1, root))
        i += 1

    return roots


def _registry_bcast_vars(repo_root: Path, target: str) -> set[str]:
    """Names auto-broadcast by generated_bcast.fpp for one target.

    Derived from the generator itself so the lint and the generated code can
    never disagree about which scalars are auto-broadcast: struct roots,
    TYPED_DECLS, FORTRAN_ARRAY_DIMS, and per-target exclusions are already
    removed by _classify_scalar_vars.
    """
    toolchain_dir = str(repo_root / "toolchain")
    if toolchain_dir not in sys.path:
        sys.path.insert(0, toolchain_dir)
    from mfc.params.generators.fortran_gen import _classify_scalar_vars

    return set().union(*_classify_scalar_vars(target))


def check_manual_registry_bcasts(repo_root: Path) -> list[str]:
    """Flag hand-written MPI_BCASTs of registry-bound namelist scalars.

    Those scalars are broadcast by the generated include (generated_bcast.fpp);
    a manual copy in m_mpi_proxy.fpp is a duplicate broadcast, or a new
    parameter bypassing the generator. Manual residue (computed variables and
    struct members) is permitted.
    """
    errors: list[str] = []

    for dirname, target in MPI_PROXY_TARGETS.items():
        proxy = repo_root / SRC_DIR / dirname / "m_mpi_proxy.fpp"
        if not proxy.exists():
            continue
        auto_broadcast = _registry_bcast_vars(repo_root, target)
        rel = proxy.relative_to(repo_root)

        for lineno, root in _extract_bcast_roots(proxy.read_text(encoding="utf-8").splitlines()):
            if root in auto_broadcast:
                errors.append(f"  {rel}:{lineno} manual MPI_BCAST of registry-bound scalar '{root}' — it is auto-broadcast via generated_bcast.fpp; remove the manual copy")

    return errors


def check_checker_input_constraints(repo_root: Path) -> list[str]:
    """Keep input-only constraints out of the Fortran m_checker files.

    Constraints between case-file parameters are enforced in
    toolchain/mfc/case_validator.py, which runs before any binary is invoked.
    Adding them to Fortran as well means the same rule is written twice, in two
    languages, and the copies drift.

    A @:PROHIBIT is allowed only inside a subroutine in
    RUNTIME_CHECKER_SUBROUTINES, or under a "! lint: runtime-check <reason>"
    comment. The marker applies to the next @:PROHIBIT, so blank lines, Fypp
    directives, and further comments may sit between the two.
    """
    errors = []

    for path in sorted((repo_root / SRC_DIR).rglob("m_checker*.fpp")):
        rel = path.relative_to(repo_root)
        subroutine = None
        exempt = False

        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            stripped = line.strip()

            match = _SUBROUTINE_DECL.match(stripped)
            if match:
                subroutine = match.group(1)
            elif stripped.startswith("end subroutine"):
                subroutine = None

            if stripped.startswith("!"):
                exempt = exempt or RUNTIME_CHECK_MARKER in stripped
                continue
            if not stripped or stripped.startswith("#"):
                # Blank lines and Fypp/preprocessor directives do not consume the marker.
                continue

            if "@:PROHIBIT" in stripped:
                if not exempt and subroutine not in RUNTIME_CHECKER_SUBROUTINES:
                    if subroutine:
                        where = f"in {subroutine}"
                        allowlist_hint = f"add '{subroutine}' to RUNTIME_CHECKER_SUBROUTINES in {Path(__file__).name}"
                    else:
                        where = "at module scope"
                        allowlist_hint = f"put it in a subroutine listed in RUNTIME_CHECKER_SUBROUTINES in {Path(__file__).name}"
                    errors.append(
                        f"{rel}:{lineno}: @:PROHIBIT {where} looks like an input-only constraint. "
                        f"Add it to a check_* method in toolchain/mfc/case_validator.py instead. "
                        f"If it genuinely needs runtime or compiler state, {allowlist_hint}, "
                        f"or mark it with '! {RUNTIME_CHECK_MARKER} <reason>'."
                    )
                exempt = False

    return errors


def check_cluster_menu_slugs(repo_root: Path) -> list[str]:
    """Keep the ``./mfc.sh load`` cluster menu in sync with toolchain/modules.

    The menu in toolchain/bootstrap/modules.sh is hand-written so it can be
    grouped and coloured by organisation. That is fine, but it drifts: it used
    to offer Summit, which has no module set, and omitted Phoenix IFX and
    Santis, which do. Compare the advertised slugs against the data file.
    """
    modules = repo_root / "toolchain" / "modules"
    script = repo_root / "toolchain" / "bootstrap" / "modules.sh"
    if not modules.exists() or not script.exists():
        return []

    # Cluster definitions in toolchain/modules are "<slug> <System Name>". The
    # "<slug>-{all,cpu,gpu}[-unload] <modules...>" lines carry the module lists.
    module_list_line = re.compile(r"-(all|cpu|gpu)(-unload)?$")
    defined = set()
    for raw_line in modules.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        slug = line.split()[0]
        if not module_list_line.search(slug):
            defined.add(slug)

    # The menu is delimited by explicit markers rather than by prose, so rewording
    # the prompt or the read cannot silently disable this check.
    text = script.read_text()
    block = re.search(r"# lint: cluster-menu-begin\n(.*?)# lint: cluster-menu-end", text, re.S)
    if block is None:
        return [f"{script.relative_to(repo_root)}: cluster-menu-begin/end markers are missing; check_cluster_menu_slugs cannot verify the menu against toolchain/modules"]
    # Slugs are advertised as "Name (slug)"; require the closing paren to be
    # followed by a separator so shell fragments like ${G} are not picked up.
    advertised = set(re.findall(r"\((\w[\w-]*)\)(?=[\s|\"']|$)", block.group(1), re.M))

    errors = []
    rel = script.relative_to(repo_root)
    for slug in sorted(advertised - defined):
        errors.append(f"{rel}: cluster menu offers '{slug}', which has no entry in toolchain/modules (selecting it loads nothing)")
    for slug in sorted(defined - advertised):
        errors.append(f"{rel}: toolchain/modules defines '{slug}', but the cluster menu does not offer it (users cannot discover it)")
    return errors


def main():
    repo_root = Path(__file__).resolve().parents[2]

    all_errors: list[str] = []
    all_errors.extend(check_raw_directives(repo_root))
    all_errors.extend(check_double_precision(repo_root))
    all_errors.extend(check_junk_code(repo_root))
    all_errors.extend(check_false_integers(repo_root))
    all_errors.extend(check_integer_wp(repo_root))
    all_errors.extend(check_junk_comments(repo_root))
    all_errors.extend(check_fypp_list_duplicates(repo_root))
    all_errors.extend(check_duplicate_lines(repo_root))
    all_errors.extend(check_hardcoded_byte_size(repo_root))
    all_errors.extend(check_manual_registry_bcasts(repo_root))
    all_errors.extend(check_checker_input_constraints(repo_root))
    all_errors.extend(check_cluster_menu_slugs(repo_root))

    if all_errors:
        print("Source lint failed:")
        for e in all_errors:
            print(e)
        sys.exit(1)


if __name__ == "__main__":
    main()
