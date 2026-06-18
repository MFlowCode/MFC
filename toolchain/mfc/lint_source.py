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
    precision_re = re.compile(
        r"\b(?:double_precision|double\s+precision|dsqrt|dexp|dlog|dble|dabs|" r"dprod|dmin|dmax|dfloat|dreal|dcos|dsin|dtan|dsign|dtanh|dsinh|dcosh)\b|" r"\breal\s*\(\s*[48]\s*\)|" r"[0-9]d0",
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


def check_allocate_deallocate_pairing(repo_root: Path) -> list[str]:
    """Flag @:ALLOCATE'd names in s_initialize_* with no matching @:DEALLOCATE in s_finalize_*.

    Only checks modules that have both an s_initialize_* and s_finalize_* subroutine.
    Files without a finalize subroutine are skipped — those allocations are
    program-lifetime by convention.
    """
    errors: list[str] = []
    src_dir = repo_root / SRC_DIR

    allocate_re = re.compile(r"@:ALLOCATE\((\w[\w%]*)")
    deallocate_re = re.compile(r"@:DEALLOCATE\((\w[\w%]*)")
    init_re = re.compile(r"^\s*(?:impure\s+)?subroutine\s+s_initialize_\w+", re.IGNORECASE)
    final_re = re.compile(r"^\s*(?:impure\s+)?subroutine\s+s_finalize_\w+", re.IGNORECASE)
    end_sub_re = re.compile(r"^\s*end\s+subroutine\b", re.IGNORECASE)

    # Known pre-existing missing deallocations
    KNOWN_MISSING = {
        "gs_min",
        "pi_infs",
        "ps_inf",
        "cvs",
        "qvs",
        "qvps",
        "Gs_vc",
        "data_in",
        "data_out",
        "data_cmplx",
        "data_cmplx_y",
        "data_cmplx_z",
        "En_real",
        "En",
        "F_src_rsx_vf",
        "flux_src_rsx_vf_l",
        "F_src_rsy_vf",
        "flux_src_rsy_vf_l",
        "F_src_rsz_vf",
        "flux_src_rsz_vf_l",
        "pres_in",
        "Del_in",
        "alpha_rho_in",
        "Rc_sf",
        "data_cmplx_gpu",
        "data_fltr_cmplx_gpu",
        "qbmm_idx",
        "qbmm_idx%ps",
        "x_cc",
        "dx",
        "y_cc",
        "dy",
        "z_cc",
        "dz",
        "dw_dx_hypo",
        "coeff_R",
        "qR_rsx_vf",
        "dqR_rsx_vf",
        "gR_x",
        "q_prim_ts2",
        "poly_coef_cbR_x",
        "d_cbR_x",
        "poly_coef_cbR_y",
        "d_cbR_y",
        "poly_coef_cbR_z",
        "d_cbR_z",
    }

    for src in _fortran_fpp_files(src_dir):
        lines = src.read_text(encoding="utf-8").splitlines()
        rel = src.relative_to(repo_root)

        # Extract allocations from s_initialize_* subroutines
        allocated: dict[str, int] = {}  # name -> line number
        deallocated: set[str] = set()
        in_init = False
        in_final = False
        has_final = False

        for i, line in enumerate(lines):
            stripped = line.strip()
            if init_re.match(stripped):
                in_init = True
                in_final = False
            elif final_re.match(stripped):
                in_final = True
                in_init = False
                has_final = True
            elif end_sub_re.match(stripped):
                in_init = False
                in_final = False

            if _is_comment_or_blank(stripped):
                continue

            if in_init:
                m = allocate_re.search(stripped)
                if m:
                    name = m.group(1)
                    if name not in allocated:
                        allocated[name] = i + 1

            if in_final:
                m = deallocate_re.search(stripped)
                if m:
                    deallocated.add(m.group(1))

        if not has_final:
            continue

        for name, lineno in allocated.items():
            if name not in deallocated and name not in KNOWN_MISSING:
                errors.append(f"  {rel}:{lineno} @:ALLOCATE({name}...) in s_initialize_* has no matching " f"@:DEALLOCATE in s_finalize_*. Fix: add @:DEALLOCATE in the finalize subroutine")

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


def main():
    repo_root = Path(__file__).resolve().parents[2]

    all_errors: list[str] = []
    all_errors.extend(check_raw_directives(repo_root))
    all_errors.extend(check_double_precision(repo_root))
    all_errors.extend(check_junk_code(repo_root))
    all_errors.extend(check_false_integers(repo_root))
    all_errors.extend(check_junk_comments(repo_root))
    all_errors.extend(check_fypp_list_duplicates(repo_root))
    all_errors.extend(check_duplicate_lines(repo_root))
    all_errors.extend(check_hardcoded_byte_size(repo_root))
    all_errors.extend(check_allocate_deallocate_pairing(repo_root))
    all_errors.extend(check_manual_registry_bcasts(repo_root))

    if all_errors:
        print("Source lint failed:")
        for e in all_errors:
            print(e)
        sys.exit(1)


if __name__ == "__main__":
    main()
