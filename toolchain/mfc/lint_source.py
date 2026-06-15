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
    """Flag @:ALLOCATE'd array names with no matching @:DEALLOCATE anywhere in src/.

    Every @:ALLOCATE must have a matching @:DEALLOCATE so that GPU device
    memory is properly released alongside host memory.
    """
    errors: list[str] = []
    src_dir = repo_root / SRC_DIR

    # Known program-lifetime allocations with no finalize subroutine
    PROGRAM_LIFETIME_EXEMPTIONS = {
        # m_helper
        "weight",
        "pb0",
        "Re_trans_T",
        # m_model
        "gpu_ntrs",
        "gpu_trs_v",
        "gpu_trs_n",
        "gpu_boundary_edge_count",
        "gpu_total_vertices",
        "gpu_boundary_v",
        # m_variables_conversion
        "gs_min",
        "pi_infs",
        "ps_inf",
        "cvs",
        "qvs",
        "qvps",
        "Gs_vc",
        # m_start_up post_process
        "data_in",
        "data_out",
        "data_cmplx",
        "data_cmplx_y",
        "data_cmplx_z",
        "En_real",
        "En",
        # m_acoustic_src
        "loc_acoustic",
        "mass_src",
        "mom_src",
        "E_src",
        "source_spatials_num_points",
        "source_spatials",
        # m_bubbles_EE
        "bub_adv_src",
        "bub_r_src",
        "bub_v_src",
        "bub_p_src",
        "bub_m_src",
        "divu",
        "ms",
        "ps",
        "rs",
        "vs",
        # m_qbmm
        "bubmoms",
        "momrhs",
        # m_cbc
        "F_src_rsx_vf",
        "flux_src_rsx_vf_l",
        "F_src_rsy_vf",
        "flux_src_rsy_vf_l",
        "F_src_rsz_vf",
        "flux_src_rsz_vf_l",
        "pres_in",
        "Del_in",
        "alpha_rho_in",
        # m_data_output
        "Rc_sf",
        # m_fftw
        "data_cmplx_gpu",
        "data_fltr_cmplx_gpu",
        # m_global_parameters
        "qbmm_idx",
        "x_cc",
        "dx",
        "y_cc",
        "dy",
        "z_cc",
        "dz",
        # m_hypoelastic
        "dw_dx_hypo",
        # m_igr
        "coeff_R",
        # m_rhs
        "qR_rsx_vf",
        "dqR_rsx_vf",
        # m_surface_tension
        "gR_x",
        # m_time_steppers
        "q_prim_ts2",
        # m_weno
        "poly_coef_cbR_x",
        "d_cbR_x",
        "poly_coef_cbR_y",
        "d_cbR_y",
        "poly_coef_cbR_z",
        "d_cbR_z",
        "divu%sf",
        "qbmm_idx%ps",
    }

    allocate_re = re.compile(r"@:ALLOCATE\((\w[\w%]*)")
    deallocate_re = re.compile(r"@:DEALLOCATE\((\w[\w%]*)")

    # Collect all deallocated names across all files
    deallocated: set[str] = set()
    for src in _fortran_fpp_files(src_dir):
        text = src.read_text(encoding="utf-8")
        for match in deallocate_re.finditer(text):
            deallocated.add(match.group(1))

    # Check each allocation against the global deallocated set
    for src in _fortran_fpp_files(src_dir):
        lines = src.read_text(encoding="utf-8").splitlines()
        rel = src.relative_to(repo_root)
        for i, line in enumerate(lines):
            stripped = line.strip()
            if _is_comment_or_blank(stripped):
                continue
            match = allocate_re.search(stripped)
            if match:
                name = match.group(1)
                if name not in deallocated and name not in PROGRAM_LIFETIME_EXEMPTIONS:
                    errors.append(
                        f"  {rel}:{i + 1} @:ALLOCATE({name}...) has no matching "
                        f"@:DEALLOCATE anywhere in src/. Fix: add @:DEALLOCATE in "
                        f"the module's s_finalize_* subroutine or annotate as program-lifetime"
                    )

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

    if all_errors:
        print("Source lint failed:")
        for e in all_errors:
            print(e)
        sys.exit(1)


if __name__ == "__main__":
    main()
