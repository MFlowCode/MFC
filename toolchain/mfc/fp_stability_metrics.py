"""Pure metrics, source-resolution, and parsing helpers for the FP-stability suite.

Leaf module: imports only stdlib + MFC_ROOT_DIR. No sibling fp_stability*
imports, so the runners/report/orchestrator modules can all depend on it.
"""

import glob
import math
import os
import re

from .common import MFC_ROOT_DIR

# Mantissa-bit levels for the VPREC sweep (C).
# 52 = full double, 23 = single, 16 = half-ish, 10 = ultra-low.
VPREC_MANTISSA_BITS = [52, 23, 16, 10]

_OUTPUT_DAT = re.compile(r"^(cons|prim)\.\d+\.\d+\.(\d+)\.dat$")


def _autodetect_compare(filenames: list) -> list:
    """Pick the D/ output files to diff for a user-supplied case: the conserved-
    variable files at the latest written time step (falling back to primitive
    files if none are written). Returns [] if the case produced no field output."""
    by_step = {}
    for f in filenames:
        m = _OUTPUT_DAT.match(os.path.basename(f))
        if m:
            by_step.setdefault(int(m.group(2)), {"cons": [], "prim": []})[m.group(1)].append(os.path.basename(f))
    if not by_step:
        return []
    last = by_step[max(by_step)]
    return sorted(last["cons"] or last["prim"])


# Stability pass/fail (stage A) is scale-free: a case must retain at least this
# many significant bits under random rounding (sig_bits = -log2(max_dev/scale)).
# 24 ~= single precision. One global floor replaces per-case absolute thresholds
# (which spanned 6 orders of magnitude purely from field scale + conditioning);
# normalising by the field scale collapses that, so a single number suffices.
MIN_SIG_BITS = 24


def _sig_bits(max_dev: float, ref_scale: float) -> float:
    """Significant bits retained = -log2(max_dev / ref_scale).

    Scale-free: dividing the deviation by the field's peak magnitude removes the
    absolute scale, leaving only the conditioning.  Zero deviation (or zero
    scale) returns 53.0 = full double precision retained.
    """
    if not (max_dev > 0) or not (ref_scale > 0):
        return 53.0
    return -math.log2(max_dev / ref_scale)


# Files to exclude from cancellation / float-max reports (runtime loaders, XALT).
_EXTERNAL_SRCS = ("xalt", "dl-init", "ld-linux", "libc.so", "libm.so")

# Matches the first "at" frame in a Valgrind stack trace: "(file.fpp:LINE)".
_VGFRAME_RE = re.compile(r"\(([^):]+\.(?:fpp|f90|F90|c|cpp))\s*:(\d+)\)")

# Fypp block directives. The duplicating ones (#:for expands to N copies, #:def
# defines a macro instantiated at multiple call sites) collapse many distinct
# generated computations onto a single .fpp source line, so a cancellation site
# inside one cannot be pinned to a unique runtime instance. #:if/#:with/#:mute
# select code but do not duplicate it, so they are tracked for balance but not flagged.
_FYPP_BLOCK_OPEN = re.compile(r"^\s*#:(for|def|block|call|if|with|mute)\b", re.IGNORECASE)
_FYPP_BLOCK_CLOSE = re.compile(r"^\s*#:end(for|def|block|call|if|with|mute)?\b", re.IGNORECASE)
_FYPP_DUPLICATING = ("for", "def", "block", "call")


def _resolve_source(fname: str) -> str:
    """Resolve a (possibly bare) source filename to an existing path, or '' if not
    found.  An absolute existing path is used as-is; otherwise the basename is
    located recursively under src/."""
    if os.path.isabs(fname) and os.path.isfile(fname):
        return fname
    candidates = glob.glob(os.path.join(MFC_ROOT_DIR, "src", "**", os.path.basename(fname)), recursive=True)
    return candidates[0] if candidates else ""


def _read_source_lines(fname: str) -> list:
    """Resolve `fname` and return its lines (with newlines), or [] if unreadable."""
    path = _resolve_source(fname)
    if not path:
        return []
    try:
        with open(path) as fh:
            return fh.readlines()
    except OSError:
        return []


def _macro_context_in_lines(lines: list, lineno: int) -> str:
    """Return the innermost code-duplicating fypp block ('#:for'/'#:def'/...) that
    encloses `lineno` (1-based) in `lines`, or None if none does.

    Used to flag cancellation sites whose .fpp line is shared across multiple
    expanded instances (a #:for body, a #:def macro used in many places), where
    line-level attribution cannot identify which instance is responsible.
    """
    stack = []
    for raw in lines[: max(0, lineno - 1)]:
        mo = _FYPP_BLOCK_OPEN.match(raw)
        if mo:
            stack.append(mo.group(1).lower())
            continue
        if _FYPP_BLOCK_CLOSE.match(raw) and stack:
            stack.pop()
    for kw in reversed(stack):
        if kw in _FYPP_DUPLICATING:
            return f"#:{kw}"
    return None


def _macro_context(fname: str, lineno: int) -> str:
    """File-backed wrapper around _macro_context_in_lines; '' path safe."""
    lines = _read_source_lines(fname)
    if not lines:
        return None
    return _macro_context_in_lines(lines, lineno)


def _dat_column(path: str):
    """Load column 1 (the field value) from an MFC .dat file, robust to a
    single-row file (np.loadtxt returns 1-D then, which [:, 1] would crash on)."""
    import numpy as np

    return np.atleast_2d(np.loadtxt(path))[:, 1]


def _max_diff_np(ref_dir: str, run_dir: str, compare_files: list) -> float:
    import numpy as np

    total = 0.0
    for fname in compare_files:
        ref_p, run_p = os.path.join(ref_dir, fname), os.path.join(run_dir, fname)
        if not os.path.exists(ref_p) or not os.path.exists(run_p):
            return float("inf")
        total = max(total, float(np.max(np.abs(_dat_column(ref_p) - _dat_column(run_p)))))
    return total


def _max_abs_np(ref_dir: str, compare_files: list) -> float:
    """Return the maximum absolute value across all reference output files."""
    import numpy as np

    total = 0.0
    for fname in compare_files:
        ref_p = os.path.join(ref_dir, fname)
        if not os.path.exists(ref_p):
            continue
        total = max(total, float(np.max(np.abs(_dat_column(ref_p)))))
    return total


def _parse_cancel_gen(gen_path: str) -> list:
    """Parse cc-gen-file TSV (file\\tline\\tsymbol) -> sorted unique [(fname, line)] for MFC sources."""
    if not os.path.isfile(gen_path):
        return []
    locs = []
    seen = set()
    with open(gen_path) as fh:
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            fname = parts[0].strip()
            if any(ext in fname for ext in _EXTERNAL_SRCS):
                continue
            if not fname.endswith((".fpp", ".f90", ".F90", ".c", ".cpp")):
                continue
            try:
                lineno = int(parts[1].strip())
            except ValueError:
                continue
            key = (fname, lineno)
            if key not in seen:
                seen.add(key)
                locs.append(key)
    return locs


def _parse_vg_error_locs(log_path: str, error_keyword: str) -> list:
    """Extract first MFC-source frame from each Valgrind error matching error_keyword."""
    if not os.path.isfile(log_path):
        return []
    locs = []
    seen = set()
    in_error = False
    with open(log_path) as fh:
        for raw in fh:
            line = re.sub(r"^==\d+== ?", "", raw)
            if error_keyword in line:
                in_error = True
                continue
            if in_error:
                if "   at " in line or "   by " in line:
                    m = _VGFRAME_RE.search(line)
                    if m:
                        fname = m.group(1)
                        if any(ext in fname for ext in _EXTERNAL_SRCS):
                            continue
                        lineno = int(m.group(2))
                        key = (fname, lineno)
                        if key not in seen:
                            seen.add(key)
                            locs.append(key)
                        in_error = False
                elif line.strip() == "":
                    in_error = False
    return locs


# Verrou exposes no per-site bit-count, but --cc-threshold-double is a severity
# filter: a site is reported only if it lost >= the threshold bits. Sweeping these
# levels and taking the highest each site survives gives a per-site "bits lost"
# severity (a lower bound - no false positives). 48 is near the full 53-bit
# double mantissa (the top of the sweep), not the mantissa width itself.
CANCEL_BIT_LEVELS = [10, 20, 30, 40, 48]


def _cancellation_severity(level_sites: list) -> dict:
    """Given [(threshold, [sites])], return {site: highest threshold it survives}
    = the per-site bits-lost severity (a lower bound)."""
    sev = {}
    for level, sites in level_sites:
        for site in sites:
            if level > sev.get(site, 0):
                sev[site] = level
    return sev


def _digits_left(bits_lost: float) -> float:
    """Approximate trustworthy decimal digits remaining after losing `bits_lost`
    bits of a double's 53-bit mantissa (~15.95 digits full)."""
    return max(0.0, (53 - bits_lost) / math.log2(10))
