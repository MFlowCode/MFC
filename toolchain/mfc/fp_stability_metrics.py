"""Pure metrics, source-resolution, and parsing helpers for the FP-stability suite.

Leaf module: imports only stdlib + MFC_ROOT_DIR + cons. No sibling fp_stability*
imports, so the runners/report/orchestrator modules can all depend on it.
"""

import glob
import math
import os
import re

from .common import MFC_ROOT_DIR
from .printer import cons

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

# Fallback absolute threshold for the dd_sym/dd_line oracle when no float-proxy-
# derived threshold is supplied (callers always pass one, so this is only a guard).
_DD_FALLBACK_THRESHOLD = 1e-12


def _sig_bits(max_dev: float, ref_scale: float) -> float:
    """Significant bits retained = -log2(max_dev / ref_scale).

    Scale-free: dividing the deviation by the field's peak magnitude removes the
    absolute scale, leaving only the conditioning.  Zero deviation (or zero
    scale) returns 53.0 = full double precision retained.
    """
    if not (max_dev > 0) or not (ref_scale > 0):
        return 53.0
    return -math.log2(max_dev / ref_scale)


def _stability_pass(max_dev: float, ref_scale: float, floor: float) -> bool:
    """A case passes when it retains at least `floor` significant bits."""
    return _sig_bits(max_dev, ref_scale) >= floor


# Matches "path/file.f90:123" or "path/file.fpp:123-456" in dd_line rddmin_summary.
_LOC_RE = re.compile(r"(\S+\.(?:f90|fpp|c|cpp|h|F90))\s*:(\d+)(?:-(\d+))?", re.IGNORECASE)

# Files to exclude from cancellation / float-max reports (runtime loaders, XALT).
_EXTERNAL_SRCS = ("xalt", "dl-init", "ld-linux", "libc.so", "libm.so")

# Matches the first "at" frame in a Valgrind stack trace: "(file.fpp:LINE)".
_VGFRAME_RE = re.compile(r"\(([^):]+\.(?:fpp|f90|F90|c|cpp))\s*:(\d+)\)")

# Fypp block directives. The duplicating ones (#:for expands to N copies, #:def
# defines a macro instantiated at multiple call sites) collapse many distinct
# generated computations onto a single .fpp source line, so a dd_line hit inside
# one cannot be pinned to a unique runtime instance. #:if/#:with/#:mute select
# code but do not duplicate it, so they are tracked for balance but not flagged.
_FYPP_BLOCK_OPEN = re.compile(r"^\s*#:(for|def|block|call|if|with|mute)\b", re.IGNORECASE)
_FYPP_BLOCK_CLOSE = re.compile(r"^\s*#:end(for|def|block|call|if|with|mute)?\b", re.IGNORECASE)
_FYPP_DUPLICATING = ("for", "def", "block", "call")

# Lines that are clearly control-flow delimiters rather than arithmetic.
# dd_line sometimes reports these when the responsible arithmetic is on the
# preceding line but shares DWARF debug info with the delimiter (e.g. loop
# boundaries in #:for-expanded code, or inlined functions at call sites).
_CONTROL_FLOW_RE = re.compile(
    r"^\s*("
    r"end\s+(do|if|select|where|forall|subroutine|function|module|program|block)\b"
    r"|do\s+\w+\s*=\s*[\w,\s]+"  # naked do-loop header (no arithmetic)
    r"|else(\s+if\s*\(.*\)\s*then)?\s*$"  # else / else if (...) then
    r"|(recursive\s+|pure\s+|elemental\s+)*subroutine\s+\w+"  # subroutine declaration
    r"|\$:END_GPU\w+"  # fypp GPU macro closers
    r"|#:end\w*"  # fypp directive closers (#:endfor, #:enddef, etc.)
    r"|\s*!\s*$"  # comment-only lines
    r"|\s*$"  # blank lines
    r")",
    re.IGNORECASE,
)


def _resolve_source(fname: str, search_whole_tree: bool = False) -> str:
    """Resolve a (possibly bare) source filename to an existing path, or '' if not
    found.  An absolute existing path is used as-is; otherwise the basename is
    located recursively under src/ (then the whole tree if `search_whole_tree`)."""
    if os.path.isabs(fname) and os.path.isfile(fname):
        return fname
    candidates = glob.glob(os.path.join(MFC_ROOT_DIR, "src", "**", os.path.basename(fname)), recursive=True)
    if not candidates and search_whole_tree:
        candidates = glob.glob(os.path.join(MFC_ROOT_DIR, "**", os.path.basename(fname)), recursive=True)
    return candidates[0] if candidates else ""


def _read_source_lines(fname: str, search_whole_tree: bool = False) -> list:
    """Resolve `fname` and return its lines (with newlines), or [] if unreadable."""
    path = _resolve_source(fname, search_whole_tree)
    if not path:
        return []
    try:
        with open(path) as fh:
            return fh.readlines()
    except OSError:
        return []


def _read_source_line(fname: str, lineno: int) -> str:
    """Return the raw source line at lineno (1-based), or '' if unavailable."""
    lines = _read_source_lines(fname)
    return lines[lineno - 1] if 0 < lineno <= len(lines) else ""


def _macro_context_in_lines(lines: list, lineno: int) -> str:
    """Return the innermost code-duplicating fypp block ('#:for'/'#:def'/...) that
    encloses `lineno` (1-based) in `lines`, or None if none does.

    Used to flag dd_line hotspots whose .fpp line is shared across multiple
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


def _ends_with_continuation(line: str) -> bool:
    """True if a free-form Fortran line ends with a continuation '&' (the last
    non-blank token before any trailing comment)."""
    code = line.split("!", 1)[0].rstrip()  # drop trailing comment (string-'!' is rare; fine here)
    return code.endswith("&")


def _statement_bounds_in_lines(lines: list, lineno: int) -> tuple:
    """Return the (start, end) 1-based physical line range of the Fortran logical
    statement containing lineno, following '&' continuations in both directions.

    A hit reported on a continuation fragment thus resolves to the whole
    statement, so the labelled location is the full expression rather than a
    mid-statement piece.
    """
    n = len(lines)
    start = lineno
    while start > 1 and _ends_with_continuation(lines[start - 2]):
        start -= 1
    end = lineno
    while end < n and _ends_with_continuation(lines[end - 1]):
        end += 1
    return start, end


def _statement_at(fname: str, lineno: int) -> tuple:
    """File-backed (start, end, text) for the logical statement at fname:lineno;
    text is the joined statement. Returns (lineno, lineno, '') if unreadable."""
    lines = _read_source_lines(fname)
    if not 0 < lineno <= len(lines):
        return lineno, lineno, ""
    start, end = _statement_bounds_in_lines(lines, lineno)
    # join physical lines, dropping the continuation '&' that may lead or trail each
    text = " ".join(line.strip().strip("&").strip() for line in lines[start - 1 : end])
    return start, end, text


def _is_arithmetic_loc(fname: str, start: int, end: int) -> bool:
    """Return True if any line in [start, end] contains non-trivial arithmetic.

    Filters out loop delimiters and fypp directive lines that dd_line sometimes
    reports when the responsible arithmetic shares DWARF info with its enclosing
    control-flow boundary (inlining, #:for template expansion, etc.).
    Returns True (keep) when uncertain so we never silently drop real hotspots.
    """
    for lineno in range(start, end + 1):
        line = _read_source_line(fname, lineno)
        if not line:
            return True  # can't read — keep to be safe
        if not _CONTROL_FLOW_RE.match(line):
            return True
    return False


def _get_source_context(fname: str, lineno: int, context: int = 2) -> str:
    """Return a annotated source snippet around lineno, or '' if file not found.

    fname may be a bare basename (e.g. 'm_weno.fpp') or a relative path.
    Searches recursively under MFC_ROOT_DIR/src/ first, then the whole tree.
    """
    lines = _read_source_lines(fname, search_whole_tree=True)
    if not lines:
        return ""
    start = max(0, lineno - context - 1)
    end = min(len(lines), lineno + context)
    rows = []
    for i, line in enumerate(lines[start:end], start=start + 1):
        marker = ">" if i == lineno else " "
        rows.append(f"{marker}{i:5d} | {line.rstrip()}")
    return "\n".join(rows)


def _max_diff_np(ref_dir: str, run_dir: str, compare_files: list) -> float:
    import numpy as np

    total = 0.0
    for fname in compare_files:
        ref_p, run_p = os.path.join(ref_dir, fname), os.path.join(run_dir, fname)
        if not os.path.exists(ref_p) or not os.path.exists(run_p):
            return float("inf")
        ref = np.loadtxt(ref_p)[:, 1]
        run = np.loadtxt(run_p)[:, 1]
        total = max(total, float(np.max(np.abs(ref - run))))
    return total


def _max_abs_np(ref_dir: str, compare_files: list) -> float:
    """Return the maximum absolute value across all reference output files."""
    import numpy as np

    total = 0.0
    for fname in compare_files:
        ref_p = os.path.join(ref_dir, fname)
        if not os.path.exists(ref_p):
            continue
        ref = np.loadtxt(ref_p)[:, 1]
        total = max(total, float(np.max(np.abs(ref))))
    return total


def _parse_cancel_gen(gen_path: str) -> list:
    """Parse cc-gen-file TSV (file\\tline\\tsymbol) → sorted unique [(fname, line)] for MFC sources."""
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
# severity (a lower bound — no false positives). 48 ~ full double mantissa.
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


def _parse_rddmin_locs(summary_path: str) -> list:
    """Extract dd_line locations from an rddmin_summary as
    [{path, start, end, macro}] dicts (path is repo-relative; macro is the
    enclosing fypp duplicating block, e.g. '#:for', or None).

    Filters out locations whose source lines are pure control-flow delimiters
    (loop boundaries, fypp directive closers, blank/comment lines).  These can
    appear when the responsible arithmetic shares DWARF debug info with an
    enclosing boundary due to inlining or #:for template expansion.
    """
    if not os.path.isfile(summary_path):
        return []
    locs = []
    skipped = []
    with open(summary_path) as fh:
        for line in fh:
            m = _LOC_RE.search(line)
            if not m:
                continue
            path = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3)) if m.group(3) else start
            try:
                rel = os.path.relpath(path, MFC_ROOT_DIR)
                if rel.startswith(".."):
                    rel = path
            except ValueError:
                rel = path
            rel = rel.replace("\\", "/")
            if _is_arithmetic_loc(path, start, end):
                locs.append({"path": rel, "start": start, "end": end, "macro": _macro_context(path, start)})
            else:
                skipped.append((rel, start, end))
    for rel, start, end in skipped:
        loc = f"{rel}:{start}" if start == end else f"{rel}:{start}-{end}"
        cons.print(f"  [dim]dd_line: skipped control-flow boundary {loc}[/dim]")
    return locs


def _parse_rddmin_syms(summary_path: str) -> list:
    """Extract symbol/function names from a dd_sym rddmin_summary.

    rddmin_summary format:
      ddmin0:\\tFail Ratio: ...\\tFail indexes: ...
      \\t<funcname>\\t<binary_path>
      ddmin1:\\t...
      \\t<funcname>\\t<binary_path>

    Lines starting with 'ddmin' are metadata; function names are on the
    indented (tab-prefixed) lines as the first tab-delimited field.
    """
    if not os.path.isfile(summary_path):
        return []
    syms = []
    with open(summary_path) as fh:
        for ln in fh:
            stripped = ln.strip()
            if not stripped or stripped.startswith("ddmin"):
                continue
            sym = stripped.split("\t")[0].strip()
            if sym:
                syms.append(sym)
    return syms


def _build_source_filter(gen_lines: list, suspect_locs: list) -> list:
    """Select the Verrou --source lines (FILE\\tLINE\\tSYMBOL) that fall on a
    suspect dd_line location.

    gen_lines come from a --gen-source run and carry the exact symbol Verrou
    requires (--source matches on file+line+symbol, not file+line alone).
    suspect_locs are (path, start, end) tuples whose path may be a repo-relative
    path while gen-source emits a basename, so matching is by basename + line.
    """
    ranges = {}
    for path, start, end in suspect_locs:
        ranges.setdefault(os.path.basename(path), []).append((start, end))
    out = []
    for raw in gen_lines:
        parts = raw.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        base = os.path.basename(parts[0].strip())
        try:
            ln = int(parts[1].strip())
        except ValueError:
            continue
        if any(s <= ln <= e for s, e in ranges.get(base, [])):
            out.append(raw if raw.endswith("\n") else raw + "\n")
    return out


def _confirm_decision(suspect_dev, dd_threshold: float):
    """Decide whether perturbing only the suspect lines reproduces the instability.

    Returns True (confirmed), False (suspect lines are inert -> attribution
    suspect, e.g. macro-collapse misattribution), or None if unmeasured.
    """
    if suspect_dev is None:
        return None
    return suspect_dev >= dd_threshold


def _rank_locs(locs: list, total: float) -> list:
    """Attach a 'share' (per-line deviation / total) to each loc dict — which
    must already carry 'share_dev' from a single-line positive control — and
    return the locs sorted by that deviation, most flagrant first.

    'total' is normally float_proxy, so share is the fraction of the full
    single-precision deviation that perturbing that one line alone reproduces.
    A non-positive total yields share=None (cannot normalize).
    """
    for loc in locs:
        dev = loc.get("share_dev")
        loc["share"] = (dev / total) if (dev is not None and total and total > 0) else None
    return sorted(locs, key=lambda loc: (loc.get("share_dev") or 0.0), reverse=True)


def _mark_cancellation(dd_line_locs: list, cancellation_locs: list) -> list:
    """Set loc['cancellation']=True for each dd_line loc whose line range covers a
    catastrophic-cancellation site (stage F), matched by basename + line.

    This pins the flagrant operation on a multi-op line to the subtraction that
    cancels, rather than just naming the line.
    """
    by_base = {}
    for fname, lineno in cancellation_locs:
        by_base.setdefault(os.path.basename(fname), set()).add(lineno)
    for loc in dd_line_locs:
        lines = by_base.get(os.path.basename(loc["path"]), set())
        loc["cancellation"] = any(ln in lines for ln in range(loc["start"], loc["end"] + 1))
    return dd_line_locs


def _cancellation_by_file(cancellation_locs: list) -> list:
    """Aggregate cancellation sites by source file → [(basename, count)] sorted by
    count (desc), ties by name.

    This is the cancellation-*origin* view (where ill-conditioning concentrates),
    as opposed to the per-line --source share, which is a *sensitivity* view
    (where reduced precision most moves the output — typically the time
    integrator / final accumulation, regardless of where error originates).
    """
    counts = {}
    for fname, _lineno in cancellation_locs:
        base = os.path.basename(fname)
        counts[base] = counts.get(base, 0) + 1
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
