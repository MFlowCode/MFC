"""`./mfc.sh context`: what to re-read to check a diff to MFC's Fortran.

Works on Fypp-expanded sources (context_expand), whose '# N "file"' markers map every expanded line back to its
source line, so a changed source line is traced to every generated copy of it.
"""

import collections
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

MARKER = re.compile(r'^#\s*(\d+)\s+"([^"]+)"')
ROUTINE = re.compile(r"^\s*(?:(?:pure|elemental|impure|recursive|module)\s+)*(?:(?:real|integer|logical|complex|type\s*\([^)]*\))\S*\s+)?(subroutine|function)\s+(\w+)", re.I)
END_ROUTINE = re.compile(r"^\s*end\s+(subroutine|function)\b", re.I)
MODULE = re.compile(r"^\s*module\s+(?!procedure\b)(\w+)\s*$", re.I)
TYPE_START = re.compile(r"^\s*type\s*(?:,[^:]*)?(?:::)?\s*(\w+)\s*$", re.I)
TYPE_END = re.compile(r"^\s*end\s*type\b", re.I)
CONTAINS = re.compile(r"^\s*contains\s*$", re.I)
DECL = re.compile(r"^\s*(real|integer|logical|complex|character|type\s*\()[^!]*::\s*(.*)$", re.I)
USE = re.compile(r"^\s*use\s+(\w+)\s*(?:,\s*only\s*:\s*(.*))?$", re.I)
CALL = re.compile(r"\bcall\s+(\w+)", re.I)
FUNC_REF = re.compile(r"\b(f_\w+)\s*\(", re.I)
ASSIGN = re.compile(r"^\s*(\w+)\s*(?:\([^=]*\))?\s*(?:%[\w%()]*)?\s*=(?!=)")


def _code(text: str) -> str:
    return text if "!$" in text else text.split("!")[0]


def _decl_names(rhs: str) -> list:
    rhs = re.sub(r"\([^()]*\)", "", re.sub(r"\([^()]*\)", "", rhs))
    return [m.group(1).lower() for tok in rhs.split(",") if (m := re.match(r"\s*(\w+)", tok))]


def read_expanded(path):
    """(text, source file, source line) for each expanded line, from the Fypp line markers."""
    rows, src, ln = [], "", 0
    for text in Path(path).read_text(errors="replace").splitlines():
        if m := MARKER.match(text):
            ln, src = int(m.group(1)), m.group(2)
            continue
        rows.append((text, src, ln))
        ln += 1
    return rows


@dataclass
class Routine:  # pylint: disable=too-many-instance-attributes
    name: str
    target: str
    module: str
    rows: list = field(default_factory=list)
    calls: set = field(default_factory=set)
    idents: set = field(default_factory=set)
    assigns: set = field(default_factory=set)
    locals: set = field(default_factory=set)
    src_file: str = ""
    src_range: tuple = (0, 0)
    src_count: collections.Counter = field(default_factory=collections.Counter)
    src_runs: collections.Counter = field(default_factory=collections.Counter)  # separate appearances = copies


@dataclass
class Index:
    routines: list
    module_vars: dict  # name -> modules declaring it
    module_uses: dict  # (target, module) -> {used module: only-list, or None for all}


def _finish(r: Routine):
    code = [_code(t) for t, _, _ in r.rows]
    body = "\n".join(code)
    r.calls = {c.lower() for c in CALL.findall(body)} | {c.lower() for c in FUNC_REF.findall(body)}
    r.idents = {w.lower() for w in re.findall(r"(?<![%\w])[a-z_]\w*", body, re.I)}
    r.assigns = {m.group(1).lower() for c in code if (m := ASSIGN.match(c))}
    r.locals = {n for c in code if (d := DECL.match(c)) for n in _decl_names(d.group(2))}
    srcs = [(s, ln) for _, s, ln in r.rows if s.endswith((".fpp", ".f90")) and not re.search(r"/src/\w+/include/", s)]
    r.src_count = collections.Counter(srcs)
    prev = None
    for _, s, ln in r.rows:
        if (s, ln) != prev:
            r.src_runs[(s, ln)] += 1
        prev = (s, ln)
    if srcs:
        r.src_file = collections.Counter(s for s, _ in srcs).most_common(1)[0][0]
        own = [ln for s, ln in srcs if s == r.src_file]
        r.src_range = (min(own), max(own))


def build_index(files):
    """Index routines, module variables and use-lists of the expanded files of each target."""
    routines, module_vars, module_uses = [], {}, collections.defaultdict(dict)
    for target, paths in files.items():
        for path in paths:
            stack, module, in_spec, in_type = [], None, False, 0
            for text, src, ln in read_expanded(path):
                code = _code(text)
                if m := MODULE.match(code):
                    module, in_spec = m.group(1).lower(), True
                elif in_spec and not stack and CONTAINS.match(code):
                    in_spec = False
                if module and not stack and (u := USE.match(code)):
                    only = {n.split("=>")[-1].strip().lower() for n in u.group(2).split(",")} if u.group(2) else None
                    module_uses[(target, module)][u.group(1).lower()] = only
                if TYPE_START.match(code):
                    in_type += 1
                elif TYPE_END.match(code):
                    in_type = max(0, in_type - 1)
                elif in_spec and not stack and not in_type and (d := DECL.match(code)) and "parameter" not in code.lower().split("::")[0]:
                    for n in _decl_names(d.group(2)):
                        module_vars.setdefault(n, set()).add(module)
                if (m := ROUTINE.match(code)) and not re.match(r"\s*end\b", code, re.I):
                    stack.append(Routine(m.group(2), target, module or ""))
                for r in stack:
                    r.rows.append((text, src, ln))
                if END_ROUTINE.match(code) and stack:
                    r = stack.pop()
                    _finish(r)
                    routines.append(r)
    return Index(routines, module_vars, dict(module_uses))


CAP = 15
KERNEL_START, KERNEL_END = "GPU_PARALLEL_LOOP(", "END_GPU_PARALLEL_LOOP"
# Routines that only read, default or broadcast inputs; listing them as writers of a global is noise.
SETUP = re.compile(r"(s_read_input_file|bcast|default|s_initialize_global_parameters_module|s_read_serial)", re.I)


@dataclass
class Diff:
    changed: dict = field(default_factory=dict)  # absolute source path -> changed line numbers (on the mapped side)
    text: dict = field(default_factory=dict)  # absolute source path -> text of its + and - lines
    new: list = field(default_factory=list)  # files the diff adds
    deleted: list = field(default_factory=list)  # files the diff removes
    checks: dict = field(default_factory=dict)  # absolute source path -> [(line, text)] the mapped side should contain


def parse_diff(diff_text, repo, side):
    """Map a unified diff onto source lines of one side ('new': the working tree, 'old': the base)."""
    d = Diff(changed=collections.defaultdict(set), text=collections.defaultdict(list), checks=collections.defaultdict(list))
    cur, old_path = None, None

    def src_path(header):
        rel = header[6:].rstrip("\t")
        return None if header.endswith("/dev/null") or not rel.startswith("src/") else str(Path(repo) / rel)

    for line in diff_text.splitlines():
        if line.startswith("--- "):
            old_path = src_path(line) if not line.endswith("/dev/null") else "/dev/null"
            if side == "old":
                cur = old_path if old_path != "/dev/null" else None
        elif line.startswith("+++ "):
            new_path = src_path(line)
            if old_path == "/dev/null" and new_path:
                d.new.append(new_path)
            if line.endswith("/dev/null") and old_path not in (None, "/dev/null"):
                d.deleted.append(old_path)
            if side == "new":
                cur = new_path
        elif line.startswith("@@") and cur:
            m = re.match(r"@@ -(\d+)(?:,(\d+))? \+(\d+)", line)
            old_ln, new_ln = int(m.group(1)), int(m.group(3))
            if m.group(2) == "0":
                old_ln += 1  # git numbers a pure insertion by the line it follows
        elif cur and (not line or line[0] in "+- "):
            kind, body = (line[0], line[1:]) if line else (" ", "")
            if kind == "-":
                d.changed[cur].add(old_ln if side == "old" else max(new_ln, 1))  # a deletion maps to the line now in its place
                if side == "old":
                    d.checks[cur].append((old_ln, body))
                old_ln += 1
            elif kind == "+":
                if side == "new":
                    d.changed[cur].add(new_ln)
                    d.checks[cur].append((new_ln, body))
                elif old_ln > 1:
                    d.changed[cur].add(old_ln - 1)  # an insertion maps to the old line it follows
                new_ln += 1
            else:
                d.checks[cur].append((new_ln if side == "new" else old_ln, body))
                old_ln, new_ln = old_ln + 1, new_ln + 1
            if kind != " ":
                d.text[cur].append(body)
    for p in d.new + d.deleted:
        d.changed.pop(p, None)
    d.changed, d.text, d.checks = dict(d.changed), dict(d.text), dict(d.checks)
    return d


def _rel(path, repo):
    try:
        return str(Path(path).relative_to(Path(repo)))
    except ValueError:
        return str(path)


def _cap(items):
    return items[:CAP] + ([f"- +{len(items) - CAP} more"] if len(items) > CAP else [])


def _inline(r0, changed, read_source, repo):
    """The routine's source, or the kernel(s) holding the change when the routine is long."""
    lo, hi = r0.src_range
    src = read_source(r0.src_file)
    mine = sorted(ln for ln in changed.get(r0.src_file, ()) if lo <= ln <= hi)
    a, b, note = lo, hi, []
    if hi - lo > 300 and mine:
        a = max((i for i in range(lo, mine[0] + 1) if KERNEL_START in src[i - 1]), default=max(lo, mine[0] - 40))
        b = min((i for i in range(mine[-1], hi + 1) if KERNEL_END in src[i - 1]), default=min(hi, mine[-1] + 40))
        note = [f"(routine is {hi - lo + 1} lines; lines {a}-{b} hold the change, the rest is at {_rel(r0.src_file, repo)}:{lo}-{hi})"]
    return [*note, "```fortran", *src[a - 1 : b], "```"]


def _globals(idx, changed_routines, repo):
    lines, changed_names = set(), {r.name.lower() for r in changed_routines}
    for r in changed_routines:
        uses = idx.module_uses.get((r.target, r.module), {})
        for v in sorted((r.idents - r.locals) & set(idx.module_vars)):
            owners = {m for m in idx.module_vars[v] if m == r.module or (m in uses and (uses[m] is None or v in uses[m]))}
            if not owners:
                continue
            writers = sorted(
                {
                    f"`{w.name}` ({_rel(w.src_file, repo)}:{w.src_range[0]})"
                    for w in idx.routines
                    if w.target == r.target and v in w.assigns and v not in w.locals and w.name.lower() not in changed_names and not SETUP.search(w.name)
                }
            )
            if writers:
                more = f", +{len(writers) - 4} more" if len(writers) > 4 else ""
                lines.add(f"- `{v}` ({', '.join(sorted(owners))}) is set by {', '.join(writers[:4])}{more}")
    return sorted(lines)


def _statement_starts(changed, read_source):
    """Map each changed line to the first line of its statement: Fypp marks only first lines of continued statements."""
    out = {}
    for f, ls in changed.items():
        try:
            src = read_source(f)
        except OSError:
            src = []
        starts = {}
        for ln in ls:
            start = ln
            while 1 < start <= len(src) and src[start - 2].rstrip().endswith("&"):
                start -= 1
            starts[ln] = start
        out[f] = starts
    return out


def _mismatches(diff, read_source, repo):
    """Hunks whose lines differ from the file being mapped onto: line numbers there would point at other code."""
    out = []
    for f, checks in diff.checks.items():
        try:
            src = read_source(f)
        except OSError:
            continue
        bad = [ln for ln, body in checks if ln > len(src) or src[ln - 1].rstrip() != body.rstrip()]
        if bad:
            out.append(f"- {_rel(f, repo)}: the diff does not match this file at line {bad[0]} ({len(bad)} lines); results for it may point at other code")
    return out


def render(idx, diff, repo, max_inline=5, read_source=None):  # pylint: disable=too-many-locals
    """Markdown: what to re-read to check the diff."""
    read_source = read_source or (lambda f: Path(f).read_text().splitlines())
    if not diff.changed and not diff.new and not diff.deleted:
        return "No changes under src/."
    starts = _statement_starts(diff.changed, read_source)
    changed = {f: set(m) | set(m.values()) for f, m in starts.items()}
    new_files = diff.new
    # A routine is its name and source file: each target's m_global_parameters has its own same-named routines.
    by_key = collections.defaultdict(list)
    for r in idx.routines:
        by_key[(r.name.lower(), r.src_file)].append(r)
    keys = sorted(k for k, rs in by_key.items() if any((f, ln) in rs[0].src_count for f, ls in changed.items() for ln in ls))
    out = ["# What to re-read for this diff", ""]
    if bad := _mismatches(diff, read_source, repo):
        out += ["## Diff does not match the source", *bad, ""]
    if new_files:
        out += ["## New files: read in full", *(f"- {_rel(p, repo)}" for p in new_files), ""]
    if diff.deleted:
        out += ["## Deleted files: check what used them", *(f"- {_rel(p, repo)}" for p in diff.deleted), ""]
    out += macro_sections(changed, diff.text, repo)
    covered = {(f, ln) for key in keys for (f, ln) in by_key[key][0].src_count}
    loose = {f: sorted(ln for ln, start in m.items() if (f, start) not in covered) for f, m in starts.items() if "/include/" not in _rel(f, repo)}
    if loose := {f: ls for f, ls in loose.items() if ls}:
        out += ["## Not in any routine (module declarations, or code this build configuration leaves out): read these lines directly"]
        out += [f"- {_rel(f, repo)}: {', '.join(map(str, ls))}" for f, ls in sorted(loose.items())] + [""]
    if not keys:
        return "\n".join(out)
    out.append("## Changed routines")
    for k, key in enumerate(keys):
        r0 = by_key[key][0]
        targets = ", ".join(sorted({r.target for r in by_key[key]}))
        out.append(f"\n### `{r0.name}` {_rel(r0.src_file, repo)}:{r0.src_range[0]}-{r0.src_range[1]} (in {targets})")
        for ln in sorted(changed.get(r0.src_file, ())):
            if (copies := r0.src_runs.get((r0.src_file, ln), 0)) > 1:
                out.append(f"- line {ln} is emitted {copies} times (Fypp loop): check every copy, each kernel has its own clauses")
                break
        out += ["(not shown inline; read the range above)"] if k >= max_inline else _inline(r0, changed, read_source, repo)
    changed_routines = [r for key in keys for r in by_key[key]]
    tgts = {r.target for r in changed_routines}
    names = {r.name.lower() for r in changed_routines}
    callers = sorted({f"- `{r.name}` {_rel(r.src_file, repo)}:{r.src_range[0]}" for r in idx.routines if r.target in tgts and r.calls & names and (r.name.lower(), r.src_file) not in keys})
    callees = sorted(
        {
            f"- `{x.name}` {_rel(x.src_file, repo)}:{x.src_range[0]}"
            for r in changed_routines
            for x in idx.routines
            if x.target == r.target and x.name.lower() in r.calls and x.name.lower() not in names
        }
    )
    out += ["\n## Callers", *(_cap(callers) or ["- none found"])]
    out += ["\n## Callees", *(_cap(callees) or ["- none found"])]
    out += ["\n## Module variables used here and set elsewhere at run time", *(_cap(_globals(idx, changed_routines, repo)) or ["- none"])]
    return "\n".join(out)


DEF = re.compile(r"^\s*#:def\s+(\w+)")


def _defs(lines):
    out, start, name = [], None, None
    for i, t in enumerate(lines, 1):
        if m := DEF.match(t):
            start, name = i, m.group(1)
        elif re.match(r"\s*#:enddef", t) and start:
            out.append((name, start, i))
            start = None
    return out


def _joined(text):
    """Lines with '&' continuations joined, so a macro call or signature is one line."""
    return re.sub(r"&\s*\n\s*&?", " ", text).splitlines()


def _callers_of(name, incs):
    """The macro and every macro that calls it, directly or indirectly."""
    found, todo = {name}, [name]
    while todo:
        cur = todo.pop()
        for g in incs:
            outer = None
            for t in g.read_text().splitlines():
                if m := DEF.match(t):
                    outer = m.group(1)
                elif outer and outer not in found and re.search(rf"\b{cur}\s*\(", t):
                    found.add(outer)
                    todo.append(outer)
    return found


def macro_sections(changed, changed_text, repo):
    """Call sites of each changed Fypp macro.

    Fypp line markers record where a macro is defined, not where it expands, so call sites are found by name,
    following wrapper macros (GPU_PARALLEL_LOOP -> OMP_PARALLEL_LOOP). When the change mentions some of the
    macro's arguments, only call sites that pass them are kept.
    """
    repo = Path(repo)
    incs = sorted(repo.glob("src/*/include/*.fpp"))
    out = []
    for f, ls in changed.items():
        if not f.endswith(".fpp") or not Path(f).exists():
            continue
        src = Path(f).read_text().splitlines()
        for name, a, b in _defs(src):
            if not any(a <= ln <= b for ln in ls):
                continue
            callers = _callers_of(name, incs)
            sig = re.match(r"\s*#:def\s+\w+\s*\(([^)]*)\)", _joined("\n".join(src[a - 1 : b]))[0])
            args = [x.split("=")[0].strip() for x in sig.group(1).split(",")] if sig else []
            touched = [x for x in args if x and not x.startswith("*") and any(re.search(rf"\b{x}\b", t) for t in changed_text.get(f, []))]
            pat = re.compile(r"\b(" + "|".join(sorted(callers)) + r")\s*\(")
            byfile = collections.Counter()
            for g in sorted(repo.glob("src/**/*.fpp")):
                if "/include/" in str(g):
                    continue
                calls = [t for t in _joined(g.read_text()) if re.search(r"[$@]:", t) and pat.search(t)]
                if touched:
                    calls = [t for t in calls if any(re.search(rf"\b{x}\s*=", t) for x in touched)]
                if calls:
                    byfile[_rel(g, repo)] = len(calls)
            via = ", ".join(sorted(callers - {name})) or "no other macro"
            scope = f" passing {', '.join(touched)}" if touched else ""
            out.append(f"## Changed macro `{name}` ({_rel(f, repo)}:{a}-{b})")
            out.append(f"Used directly or through {via}: {sum(byfile.values())} call sites{scope} in {len(byfile)} files.")
            out += [*_cap([f"- {k}: {v}" for k, v in byfile.most_common()]), ""]
    return out


SOLVER_TARGETS = ("pre_process", "simulation", "post_process")


def default_base(repo):
    """The merge base with origin/master, or HEAD (with a note) when there is no such remote."""
    r = subprocess.run(["git", "merge-base", "HEAD", "origin/master"], cwd=repo, capture_output=True, text=True, check=False)
    if r.returncode == 0 and r.stdout.strip():
        return r.stdout.strip(), None
    return "HEAD", "no origin/master: showing uncommitted changes against HEAD"


def git_diff(repo, base):
    """`git diff` of src/ against base, immune to user diff settings (prefixes, colour, external tools)."""
    cmd = ["git", "-c", "diff.noprefix=false", "diff", "--no-color", "--no-ext-diff", "--src-prefix=a/", "--dst-prefix=b/", base, "--", "src"]
    return subprocess.run(cmd, cwd=repo, capture_output=True, text=True, check=True).stdout


def untracked_sources(repo):
    out = subprocess.run(["git", "ls-files", "--others", "--exclude-standard", "--", "src"], cwd=repo, capture_output=True, text=True, check=False).stdout
    return [str(Path(repo) / p) for p in out.split()]


def context():
    """Entry point of `./mfc.sh context`."""

    from . import common, context_expand  # pylint: disable=import-outside-toplevel
    from .state import ARG  # pylint: disable=import-outside-toplevel

    repo = Path(ARG("root") or common.MFC_ROOT_DIR).resolve()
    side = "old" if ARG("expand") == "base" else "new"
    note = None
    if ARG("diff"):
        diff = Path(ARG("diff")).read_text()
    else:
        base, note = (ARG("base"), None) if ARG("base") else default_base(repo)
        diff = git_diff(repo, base)
    d = parse_diff(diff, repo, side)
    if not ARG("diff"):
        d.new += [p for p in untracked_sources(repo) if p not in d.new]
    targets = [t for t in ARG("targets") if t in SOLVER_TARGETS] or list(SOLVER_TARGETS)
    files, errors = context_expand.expand(repo, targets, repo / "build" / "context") if d.changed else ({}, [])
    report = render(build_index(files), d, repo)
    if note:
        report = f"({note})\n\n{report}"
    if errors:
        report += "\n\n## Not expanded\n" + "\n".join(f"- {e}" for e in errors)
    if ARG("output"):
        Path(ARG("output")).write_text(report + "\n")
    else:
        sys.stdout.write(report + "\n")
