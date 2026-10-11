"""`./mfc.sh context`: what to re-read to check a diff to MFC's Fortran.

Works on Fypp-expanded sources (context_expand), whose '# N "file"' markers map every expanded line back to its
source line, so a changed source line is traced to every generated copy of it.
"""

import collections
import re
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
    srcs = [(s, l) for _, s, l in r.rows if s.endswith((".fpp", ".f90")) and "/include/" not in s]
    r.src_count = collections.Counter(srcs)
    if srcs:
        r.src_file = collections.Counter(s for s, _ in srcs).most_common(1)[0][0]
        own = [l for s, l in srcs if s == r.src_file]
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


def parse_diff(diff_text, repo, side):
    """Changed source lines per file (side 'new' or 'old'), changed text per file, and files the diff adds."""
    changed, text, new = collections.defaultdict(set), collections.defaultdict(list), []
    cur, old_is_null = None, False
    for line in diff_text.splitlines():
        if line.startswith("--- "):
            old_is_null = line.endswith("/dev/null")
            if side == "old":
                cur = None if old_is_null else str(Path(repo) / line[6:])
        elif line.startswith("+++ "):
            path = None if line.endswith("/dev/null") else str(Path(repo) / line[6:])
            if old_is_null and path:
                new.append(path)
            if side == "new":
                cur = path
            if cur and "/src/" not in cur:
                cur = None
        elif line.startswith("@@") and cur:
            m = re.search(r"\+(\d+)(?:,(\d+))?" if side == "new" else r"-(\d+)(?:,(\d+))?", line)
            start, n = int(m.group(1)), int(m.group(2) or 1)
            changed[cur].update(range(start, start + max(n, 1)))
        elif cur and line[:1] in "+-":
            text[cur].append(line[1:])
    new = [p for p in new if "/src/" in p]
    for p in new:
        changed.pop(p, None)
    return dict(changed), dict(text), new


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
    mine = sorted(l for l in changed.get(r0.src_file, ()) if lo <= l <= hi)
    a, b, note = lo, hi, []
    if hi - lo > 300 and mine:
        a = max((i for i in range(lo, mine[0] + 1) if KERNEL_START in src[i - 1]), default=max(lo, mine[0] - 40))
        b = min((i for i in range(mine[-1], hi + 1) if KERNEL_END in src[i - 1]), default=min(hi, mine[-1] + 40))
        note = [f"(routine is {hi - lo + 1} lines; lines {a}-{b} hold the change, the rest is at {_rel(r0.src_file, repo)}:{lo}-{hi})"]
    return [*note, "```fortran", *src[a - 1 : b], "```"]


def _globals(idx, names, by_name, repo):
    lines = set()
    for n in names:
        r = by_name[n][0]
        uses = idx.module_uses.get((r.target, r.module), {})
        for v in sorted((r.idents - r.locals) & set(idx.module_vars)):
            owners = {m for m in idx.module_vars[v] if m == r.module or (m in uses and (uses[m] is None or v in uses[m]))}
            if not owners:
                continue
            writers = sorted({f"`{w.name}` ({_rel(w.src_file, repo)}:{w.src_range[0]})" for w in idx.routines if w.target == r.target and v in w.assigns and v not in w.locals and w.name.lower() not in names and not SETUP.search(w.name)})
            if writers:
                more = f", +{len(writers) - 4} more" if len(writers) > 4 else ""
                lines.add(f"- `{v}` ({', '.join(sorted(owners))}) is set by {', '.join(writers[:4])}{more}")
    return sorted(lines)


def render(idx, changed, changed_text, new_files, repo, max_inline=5, read_source=None):  # pylint: disable=too-many-arguments,too-many-locals
    """Markdown: what to re-read to check the diff."""
    read_source = read_source or (lambda f: Path(f).read_text().splitlines())
    if not changed and not new_files:
        return "No changes under src/."
    by_name = collections.defaultdict(list)
    for r in idx.routines:
        by_name[r.name.lower()].append(r)
    names = sorted({r.name.lower() for r in idx.routines if any((f, l) in r.src_count for f, ls in changed.items() for l in ls)})
    out = ["# What to re-read for this diff", ""]
    if new_files:
        out += ["## New files: read in full", *(f"- {_rel(p, repo)}" for p in new_files), ""]
    out += macro_sections(changed, changed_text, repo)
    if not names:
        return "\n".join(out)
    out.append("## Changed routines")
    for k, n in enumerate(names):
        r0 = by_name[n][0]
        targets = ", ".join(sorted({r.target for r in by_name[n]}))
        out.append(f"\n### `{r0.name}` {_rel(r0.src_file, repo)}:{r0.src_range[0]}-{r0.src_range[1]} (in {targets})")
        for l in sorted(changed.get(r0.src_file, ())):
            if (copies := r0.src_count.get((r0.src_file, l), 0)) > 1:
                out.append(f"- line {l} is emitted {copies} times (Fypp loop): check every copy, each kernel has its own clauses")
                break
        out += ["(not shown inline; read the range above)"] if k >= max_inline else _inline(r0, changed, read_source, repo)
    tgts = {r.target for n in names for r in by_name[n]}
    callers = sorted({f"- `{r.name}` {_rel(r.src_file, repo)}:{r.src_range[0]}" for r in idx.routines if r.calls & set(names) and r.name.lower() not in names})
    callees = sorted({f"- `{c}` {_rel(x.src_file, repo)}:{x.src_range[0]}" for n in names for r in by_name[n] for c in r.calls if c not in names for x in by_name.get(c, [])[:1] if x.target in tgts})
    out += ["\n## Callers", *(_cap(callers) or ["- none found"])]
    out += ["\n## Callees", *(_cap(callees) or ["- none found"])]
    out += ["\n## Module variables used here and set elsewhere at run time", *(_cap(_globals(idx, names, by_name, repo)) or ["- none"])]
    return "\n".join(out)


def macro_sections(changed, changed_text, repo):  # pylint: disable=unused-argument
    return []
