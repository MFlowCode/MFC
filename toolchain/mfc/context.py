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
