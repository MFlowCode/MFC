"""
Fortran parameter code generator.

Generates namelist fragments and simple scalar declaration fragments
per target (pre/sim/post). Output consumed by generate.py.

Output format matches ffmt (the MFC Fortran formatter) so that
./mfc.sh format is idempotent on these generated files.
"""

import re
from typing import List

import mfc.params.definitions  # noqa: F401 - triggers registry population

from ..namelist_targets import CASE_OPT_EXCLUDE, NAMELIST_VARS
from ..registry import REGISTRY
from ..schema import ParamDef, ParamType

# ffmt collapses the two header lines into one and uses ASCII hyphen
_HEADER = "! AUTO-GENERATED - do not edit directly. Regenerate: ./mfc.sh generate\n!\n"

# ffmt formats Fortran with max 130-char lines and 4-space continuation indent
_MAX_LINE = 130
_FIRST_PREFIX = "namelist /user_inputs/ "
_CONT_PREFIX = "    & "
_CONT2_PREFIX = "        & "  # second-level continuation (inside Fypp #:if block)

# ffmt aligns '::' to a fixed column; widest type is character(LEN=path_len) = 23 chars
_DECL_COL = 24  # pad type string to this width before '::'


def get_namelist_var(param_name: str) -> str:
    """Return the Fortran namelist root for a parameter name."""
    m = re.match(r"^([a-zA-Z_]\w*)\(", param_name)
    if m and "%" in param_name:
        return m.group(1)
    if "%" in param_name:
        return param_name.split("%", maxsplit=1)[0]
    return param_name


def fortran_type_decl(param: ParamDef) -> str:
    """Return the Fortran type string for a parameter."""
    mapping = {
        ParamType.INT: "integer",
        ParamType.REAL: "real(wp)",
        ParamType.LOG: "logical",
        ParamType.ANALYTIC_INT: "integer",
        ParamType.ANALYTIC_REAL: "real(wp)",
    }
    if param.param_type == ParamType.STR:
        return f"character(LEN={param.str_len})"
    return mapping[param.param_type]


def _is_simple_scalar(name: str) -> bool:
    """Return True if name has no '%' and no '(' - i.e. a plain simple variable."""
    return "%" not in name and "(" not in name


def _vars_for_target(target: str) -> List[str]:
    """Return sorted list of namelist variable names for the given target."""
    return sorted(v for v, ts in NAMELIST_VARS.items() if target in ts)


def _pack_namelist(vars_list: List[str], first_prefix: str, cont_prefix: str, max_line: int) -> List[str]:
    """
    Pack a list of variable names into Fortran namelist continuation lines.

    Returns a list of lines WITHOUT trailing newlines.
    All lines except the last end with ', &'.
    """
    if not vars_list:
        return []

    lines: List[str] = []
    prefix = first_prefix
    current_vars: List[str] = []
    current_len = len(prefix)

    for var in vars_list:
        additional = len(var) + (2 if current_vars else 0)
        if current_vars and current_len + additional + 3 > max_line:
            # Flush with continuation marker
            lines.append(prefix + ", ".join(current_vars) + ", &")
            prefix = cont_prefix
            current_vars = [var]
            current_len = len(cont_prefix) + len(var)
        else:
            current_vars.append(var)
            current_len += additional

    if current_vars:
        lines.append(prefix + ", ".join(current_vars))

    return lines


def _format_namelist(vars_list: List[str]) -> str:
    """Format vars as a Fortran namelist statement block (no trailing newline)."""
    lines = _pack_namelist(vars_list, _FIRST_PREFIX, _CONT_PREFIX, _MAX_LINE)
    return "\n".join(lines)


def generate_namelist_fpp(target: str) -> str:
    """Generate the namelist /user_inputs/ statement for a target."""
    assert target in ("pre", "sim", "post")
    all_vars = _vars_for_target(target)

    if target != "sim":
        return _HEADER + _format_namelist(all_vars) + "\n"

    # For sim: split into normal vars and case-opt-excluded vars
    normal = [v for v in all_vars if v not in CASE_OPT_EXCLUDE]
    opt = sorted(v for v in CASE_OPT_EXCLUDE if v in NAMELIST_VARS and "sim" in NAMELIST_VARS[v])

    # Normal vars: last line gets ', &' since opt vars follow
    nl_lines = _pack_namelist(normal, _FIRST_PREFIX, _CONT_PREFIX, _MAX_LINE)
    nl_lines[-1] += ", &"

    # Opt vars: pack using cont_prefix for first line, cont2_prefix for subsequent
    opt_lines = _pack_namelist(opt, _CONT_PREFIX, _CONT2_PREFIX, _MAX_LINE)

    all_lines = [_HEADER.rstrip()] + nl_lines + ["#:if not MFC_CASE_OPTIMIZATION"] + opt_lines + ["#:endif"]
    return "\n".join(all_lines) + "\n"


def generate_decls_fpp(target: str) -> str:
    """Generate simple scalar Fortran variable declarations for a target.

    Column-aligns '::' to match ffmt output (type padded to _DECL_COL chars).
    """
    assert target in ("pre", "sim", "post")
    all_params = REGISTRY.all_params
    vars_for_target = _vars_for_target(target)
    lines = [_HEADER.rstrip()]
    for name in vars_for_target:
        if not _is_simple_scalar(name):
            continue
        param = all_params.get(name)
        if param is None:
            continue
        type_str = fortran_type_decl(param)
        padded = type_str.ljust(_DECL_COL)
        lines.append(f"{padded}:: {name}")
    return "\n".join(lines) + "\n"
