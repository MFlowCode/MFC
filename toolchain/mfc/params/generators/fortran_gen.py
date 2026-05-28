"""Fortran parameter code generator — namelist and scalar decl fragments per target."""

import re
from typing import List

import mfc.params.definitions  # noqa: F401 - triggers registry population

from ..namelist_targets import CASE_OPT_EXCLUDE, NAMELIST_VARS
from ..registry import REGISTRY
from ..schema import ParamDef, ParamType

_HEADER = "! AUTO-GENERATED - do not edit directly. Regenerate: cmake reconfigure\n!\n"

_MAX_LINE = 130
_FIRST_PREFIX = "namelist /user_inputs/ "
_CONT_PREFIX = "    & "
_CONT2_PREFIX = "        & "  # inside #:if block

_DECL_COL = 24  # '::' column, matches ffmt alignment


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
    return "%" not in name and "(" not in name


def _vars_for_target(target: str) -> List[str]:
    return sorted(v for v, ts in NAMELIST_VARS.items() if target in ts)


def _pack_namelist(vars_list: List[str], first_prefix: str, cont_prefix: str, max_line: int) -> List[str]:
    """Pack variable names into Fortran continuation lines; all but last end with ', &'."""
    if not vars_list:
        return []
    lines: List[str] = []
    prefix = first_prefix
    current_vars: List[str] = []
    current_len = len(prefix)
    for var in vars_list:
        additional = len(var) + (2 if current_vars else 0)
        if current_vars and current_len + additional + 3 > max_line:
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
    return "\n".join(_pack_namelist(vars_list, _FIRST_PREFIX, _CONT_PREFIX, _MAX_LINE))


def generate_namelist_fpp(target: str) -> str:
    """Return the namelist /user_inputs/ statement for a target as a string."""
    assert target in ("pre", "sim", "post")
    all_vars = _vars_for_target(target)

    if target != "sim":
        return _HEADER + _format_namelist(all_vars) + "\n"

    normal = [v for v in all_vars if v not in CASE_OPT_EXCLUDE]
    opt = sorted(v for v in CASE_OPT_EXCLUDE if v in NAMELIST_VARS and "sim" in NAMELIST_VARS[v])

    nl_lines = _pack_namelist(normal, _FIRST_PREFIX, _CONT_PREFIX, _MAX_LINE)
    nl_lines[-1] += ", &"
    opt_lines = _pack_namelist(opt, _CONT_PREFIX, _CONT2_PREFIX, _MAX_LINE)

    parts = [_HEADER.rstrip()] + nl_lines + ["#:if not MFC_CASE_OPTIMIZATION"] + opt_lines + ["#:endif"]
    return "\n".join(parts) + "\n"


def generate_decls_fpp(target: str) -> str:
    """Return simple scalar Fortran declarations for a target as a string."""
    assert target in ("pre", "sim", "post")
    all_params = REGISTRY.all_params
    lines = [_HEADER.rstrip()]
    for name in _vars_for_target(target):
        if not _is_simple_scalar(name):
            continue
        param = all_params.get(name)
        if param is None:
            continue
        lines.append(f"{fortran_type_decl(param).ljust(_DECL_COL)}:: {name}")
    return "\n".join(lines) + "\n"
