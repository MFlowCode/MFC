"""
Fortran parameter code generator.

Generates namelist fragments and simple scalar declaration fragments
per target (pre/sim/post). Output consumed by generate.py.
"""

import re
from pathlib import Path
from typing import List, Tuple

import mfc.params.definitions  # noqa: F401 — triggers registry population

from ..namelist_targets import CASE_OPT_EXCLUDE, NAMELIST_VARS
from ..registry import REGISTRY
from ..schema import ParamDef, ParamType

_HEADER = (
    "! AUTO-GENERATED — do not edit directly.\n"
    "! Regenerate: ./mfc.sh generate\n"
    "!\n"
)


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
    """Return True if name has no '%' and no '(' — i.e. a plain simple variable."""
    return "%" not in name and "(" not in name


def _vars_for_target(target: str) -> List[str]:
    """Return sorted list of namelist variable names for the given target."""
    return sorted(v for v, ts in NAMELIST_VARS.items() if target in ts)


def _format_namelist(vars_list: List[str]) -> str:
    """
    Format a list of variable names as a Fortran namelist continuation block.

    The first line starts with 'namelist /user_inputs/ '.
    Continuation lines use '    & ' prefix (4 spaces + ampersand + space).
    All lines except the last end with ', &' for continuation.
    Lines are wrapped at 80 characters.
    """
    if not vars_list:
        return ""

    FIRST_PREFIX = "namelist /user_inputs/ "
    CONT_PREFIX = "    & "
    MAX_WIDTH = 80

    lines: List[str] = []
    prefix = FIRST_PREFIX
    current_vars: List[str] = []
    current_len = len(prefix)

    for var in vars_list:
        # Each var takes len(var) + 2 for ", " separator (except possibly last)
        additional = len(var) + (2 if current_vars else 0)
        if current_vars and current_len + additional > MAX_WIDTH:
            # Flush current line with continuation
            lines.append(prefix + ", ".join(current_vars) + ", &")
            prefix = CONT_PREFIX
            current_vars = [var]
            current_len = len(CONT_PREFIX) + len(var)
        else:
            current_vars.append(var)
            current_len += additional

    # Emit the last line (no trailing continuation — caller adds if needed)
    if current_vars:
        lines.append(prefix + ", ".join(current_vars))

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

    lines = [_HEADER.rstrip()]

    # Normal vars always end with ", &" because opt vars follow (inside the #:if block)
    nl = _format_namelist(normal)
    lines.append(nl + ", &")
    lines.append("#:if not MFC_CASE_OPTIMIZATION")

    # Opt vars: each line except the last ends with ", &"
    for i, var in enumerate(opt):
        is_last = (i == len(opt) - 1)
        if is_last:
            lines.append(f"    & {var}")
        else:
            lines.append(f"    & {var}, &")

    lines.append("#:endif")

    return "\n".join(lines) + "\n"


def generate_decls_fpp(target: str) -> str:
    """Generate simple scalar Fortran variable declarations for a target."""
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
        lines.append(f"{fortran_type_decl(param)} :: {name}")
    return "\n".join(lines) + "\n"


def get_generated_files(include_dir: Path) -> List[Tuple[Path, str]]:
    """Return (output_path, content) for all six generated .fpp files."""
    return [
        (include_dir / "generated_namelist_pre.fpp",  generate_namelist_fpp("pre")),
        (include_dir / "generated_namelist_sim.fpp",  generate_namelist_fpp("sim")),
        (include_dir / "generated_namelist_post.fpp", generate_namelist_fpp("post")),
        (include_dir / "generated_decls_pre.fpp",     generate_decls_fpp("pre")),
        (include_dir / "generated_decls_sim.fpp",     generate_decls_fpp("sim")),
        (include_dir / "generated_decls_post.fpp",    generate_decls_fpp("post")),
    ]
