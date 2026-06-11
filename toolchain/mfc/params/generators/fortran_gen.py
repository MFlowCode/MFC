"""Fortran parameter code generator — namelist and scalar decl fragments per target."""

from pathlib import Path
from typing import List, Tuple

from ..definitions import CASE_OPT_PARAMS, FORTRAN_ARRAY_DIMS, NAMELIST_VARS, TYPED_DECLS  # noqa: F401 - triggers registry population
from ..registry import REGISTRY
from ..schema import ParamDef, ParamType

TARGETS = [("pre", "pre_process"), ("sim", "simulation"), ("post", "post_process")]
TARGET_FROM_DIR = {full: short for short, full in TARGETS}

_HEADER = "! AUTO-GENERATED - do not edit directly. Regenerate: cmake reconfigure\n!\n"
_MAX_LINE = 130
_FIRST_PREFIX = "namelist /user_inputs/ "
_CONT_PREFIX = "    & "
_CONT2_PREFIX = "        & "  # inside #:if block
_DECL_COL = 24  # '::' column for scalars, matches ffmt alignment
_ARRAY_DECL_COL = 36  # '::' column for array decls

_FORTRAN_TYPES = {
    ParamType.INT: "integer",
    ParamType.ANALYTIC_INT: "integer",
    ParamType.LOG: "logical",
}

_REAL_TYPES = (ParamType.REAL, ParamType.ANALYTIC_REAL)

_VALID_TARGETS = ("pre", "sim", "post")


def _check_target(target: str) -> None:
    if target not in _VALID_TARGETS:
        raise ValueError(f"Unknown target {target!r}; expected one of {_VALID_TARGETS}")


def get_namelist_var(name: str) -> str:
    """Return the Fortran namelist root for a parameter name."""
    if "(" in name and "%" in name:
        return name.split("(", 1)[0]
    if "%" in name:
        return name.split("%", 1)[0]
    return name


def fortran_type_decl(param: ParamDef) -> str:
    if param.param_type == ParamType.STR:
        return f"character(LEN={param.str_len})"
    if param.param_type in _REAL_TYPES:
        return f"real({'stp' if param.storage_precision else 'wp'})"
    return _FORTRAN_TYPES[param.param_type]


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
    _check_target(target)
    all_vars = _vars_for_target(target)

    if target != "sim":
        return _HEADER + _format_namelist(all_vars) + "\n"

    normal = [v for v in all_vars if v not in CASE_OPT_PARAMS]
    opt = sorted(v for v in CASE_OPT_PARAMS if v in NAMELIST_VARS and "sim" in NAMELIST_VARS[v])
    nl_lines = _pack_namelist(normal, _FIRST_PREFIX, _CONT_PREFIX, _MAX_LINE)
    if opt and nl_lines:
        opt_lines = _pack_namelist(opt, _CONT_PREFIX, _CONT2_PREFIX, _MAX_LINE)
        nl_with_cont = nl_lines[:]
        nl_with_cont[-1] += ", &"
        parts = [_HEADER.rstrip(), "#:if MFC_CASE_OPTIMIZATION"] + nl_lines + ["#:else"] + nl_with_cont + opt_lines + ["#:endif"]
    else:
        parts = [_HEADER.rstrip()] + nl_lines
    return "\n".join(parts) + "\n"


def generate_decls_fpp(target: str) -> str:
    """Return Fortran declarations (scalars + known arrays) for a target."""
    _check_target(target)
    lines = [_HEADER.rstrip()]
    for name in _vars_for_target(target):
        if not _is_simple_scalar(name):
            continue
        if target == "sim" and name in CASE_OPT_PARAMS:
            continue
        # TYPED_DECLS handles these; skip here to avoid double-emission.
        if name in TYPED_DECLS:
            continue
        if name in FORTRAN_ARRAY_DIMS:
            member = REGISTRY.all_params.get(f"{name}(1)")
            if member is None:
                raise ValueError(f"FORTRAN_ARRAY_DIMS[{name!r}] has no {name}(1) in the registry. " "Register at least one indexed variant (e.g. _r(f'{name}(1)', ...)).")
            ftype = fortran_type_decl(member)
            dim = FORTRAN_ARRAY_DIMS[name]
            lines.append(f"{(ftype + ', dimension(' + dim + ')').ljust(_ARRAY_DECL_COL)}:: {name}")
            continue
        param = REGISTRY.all_params.get(name)
        if param is None:
            continue
        if any(k.startswith(f"{name}(") for k in REGISTRY.all_params):
            raise ValueError(f"{name!r} has indexed variants (e.g. {name}(1)) but is missing from " "FORTRAN_ARRAY_DIMS. Add it there with its Fortran dimension expression.")
        lines.append(f"{fortran_type_decl(param).ljust(_DECL_COL)}:: {name}")
    for name, (ftype, dim, gpu, desc) in TYPED_DECLS.items():
        if name not in NAMELIST_VARS or target not in NAMELIST_VARS[name]:
            continue
        decl = f"{ftype}, dimension({dim})" if dim else ftype
        padded = decl.ljust(_ARRAY_DECL_COL)
        if not padded.endswith(" "):
            padded += " "
        doc = f" !< {desc}" if desc else ""
        lines.append(f"{padded}:: {name}{doc}")
        if gpu and target == "sim":
            lines.append(f"$:GPU_DECLARE(create='[{name}]')")
    return "\n".join(lines) + "\n"


def generate_constants_fpp() -> str:
    """Named integer constants for enumerated parameters; identical for all targets."""
    from ..definitions import CONSTRAINTS

    lines = [_HEADER.rstrip()]
    for param in sorted(CONSTRAINTS):
        names = CONSTRAINTS[param].get("names")
        if not names:
            continue
        for name, value in sorted(names.items(), key=lambda kv: kv[1]):
            lines.append(f"integer, parameter :: {param}_{name} = {value}")
    return "\n".join(lines) + "\n"


# case.py-computed extras that are not CASE_OPT_PARAMS but appear in the
# case-optimization declaration block (injected as Fypp #:set variables
# by toolchain/mfc/case.py, not from the namelist registry).
# Order matches the reference manual block.
CASE_OPT_EXTRA_LINES = [
    ("num_dims", "integer", "Number of spatial dimensions"),
    ("num_vels", "integer", "Number of velocity components (different from num_dims for mhd)"),
    ("weno_polyn", "integer", "Degree of the WENO polynomials"),
    ("muscl_polyn", "integer", "Degree of the MUSCL polynomials"),
    ("weno_num_stencils", "integer", "Number of stencils for WENO reconstruction"),
    ("wenojs", "logical", "WENO-JS (default)"),
]

_CASE_OPT_DECL_COL = 24  # '::' alignment for case-opt declarations


def generate_case_opt_decls_fpp() -> str:
    """Return the case-optimization declaration block for src/simulation/m_global_parameters.fpp.

    Emits one #:if MFC_CASE_OPTIMIZATION / #:else / #:endif block containing:
    - The CASE_OPT_EXTRA_LINES (case.py-computed: num_dims, num_vels, weno_polyn,
      muscl_polyn, weno_num_stencils, wenojs) — stored as a literal list because
      they are not namelist-registry parameters.
    - Every CASE_OPT_PARAMS entry (excluding nb, which lives in a separate block)
      driven by the registry for its Fortran type.
    """
    params_to_emit = sorted(CASE_OPT_PARAMS - {"nb"})

    if_lines = []
    else_lines = []

    # Emit extras first (case.py-computed, not in registry)
    for name, ftype, desc in CASE_OPT_EXTRA_LINES:
        if ftype == "logical":
            rhs = f"(${{{name}}}$ /= 0)"
            if_lines.append(f"    {ftype}, parameter :: {name} = {rhs}  !< {desc}")
        else:
            if_lines.append(f"    {ftype}, parameter :: {name} = ${{{name}}}$  !< {desc}")
        else_lines.append(f"    {ftype.ljust(_CASE_OPT_DECL_COL)}:: {name}")

    # Emit registry-driven CASE_OPT_PARAMS
    for name in params_to_emit:
        pdef = REGISTRY.get_param_def(name)
        if pdef is None:
            raise ValueError(f"CASE_OPT_PARAMS entry {name!r} not found in registry")
        desc = pdef.description or ""
        if pdef.param_type in _REAL_TYPES:
            ftype = "real(wp)"
            if_lines.append(f"    {ftype}, parameter :: {name} = ${{{name}}}$  !< {desc}")
            else_lines.append(f"    {ftype.ljust(_CASE_OPT_DECL_COL)}:: {name}")
        elif pdef.param_type == ParamType.LOG:
            ftype = "logical"
            rhs = f"(${{{name}}}$ /= 0)"
            if_lines.append(f"    {ftype}, parameter :: {name} = {rhs}  !< {desc}")
            else_lines.append(f"    {ftype.ljust(_CASE_OPT_DECL_COL)}:: {name}")
        else:
            ftype = "integer"
            if_lines.append(f"    {ftype}, parameter :: {name} = ${{{name}}}$  !< {desc}")
            else_lines.append(f"    {ftype.ljust(_CASE_OPT_DECL_COL)}:: {name}")

    parts = [_HEADER.rstrip(), "#:if MFC_CASE_OPTIMIZATION"]
    parts.extend(if_lines)
    parts.append("#:else")
    parts.extend(else_lines)
    parts.append("#:endif")
    return "\n".join(parts) + "\n"


def resolve_namelist_content(fpp_path: Path) -> str:
    """Return the namelist content for an fpp file.

    If the file delegates to a generated include, returns the generated content.
    Otherwise returns the file's raw text.
    """
    text = fpp_path.read_text()
    if "#:include 'generated_namelist.fpp'" not in text:
        return text
    short = TARGET_FROM_DIR.get(fpp_path.parent.name)
    if short is None:
        raise ValueError(f"Cannot determine MFC target from path: {fpp_path}")
    return generate_namelist_fpp(short)


def get_generated_files(build_dir: Path) -> List[Tuple[Path, str]]:
    """Return (path, content) for all 10 generated .fpp files under build_dir.

    Paths match the cmake include directory structure:
      build_dir/include/{full_target}/generated_{namelist,decls,constants}.fpp
    The simulation target also gets generated_case_opt_decls.fpp.
    """
    result = []
    for short, full in TARGETS:
        inc = build_dir / "include" / full
        result.append((inc / "generated_namelist.fpp", generate_namelist_fpp(short)))
        result.append((inc / "generated_decls.fpp", generate_decls_fpp(short)))
        result.append((inc / "generated_constants.fpp", generate_constants_fpp()))
    sim_inc = build_dir / "include" / "simulation"
    result.append((sim_inc / "generated_case_opt_decls.fpp", generate_case_opt_decls_fpp()))
    return result
