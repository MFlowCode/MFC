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


# Struct roots in NAMELIST_VARS whose member-level broadcasts are irregular
# (per-target member subsets, grouped array members, nested loops, etc.).
# These are kept in the manual residue of m_mpi_proxy.fpp.
_STRUCT_ROOTS = frozenset({"bc_x", "bc_y", "bc_z", "x_domain", "y_domain", "z_domain", "x_output", "y_output", "z_output"})

# Variables excluded from broadcast generation (derived post-broadcast or non-namelist).
_BCAST_EXCLUDE = frozenset({"muscl_eps"})

# Post-process scalars that are namelist-bound but consumed on rank 0 only (reading/init).
# Broadcasting them would be harmless but changes the existing call set, which we preserve.
_POST_BCAST_EXCLUDE = frozenset({"avg_state", "cfl_target", "igr_order", "num_bc_patches", "recon_type", "sigR"})

# TYPED_DECLS entries kept entirely in the manual residue: their member-broadcast structure
# is irregular (non-uniform subsets, size() arrays, nested loops, complex guards).


def _mpi_type_for(ptype: ParamType) -> str:
    """Return the Fortran MPI type constant for a ParamType."""
    if ptype in _REAL_TYPES:
        return "mpi_p"
    if ptype in (ParamType.INT, ParamType.ANALYTIC_INT):
        return "MPI_INTEGER"
    if ptype == ParamType.LOG:
        return "MPI_LOGICAL"
    if ptype == ParamType.STR:
        return "MPI_CHARACTER"
    raise ValueError(f"No MPI type mapping for {ptype!r}")


def _bcast_scalar(name: str, mpi_type: str, count: str = "1") -> str:
    return f"        call MPI_BCAST({name}, {count}, {mpi_type}, 0, MPI_COMM_WORLD, ierr)"


def _classify_scalar_vars(target: str) -> Tuple[List[str], List[str], List[str], List[str]]:
    """Return (int_vars, log_vars, real_vars, case_opt_vars) for class-(a) scalars.

    case_opt_vars are sim CASE_OPT_PARAMS (wrapped in #:if not MFC_CASE_OPTIMIZATION).
    All four lists contain variable names suitable for a 1-element MPI_BCAST.
    Struct roots, TYPED_DECLS, FORTRAN_ARRAY_DIMS, case_dir, and _BCAST_EXCLUDE are removed.
    """
    int_vars: List[str] = []
    log_vars: List[str] = []
    real_vars: List[str] = []
    case_opt_vars: List[str] = []  # (name, mpi_type) pairs for sim case-opt section

    for name in sorted(NAMELIST_VARS):
        if target not in NAMELIST_VARS[name]:
            continue
        if name in _BCAST_EXCLUDE:
            continue
        if target == "post" and name in _POST_BCAST_EXCLUDE:
            continue
        if name in _STRUCT_ROOTS:
            continue
        if name in TYPED_DECLS:
            continue
        if name in FORTRAN_ARRAY_DIMS:
            continue
        if name == "case_dir":
            continue

        pdef = REGISTRY.all_params.get(name)
        if pdef is None:
            continue  # no registry entry — skip (e.g. 'G' in post is rank-0-only)

        if target == "sim" and name in CASE_OPT_PARAMS:
            case_opt_vars.append(name)
            continue

        mpi_type = _mpi_type_for(pdef.param_type)
        if mpi_type == "MPI_INTEGER":
            int_vars.append(name)
        elif mpi_type == "MPI_LOGICAL":
            log_vars.append(name)
        else:
            real_vars.append(name)

    return sorted(int_vars), sorted(log_vars), sorted(real_vars), sorted(case_opt_vars)


def _emit_bcast_group(lines: List[str], vars_list: List[str], mpi_type: str) -> None:
    """Append one MPI_BCAST per variable in vars_list to lines."""
    for name in vars_list:
        lines.append(_bcast_scalar(name, mpi_type))


def _emit_fluid_pp(lines: List[str], target: str) -> None:
    """Emit the fluid_pp(i) member-loop broadcast block.

    Members broadcast: all REAL registry members of fluid_pp (gamma, pi_inf, G, cv, qv,
    qvp) derived from physical_parameters.  Sim additionally: Re(1) with count=2.
    mul0/ss/pv/gamma_v/M_v/mu_v/k_v/cp_v/D_v were removed from the Fortran type by
    upstream #1085/#1093 and are no longer registered.
    """
    # Walk the registry for fluid_pp REAL members (Re handled separately; exclude).
    fp_real_members = sorted(k.split("%", 1)[1] for k in REGISTRY.all_params if k.startswith("fluid_pp(1)%") and not k.startswith("fluid_pp(1)%Re("))
    lines.append("        do i = 1, num_fluids_max")
    for mem in fp_real_members:
        lines.append(f"            call MPI_BCAST(fluid_pp(i)%{mem}, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)")
    if target == "sim":
        lines.append("            call MPI_BCAST(fluid_pp(i)%Re(1), 2, mpi_p, 0, MPI_COMM_WORLD, ierr)")
    lines.append("        end do")


def _emit_bub_pp(lines: List[str]) -> None:
    """Emit the bub_pp member broadcast block (all targets, under bubbles guard).

    All bub_pp members in the registry are REAL — emit them all sorted.
    """
    bub_members = sorted(k.split("%", 1)[1] for k in REGISTRY.all_params if k.startswith("bub_pp%"))
    lines.append("        if (bubbles_euler .or. bubbles_lagrange) then")
    for mem in bub_members:
        lines.append(f"            call MPI_BCAST(bub_pp%{mem}, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)")
    lines.append("        end if")


def _emit_lag_params(lines: List[str]) -> None:
    """Emit the lag_params member broadcast block (sim-only, under bubbles_lagrange guard).

    All registered lag_params members are broadcast.  T0/Thost/c0/rho0/x0 were removed
    from the Fortran type by upstream #1085/#1093 and are no longer in the registry.
    """
    # Walk the registry for lag_params members, split by type.
    lag_log = sorted(k.split("%", 1)[1] for k in REGISTRY.all_params if k.startswith("lag_params%") and REGISTRY.all_params[k].param_type == ParamType.LOG)
    lag_int = sorted(k.split("%", 1)[1] for k in REGISTRY.all_params if k.startswith("lag_params%") and REGISTRY.all_params[k].param_type in (ParamType.INT, ParamType.ANALYTIC_INT))
    lag_real = sorted(k.split("%", 1)[1] for k in REGISTRY.all_params if k.startswith("lag_params%") and REGISTRY.all_params[k].param_type in _REAL_TYPES)
    lines.append("        if (bubbles_lagrange) then")
    for mem in sorted(lag_log):
        lines.append(f"            call MPI_BCAST(lag_params%{mem}, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)")
    for mem in sorted(lag_int):
        lines.append(f"            call MPI_BCAST(lag_params%{mem}, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)")
    for mem in sorted(lag_real):
        lines.append(f"            call MPI_BCAST(lag_params%{mem}, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)")
    lines.append("        end if")


def _emit_chem_params(lines: List[str]) -> None:
    """Emit the chem_params member broadcast block (sim-only, under chemistry guard).

    Broadcasts the registry's LOG and INT chem_params members (the only kinds registered today; extend the type split if REAL members appear).
    """
    chem_members = sorted(k.split("%", 1)[1] for k in REGISTRY.all_params if k.startswith("chem_params%"))
    chem_log = [m for m in chem_members if REGISTRY.all_params[f"chem_params%{m}"].param_type == ParamType.LOG]
    chem_int = [m for m in chem_members if REGISTRY.all_params[f"chem_params%{m}"].param_type == ParamType.INT]
    lines.append("        if (chemistry) then")
    for mem in sorted(chem_log):
        lines.append(f"            call MPI_BCAST(chem_params%{mem}, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)")
    for mem in sorted(chem_int):
        lines.append(f"            call MPI_BCAST(chem_params%{mem}, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)")
    lines.append("        end if")


def _emit_fortran_array_dims(lines: List[str], target: str) -> None:
    """Emit broadcasts for FORTRAN_ARRAY_DIMS entries belonging to this target.

    Each entry is a 1D array; we emit: MPI_BCAST(name(1), dim, mpi_type, ...).
    The element type is determined by the registered member param_type.
    """
    for name in sorted(FORTRAN_ARRAY_DIMS):
        if name not in NAMELIST_VARS or target not in NAMELIST_VARS[name]:
            continue
        dim = FORTRAN_ARRAY_DIMS[name]
        # Determine element type from registry (use the (1) example entry)
        member_key = f"{name}(1)"
        pdef = REGISTRY.all_params.get(member_key)
        if pdef is None:
            raise ValueError(f"No registry entry for {member_key!r} (needed for FORTRAN_ARRAY_DIMS broadcast)")
        mpi_type = _mpi_type_for(pdef.param_type)
        lines.append(f"        call MPI_BCAST({name}(1), {dim}, {mpi_type}, 0, MPI_COMM_WORLD, ierr)")


def generate_bcast_fpp(target: str) -> str:
    """Return the generated MPI broadcast statements for a target as a Fortran string.

    Generates class-(a) namelist scalar broadcasts (sorted, MPI type from registry)
    and class-(b) struct/array broadcasts for cleanly walkable TYPED_DECLS entries.
    The manual residue (bc_x members, patch_ib, patch_icpp, simplex_params, etc.)
    stays in the hand-written m_mpi_proxy.fpp.

    For sim, CASE_OPT_PARAMS scalars are wrapped in ``#:if not MFC_CASE_OPTIMIZATION``.
    The ``chem_wrt_Y`` array broadcast for post is intentionally ADDED here (latent bug fix:
    it is namelist-bound and consumed on all ranks but was previously not broadcast).
    """
    _check_target(target)
    lines: List[str] = [_HEADER.rstrip()]

    # -- case_dir (STR, special len() form) --
    lines.append("        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)")
    lines.append("")

    # -- Class (a): simple scalar broadcasts --
    int_vars, log_vars, real_vars, case_opt_vars = _classify_scalar_vars(target)

    if int_vars:
        lines.append("        ! Integer scalars")
        _emit_bcast_group(lines, int_vars, "MPI_INTEGER")
        lines.append("")

    if log_vars:
        lines.append("        ! Logical scalars")
        _emit_bcast_group(lines, log_vars, "MPI_LOGICAL")
        lines.append("")

    if real_vars:
        lines.append("        ! Real scalars")
        _emit_bcast_group(lines, real_vars, "mpi_p")
        lines.append("")

    # Sim case-optimization guard
    if target == "sim" and case_opt_vars:
        lines.append("        ! Case-optimization scalars (absent when constants are baked in)")
        lines.append("        #:if not MFC_CASE_OPTIMIZATION")
        # Classify by type
        co_int = sorted(v for v in case_opt_vars if _mpi_type_for(REGISTRY.all_params[v].param_type) == "MPI_INTEGER")
        co_log = sorted(v for v in case_opt_vars if _mpi_type_for(REGISTRY.all_params[v].param_type) == "MPI_LOGICAL")
        co_real = sorted(v for v in case_opt_vars if _mpi_type_for(REGISTRY.all_params[v].param_type) == "mpi_p")
        for name in co_int:
            lines.append("    " + _bcast_scalar(name, "MPI_INTEGER"))
        for name in co_log:
            lines.append("    " + _bcast_scalar(name, "MPI_LOGICAL"))
        for name in co_real:
            lines.append("    " + _bcast_scalar(name, "mpi_p"))
        lines.append("        #:endif")
        lines.append("")

    # -- Class (b): FORTRAN_ARRAY_DIMS (simple 1D arrays) --
    array_dim_names = [n for n in FORTRAN_ARRAY_DIMS if n in NAMELIST_VARS and target in NAMELIST_VARS[n]]
    if array_dim_names:
        lines.append("        ! Array broadcasts (dimension from FORTRAN_ARRAY_DIMS)")
        _emit_fortran_array_dims(lines, target)
        lines.append("")

    # -- Class (b): TYPED_DECLS structs (cleanly walkable) --
    if "fluid_pp" in NAMELIST_VARS and target in NAMELIST_VARS["fluid_pp"]:
        lines.append("        ! fluid_pp member loop")
        _emit_fluid_pp(lines, target)
        lines.append("")

    if "bub_pp" in NAMELIST_VARS and target in NAMELIST_VARS["bub_pp"]:
        lines.append("        ! bub_pp members (under bubbles guard)")
        _emit_bub_pp(lines)
        lines.append("")

    if target == "sim":
        if "lag_params" in NAMELIST_VARS and "sim" in NAMELIST_VARS["lag_params"]:
            lines.append("        ! lag_params members (under bubbles_lagrange guard)")
            _emit_lag_params(lines)
            lines.append("")
        if "chem_params" in NAMELIST_VARS and "sim" in NAMELIST_VARS["chem_params"]:
            lines.append("        ! chem_params members (under chemistry guard)")
            _emit_chem_params(lines)
            lines.append("")

    return "\n".join(lines) + "\n"


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
    """Return (path, content) for all 15 generated .fpp files under build_dir.

    Paths match the cmake include directory structure:
      build_dir/include/{full_target}/generated_{namelist,decls,constants,case_opt_decls,bcast}.fpp
    Every target gets generated_case_opt_decls.fpp: real content for simulation,
    a header-only stub for the others (Fypp resolves #:include at parse time, so
    the file must exist for every target even inside a dead conditional).
    Every target gets generated_bcast.fpp with its MPI broadcast statements.
    """
    result = []
    for short, full in TARGETS:
        inc = build_dir / "include" / full
        result.append((inc / "generated_namelist.fpp", generate_namelist_fpp(short)))
        result.append((inc / "generated_decls.fpp", generate_decls_fpp(short)))
        result.append((inc / "generated_constants.fpp", generate_constants_fpp()))
    for short, full in TARGETS:
        inc = build_dir / "include" / full
        content = generate_case_opt_decls_fpp() if short == "sim" else _HEADER + "! (no case-optimization declarations for this target)\n"
        result.append((inc / "generated_case_opt_decls.fpp", content))
    for short, full in TARGETS:
        inc = build_dir / "include" / full
        result.append((inc / "generated_bcast.fpp", generate_bcast_fpp(short)))
    return result
