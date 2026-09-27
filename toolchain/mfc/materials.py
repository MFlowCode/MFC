"""Resolve external analytic-EOS material files into ordinary case parameters."""

import math
import os
import re
from pathlib import Path

import yaml

from . import eos
from .common import MFCException
from .params.eos_families import EOS_FAMILIES
from .printer import cons

_MATERIAL_KEY = re.compile(r"^fluid_pp\(([1-9][0-9]*)\)%material_file$")
_EOS = {family.suffix: family for family in EOS_FAMILIES if family.state_dependent}


def _material_path(filename: str, case_dir: str) -> Path:
    path = Path(filename).expanduser()
    candidates = [path] if path.is_absolute() else [path, Path(case_dir) / path]
    public_dir = os.environ.get("MFC_PUBLIC_MATERIAL_DIR")
    if public_dir and not path.is_absolute():
        candidates.append(Path(public_dir).expanduser() / path)
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise MFCException(f"Material file '{filename}' not found. Searched: {', '.join(map(str, candidates))}.")


def _read_material(path: Path) -> tuple[str, dict, float | None]:
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise MFCException(f"Cannot read material file '{path}': {exc}") from exc
    if not isinstance(data, dict) or set(data) != {"material", "provenance"}:
        raise MFCException(f"Material file '{path}' requires material and provenance mappings")
    material, provenance = data["material"], data["provenance"]
    if not isinstance(material, dict) or set(material) != {"name", "eos"} or not isinstance(material["name"], str) or not material["name"].strip():
        raise MFCException(f"Material file '{path}' requires material.name and material.eos")
    if not isinstance(provenance, dict) or set(provenance) != {"citation", "release_status"} or any(not isinstance(v, str) or not v.strip() for v in provenance.values()):
        raise MFCException(f"Material file '{path}' requires nonempty provenance.citation and provenance.release_status")
    eos = material["eos"]
    if not isinstance(eos, dict) or set(eos) != {"family", "parameters"} or not isinstance(eos["family"], str) or eos["family"] not in _EOS:
        raise MFCException(f"Material file '{path}' requires a supported analytic eos.family and eos.parameters")
    family, params = eos["family"], eos["parameters"]
    spec = _EOS[family]
    if not isinstance(params, dict):
        raise MFCException(f"Material file '{path}' eos.parameters must be a mapping")
    normalized = {}
    for name, value in params.items():
        key = str(name).lower()
        if key in normalized:
            raise MFCException(f"Material file '{path}' repeats EOS parameter '{key}'")
        normalized[key] = value
    params = normalized
    required = {name for name, _ in spec.required}
    allowed = required | {name for name, _ in spec.optional} | {"cv", "qv"}
    if family == "jwl":
        allowed.add("q")
    if not required <= params.keys() or params.keys() - allowed:
        raise MFCException(f"Material file '{path}' {family} parameters require {sorted(required)}; optional: {sorted(allowed - required)}")
    coefficients = {}
    for name, raw_value in params.items():
        try:
            value = float(raw_value)
            if isinstance(raw_value, bool) or not math.isfinite(value):
                raise ValueError
        except (TypeError, ValueError) as exc:
            raise MFCException(f"Material file '{path}' requires finite numeric EOS parameters") from exc
        coefficients[name if name in ("cv", "qv", "q") else f"{spec.prefix}_{name}"] = value
    q = coefficients.pop("q", None)
    if q is not None and q <= 0:
        raise MFCException(f"Material file '{path}' requires a positive JWL Q")
    return family, coefficients, q


def _reactant_qv(params: dict, q: float, rho0: float) -> float:
    """Reactant qv that puts the unreacted state (rho0, p = 0) at energy Q on the products' JWL scale."""
    family = _EOS.get(params.get("fluid_pp(1)%eos"))
    try:
        if family is None:  # stiffened or ideal gas: a constant Pi, absent for an ideal gas
            pi_inf = float(params.get("fluid_pp(1)%pi_inf") or 0.0)
        else:
            pi_inf = eos.family_coefficients(family, params.get, 1, rho0)[1]
    except (TypeError, ValueError, OverflowError, ZeroDivisionError) as exc:
        raise MFCException(f"reactive_burn with JWL Q cannot evaluate the reactant EOS at rho0 = {rho0}: {exc}") from exc
    return q - pi_inf / rho0


def resolve_materials(params: dict, case_dir: str) -> dict:
    """Expand fluid_pp(i)%material_file before schema validation or namelist generation."""
    resolved = dict(params)
    reactive = params.get("reactive_burn", "F") == "T"
    products_q = None
    for directive, filename in params.items():
        match = _MATERIAL_KEY.fullmatch(directive)
        if match is None:
            continue
        if not isinstance(filename, str) or not filename.strip():
            raise MFCException(f"{directive} requires a material YAML path")
        phase = match.group(1)
        if isinstance(params.get("num_fluids"), int) and int(phase) > params["num_fluids"]:
            raise MFCException(f"{directive} refers to a fluid beyond num_fluids")
        family, coefficients, q = _read_material(_material_path(filename, case_dir))
        if q is not None:
            if not reactive:
                cons.print("[yellow]Warning:[/yellow] material Q is metadata without reactive_burn; use qv to set a runtime energy offset.")
            elif phase != "2":
                raise MFCException(f"{directive} contains Q, which belongs to the products (fluid 2) of a reactive burn")
            elif "qv" in coefficients or any(f"fluid_pp({k})%qv" in params for k in (1, 2)):
                raise MFCException(f"{directive} contains Q, which sets both reactant and product qv; remove them from the case and material files")
            else:
                products_q = (q, coefficients["jwl_rho0"])
        for name, value in {"eos": family, **coefficients}.items():
            key = f"fluid_pp({phase})%{name}"
            if key in resolved and (name != "eos" or resolved[key] not in (family, _EOS[family].value)):
                raise MFCException(f"{directive} conflicts with {key}")
            resolved[key] = value
        del resolved[directive]
    if products_q is not None:
        resolved["fluid_pp(1)%qv"] = _reactant_qv(resolved, *products_q)
        resolved["fluid_pp(2)%qv"] = 0.0
    return resolved
