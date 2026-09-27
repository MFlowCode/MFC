"""Resolve external analytic-EOS material files into ordinary case parameters."""

import math
import os
import re
import typing
from pathlib import Path

import yaml

from . import eos
from .common import MFCException
from .params.eos_families import EOS_FAMILIES

_MATERIAL_KEY = re.compile(r"^fluid_pp\(([1-9][0-9]*)\)%material_file$")
_EOS = {family.suffix: family for family in EOS_FAMILIES if family.state_dependent}
_LAYOUT = {"material": {"name", "eos"}, "provenance": {"citation", "release_status"}}
_PHASE_KEYS = ("cv", "qv", "q")  # per-fluid, not family-prefixed; q is JWL's detonation energy Q


def _material_path(filename: str, case_dir: str) -> Path:
    path, public_dir = Path(str(filename)).expanduser(), os.environ.get("MFC_PUBLIC_MATERIAL_DIR")
    candidates = [path] if path.is_absolute() else [path, Path(case_dir) / path] + ([Path(public_dir).expanduser() / path] if public_dir else [])
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise MFCException(f"Material file '{filename}' not found. Searched: {', '.join(map(str, candidates))}.")


def _number(value) -> float:
    if isinstance(value, bool) or not math.isfinite(number := float(value)):
        raise ValueError
    return number


def _read_material(path: Path) -> typing.Tuple[str, dict, typing.Optional[float]]:
    """Family, fluid_pp coefficients and JWL Q of one material file; any malformed field is an error."""
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        entry = data["material"]["eos"]
        family, spec = entry["family"], _EOS[entry["family"]]
        texts = (data["material"]["name"], *data["provenance"].values())
        if {key: set(value) for key, value in data.items()} != _LAYOUT or set(entry) != {"family", "parameters"} or not all(isinstance(t, str) and t.strip() for t in texts):
            raise ValueError
        params = {str(name).lower(): _number(value) for name, value in entry["parameters"].items()}
    except (OSError, UnicodeError, yaml.YAMLError, TypeError, KeyError, ValueError, AttributeError) as exc:
        raise MFCException(
            f"Material file '{path}' requires material.name, material.eos.family ({', '.join(_EOS)}), finite numeric material.eos.parameters and a nonempty provenance citation and release_status"
        ) from exc
    if len(params) != len(entry["parameters"]):
        raise MFCException(f"Material file '{path}' names a parameter twice")
    # Names the family lacks, and Q outside JWL, reach the schema as unknown parameters.
    coefficients = {name if name in _PHASE_KEYS else f"{spec.prefix}_{name}": value for name, value in params.items()}
    q = coefficients.pop("q", None) if family == "jwl" else None
    return family, coefficients, q


def _reactant_qv(params: dict, q: float, rho0: float) -> float:
    """Reactant qv that puts the unreacted state (rho0, p = 0) at energy Q on the products' JWL scale."""
    eos_id = params.get("fluid_pp(1)%eos")
    family = next((f for f in _EOS.values() if eos_id in (f.suffix, f.value)), None)
    try:  # a stiffened gas has a constant Pi, an ideal gas none
        pi_inf = eos.family_coefficients(family, params.get, 1, rho0)[1] if family else params.get("fluid_pp(1)%pi_inf") or 0.0
        return q - pi_inf / rho0
    except (TypeError, ValueError, OverflowError, ZeroDivisionError) as exc:
        raise MFCException(f"reactive_burn with JWL Q cannot evaluate the reactant EOS at rho0 = {rho0}: {exc}") from exc


def resolve_materials(params: dict, case_dir: str) -> dict:
    """Expand fluid_pp(i)%material_file before schema validation or namelist generation."""
    resolved, products_q = dict(params), None
    for directive, filename in params.items():
        if (match := _MATERIAL_KEY.fullmatch(directive)) is None:
            continue
        phase = match.group(1)
        family, coefficients, q = _read_material(_material_path(filename, case_dir))
        # Q is the fit's detonation energy: a blast started from products already carries it, a burn releases it.
        if q is not None and params.get("reactive_burn", "F") == "T":
            if phase != "2":
                raise MFCException(f"{directive} contains Q, which belongs to the products (fluid 2) of a reactive burn")
            products_q = (q, coefficients["jwl_rho0"])
        for name, value in {"eos": family, **coefficients}.items():
            key = f"fluid_pp({phase})%{name}"
            if key in resolved and (name != "eos" or resolved[key] not in (family, _EOS[family].value)):
                raise MFCException(f"{directive} conflicts with {key}")
            resolved[key] = value
        del resolved[directive]
    if products_q is not None:
        if any(f"fluid_pp({k})%qv" in resolved for k in (1, 2)):
            raise MFCException("reactive_burn with JWL Q sets both reactant and product qv; remove them from the case and material files")
        resolved["fluid_pp(1)%qv"], resolved["fluid_pp(2)%qv"] = _reactant_qv(resolved, *products_q), 0.0
    return resolved
