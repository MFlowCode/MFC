"""Optional tuning overrides; the source files supply the defaults except FUSION_STAGING_BYTES."""

import json
import os
from pathlib import Path

from .state import ARG

# Default fused-kernel staging buffer; case.fpp hands it to m_fused_rhs_openacc.fpp.
FUSION_STAGING_BYTES = 32 * 1024


def profile() -> str:
    """Path of the tuning profile in MFC_TUNING_FILE; empty selects the source defaults."""
    return os.environ.get("MFC_TUNING_FILE", "")


def values(case):
    # A build without --fused-rhs or --fused-bubble-source keeps the source defaults.
    if not (ARG("fused_rhs", False) or ARG("fused_bubble_source", False)):
        return {}
    settings = {}
    profile_path = profile()
    if profile_path:
        data = json.loads(Path(profile_path).read_text())
        settings.update(data.get("values", {}))
    for name, value in os.environ.items():
        if name.startswith("MFC_FUSION_") or name in ("MFC_BUBBLE_VECTOR_LENGTH", "MFC_STATE_ALIGNMENT_BYTES"):
            settings[name] = [int(x) for x in value.split(",")] if name.endswith("CHEM_LANES") else int(value)
    for name, value in settings.items():
        if not (name.startswith("MFC_FUSION_") or name in ("MFC_BUBBLE_VECTOR_LENGTH", "MFC_STATE_ALIGNMENT_BYTES")):
            raise ValueError(f"Unknown tuning setting: {name}")
        if not all(type(x) is int and x > 0 for x in (value if isinstance(value, list) else [value])):
            raise ValueError(f"Invalid tuning setting: {name}={value}")
    return settings


def fusion_staging_bytes(case) -> int:
    return values(case).get("MFC_FUSION_STAGING_BYTES", FUSION_STAGING_BYTES)


def definitions(case):
    return [f"{name}={tuple(value) if isinstance(value, list) else value}" for name, value in sorted(values(case).items())]


# Solver modules that MFC_REGISTER_CAPS can name. The simulation target also
# compiles src/common's solver modules (m_eos, m_variables_conversion), so both
# directories are valid.
_SOURCE_DIRECTORIES = (
    Path(__file__).resolve().parents[2] / "src" / "simulation",
    Path(__file__).resolve().parents[2] / "src" / "common",
)


def register_caps(case=None):
    """Register caps that MFC_REGISTER_CAPS sets per solver module, as module=cap pairs."""
    caps = {}
    for entry in os.environ.get("MFC_REGISTER_CAPS", "").split(","):
        module, separator, value = entry.strip().partition("=")
        if not separator:
            if module:
                raise ValueError(f"Invalid register cap {module} in MFC_REGISTER_CAPS, expected module=cap")
            continue
        if module in caps:
            raise ValueError(f"Duplicate register cap for {module} in MFC_REGISTER_CAPS")
        if not module.isidentifier() or not any((directory / f"{module}.fpp").is_file() for directory in _SOURCE_DIRECTORIES):
            raise ValueError(f"MFC_REGISTER_CAPS names unknown solver module {module}")
        if not value.isdecimal() or not 0 < int(value) < 256:
            raise ValueError(f"Invalid register cap {value} for {module} in MFC_REGISTER_CAPS")
        caps[module] = int(value)
    if caps and case is not None:
        _validate_register_caps(caps, case)
    return caps


def _validate_register_caps(caps, case):
    """Refuse caps nvlink will not link and warn about caps that cannot bind."""
    if case.params.get("chemistry", "F") == "T" and any(caps.get(m) for m in ("m_riemann_state", "m_riemann_solver_hllc", "m_riemann_solver_hll")):
        raise ValueError(
            "MFC_REGISTER_CAPS cannot cap the Riemann solvers on chemistry cases: the solver calls"
            " s_compute_chemistry_average_state (254 registers, m_riemann_state), which stays above"
            " any cap that still buys occupancy, and viscous chemistry adds the generated mixture"
            " viscosity from m_thermochem (228 registers). nvlink refuses the link."
        )
    if int(case.params.get("model_eqns", 2)) == 3 and "m_variables_conversion" in caps:
        raise ValueError(
            "MFC_REGISTER_CAPS=m_variables_conversion on the 6-equation model: the solver executes" " s_phase_pressure_on_isentrope, whose spill under a cap is unmeasured. Drop the entry."
        )
    if (
        ARG("fused_rhs", False)
        and case.params.get("riemann_solver", 2) in (1, 2)
        and case.params.get("mhd", "F") != "T"
        and case.params.get("hypoelasticity", "F") != "T"
        and case.params.get("viscous", "F") != "T"
    ):
        print(
            "[MFC] Warning: MFC_REGISTER_CAPS is set, but the fused lane carries this case's Riemann"
            " solve, so the solver caps do not bind here. Harmless: only the shared conversion"
            " kernels are rescheduled."
        )
