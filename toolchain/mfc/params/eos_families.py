"""The equation-of-state family registry: one entry per family, everything mechanical derives from it.

See docs/superpowers/specs/2026-09-12-eos-family-registry-design.md. `required`/`optional` are
ordered (parameter suffix, Doxygen math symbol) pairs; the full parameter name is
`fluid_pp(i)%<prefix>_<suffix>`. `eos_coeffs` maps a state-dependent family's `eos_coeffs` fields
(src/common/m_derived_types.fpp) to where their value comes from: a parameter (`Param`), a fixed
Fortran literal (`FortranLiteral`), or a value a Fortran helper computes (`Computed`) — not every field is
a copy of a parameter, so a flat name->name dict cannot express all three.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Union


@dataclass(frozen=True)
class Param:
    """eos_coeffs field copied from fluid_pp(i)%<prefix>_<suffix>."""

    suffix: str


@dataclass(frozen=True)
class FortranLiteral:
    """eos_coeffs field set to a fixed Fortran literal; no parameter backs it (e.g. JWL's gruneisen_a)."""

    fortran: str


@dataclass(frozen=True)
class Computed:
    """eos_coeffs field derived by a Fortran helper, not copied or literal (e.g. mu_max)."""

    fortran_call: str


EosCoeffSource = Union[Param, FortranLiteral, Computed]


@dataclass(frozen=True)
class EosFamily:
    value: int  # fluid_pp(i)%eos selector
    suffix: str  # Fortran constant is eos_<suffix>
    label: str  # _EOS_VALUE_LABELS display text
    prefix: str | None  # parameter prefix (e.g. "mg"); None for stiffened_gas/ideal_gas
    required: tuple[tuple[str, str], ...] = ()
    optional: tuple[tuple[str, str], ...] = ()
    state_dependent: bool = False
    isentropic_reference: bool = False  # reference curve is itself an isentrope (JWL, Vinet)
    coefficients_fn: str | None = None  # name of the mirror function in toolchain/mfc/eos.py
    eos_coeffs: dict[str, EosCoeffSource] = field(default_factory=dict)


EOS_FAMILIES: tuple[EosFamily, ...] = (
    EosFamily(value=1, suffix="stiffened_gas", label="stiffened-gas", prefix=None),
    EosFamily(value=2, suffix="ideal_gas", label="ideal-gas", prefix=None),
    EosFamily(
        value=3,
        suffix="mie_gruneisen",
        label="Mie-Gruneisen",
        prefix="mg",
        required=(
            ("rho0", r"\f$\rho_{0,k}\f$"),
            ("c0", r"\f$c_{0,k}\f$"),
            ("s", r"\f$s_k\f$"),
            ("gruneisen", r"\f$\Gamma_{G,k}\f$"),
        ),
        optional=(
            ("gruneisen_a", r"\f$a_k\f$"),
            ("t0", r"\f$T_{0,k}\f$"),
            ("s2", r"\f$s_{2,k}\f$"),
            ("s3", r"\f$s_{3,k}\f$"),
        ),
        state_dependent=True,
        coefficients_fn="eos_coefficients",
        eos_coeffs={
            "c0": Param("c0"),
            "s": Param("s"),
            "s2": Param("s2"),
            "s3": Param("s3"),
            "rho0": Param("rho0"),
            "t0": Param("t0"),
            "gruneisen0": Param("gruneisen"),
            "gruneisen_a": Param("gruneisen_a"),
            "mu_max": Computed("f_hugoniot_compression_limit(mg_c0, mg_s, mg_s2, mg_s3)"),
        },
    ),
    EosFamily(
        value=4,
        suffix="jwl",
        label="JWL",
        prefix="jwl",
        required=(
            ("a", r"\f$A_k\f$"),
            ("b", r"\f$B_k\f$"),
            ("r1", r"\f$R_{1,k}\f$"),
            ("r2", r"\f$R_{2,k}\f$"),
            ("omega", r"\f$\omega_k\f$"),
            ("rho0", r"\f$\rho_{0,k}\f$"),
        ),
        optional=(("t0", r"\f$T_{0,k}\f$"),),
        state_dependent=True,
        isentropic_reference=True,
        coefficients_fn="jwl_coefficients",
        eos_coeffs={
            "a": Param("a"),
            "b": Param("b"),
            "r1": Param("r1"),
            "r2": Param("r2"),
            "rho0": Param("rho0"),
            "t0": Param("t0"),
            "gruneisen0": Param("omega"),
            "gruneisen_a": FortranLiteral("0._wp"),  # JWL has no gruneisen_a parameter
        },
    ),
    EosFamily(
        value=5,
        suffix="vinet",
        label="Vinet",
        prefix="vinet",
        required=(
            ("k0", r"\f$K_{0,k}\f$"),
            ("k0p", r"\f$K'_{0,k}\f$"),
            ("rho0", r"\f$\rho_{0,k}\f$"),
            ("gruneisen", r"\f$\Gamma_{G,k}\f$"),
        ),
        optional=(
            ("gruneisen_a", r"\f$a_k\f$"),
            ("t0", r"\f$T_{0,k}\f$"),
        ),
        state_dependent=True,
        isentropic_reference=True,
        coefficients_fn="vinet_coefficients",
        eos_coeffs={
            "k0": Param("k0"),
            "k0p": Param("k0p"),
            "rho0": Param("rho0"),
            "t0": Param("t0"),
            "gruneisen0": Param("gruneisen"),
            "gruneisen_a": Param("gruneisen_a"),
        },
    ),
)


# The `case default` arm of s_initialize_eos_module: what a fluid whose family is not
# state-dependent gets for each family-dispatched field. Not a family property, so it cannot
# live on an EosFamily. dflt_real marks "unused"; gruneisen_a must be a live zero because
# s_reference_curve reads it unconditionally.
EOS_COEFF_DEFAULT = "dflt_real"
EOS_COEFF_DEFAULTS: dict[str, str] = {"gruneisen_a": "0._wp"}
