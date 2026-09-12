"""R1: proves the EOS_FAMILIES registry changed nothing observable.

R5 (below the R1 tests): proves the registry agrees with the parts of the codebase that are
still hand-written by design — physical_parameters, s_initialize_eos_module, eos.py and
s_reference_curve — so a family that validates cleanly but computes with dflt_real coefficients
is caught here rather than at runtime.

Expected values are hard-coded (captured from the tree before the registry existed), not derived
from the registry itself — asserting the registry against itself would be circular.
"""

import re

from .. import eos as eos_module
from ..params import REGISTRY
from ..params.definitions import CONSTRAINTS
from ..params.eos_families import EOS_COEFF_DEFAULTS, EOS_FAMILIES
from ..params.generators.fortran_gen import generate_eos_fpp
from ..params.namelist_parser import get_mfc_root

_EXPECTED_NAMES = {"stiffened_gas": 1, "ideal_gas": 2, "mie_gruneisen": 3, "jwl": 4, "vinet": 5}
_EXPECTED_LABELS = {1: "stiffened-gas", 2: "ideal-gas", 3: "Mie-Gruneisen", 4: "JWL", 5: "Vinet"}
_EXPECTED_CHOICES = [1, 2, 3, 4, 5]

# Captured from `sorted(n for n in REGISTRY.all_params if n.startswith("fluid_pp(1)%"))`
# on the unmodified tree (HEAD 949a131c), before eos_families.py existed.
_EXPECTED_FLUID_PP1_PARAMS = [
    "fluid_pp(1)%G",
    "fluid_pp(1)%K",
    "fluid_pp(1)%Re(1)",
    "fluid_pp(1)%Re(2)",
    "fluid_pp(1)%cv",
    "fluid_pp(1)%eos",
    "fluid_pp(1)%gamma",
    "fluid_pp(1)%hb_m",
    "fluid_pp(1)%jwl_a",
    "fluid_pp(1)%jwl_b",
    "fluid_pp(1)%jwl_omega",
    "fluid_pp(1)%jwl_r1",
    "fluid_pp(1)%jwl_r2",
    "fluid_pp(1)%jwl_rho0",
    "fluid_pp(1)%jwl_t0",
    "fluid_pp(1)%mg_c0",
    "fluid_pp(1)%mg_gruneisen",
    "fluid_pp(1)%mg_gruneisen_a",
    "fluid_pp(1)%mg_rho0",
    "fluid_pp(1)%mg_s",
    "fluid_pp(1)%mg_s2",
    "fluid_pp(1)%mg_s3",
    "fluid_pp(1)%mg_t0",
    "fluid_pp(1)%mu_bulk",
    "fluid_pp(1)%mu_max",
    "fluid_pp(1)%mu_min",
    "fluid_pp(1)%nn",
    "fluid_pp(1)%non_newtonian",
    "fluid_pp(1)%pi_inf",
    "fluid_pp(1)%qv",
    "fluid_pp(1)%qvp",
    "fluid_pp(1)%tau0",
    "fluid_pp(1)%vinet_gruneisen",
    "fluid_pp(1)%vinet_gruneisen_a",
    "fluid_pp(1)%vinet_k0",
    "fluid_pp(1)%vinet_k0p",
    "fluid_pp(1)%vinet_rho0",
    "fluid_pp(1)%vinet_t0",
]


def test_eos_names_unchanged():
    assert CONSTRAINTS["fluid_pp(1)%eos"]["names"] == _EXPECTED_NAMES


def test_eos_value_labels_unchanged():
    assert CONSTRAINTS["fluid_pp(1)%eos"]["value_labels"] == _EXPECTED_LABELS


def test_eos_choices_unchanged():
    assert CONSTRAINTS["fluid_pp(1)%eos"]["choices"] == _EXPECTED_CHOICES


def test_fluid_pp_parameter_set_unchanged():
    names = sorted(n for n in REGISTRY.all_params if n.startswith("fluid_pp(1)%"))
    assert names == _EXPECTED_FLUID_PP1_PARAMS


# R5: agreement tests. Each closes one leg of the web docs/superpowers/specs/2026-09-12-eos-
# family-registry-design.md maps between the registry and the hand-written code it drives.

_M_EOS = (get_mfc_root() / "src" / "common" / "m_eos.fpp").read_text()
_M_DERIVED_TYPES = (get_mfc_root() / "src" / "common" / "m_derived_types.fpp").read_text()

_STATE_DEPENDENT = [f for f in EOS_FAMILIES if f.state_dependent]


def _physical_parameters_fields():
    """Field names declared on `type physical_parameters`, splitting comma-joined `::` lines
    (e.g. `real(wp) :: mg_s2, mg_s3` is two fields, not one)."""
    body = re.search(r"type physical_parameters\b(.*?)end type physical_parameters", _M_DERIVED_TYPES, re.S).group(1)
    fields = set()
    for line in body.splitlines():
        code = line.split("!", 1)[0]
        if "::" not in code:
            continue
        for raw in code.split("::", 1)[1].split(","):
            name = raw.strip()
            if name:
                fields.add(name)
    return fields


def test_every_registry_parameter_has_a_physical_parameters_field():
    fields = _physical_parameters_fields()
    for family in EOS_FAMILIES:
        for suffix, _math in family.required + family.optional:
            name = f"{family.prefix}_{suffix}"
            assert name in fields, f"{name} (family {family.suffix}) has no physical_parameters field"


def _init_module_text():
    """s_initialize_eos_module's body, plus the generated macro bodies it calls (`@:EOS_INIT_COEFFS`
    / `@:EOS_INIT_REFERENCE_STATE`) — the parameters are read via those macros, not literally in
    m_eos.fpp, since R4 moved them into generated_eos.fpp."""
    body = re.search(r"impure subroutine s_initialize_eos_module\(\)(.*?)end subroutine s_initialize_eos_module", _M_EOS, re.S).group(1)
    return body + "\n" + generate_eos_fpp()


def test_every_required_parameter_is_read_by_the_init():
    text = _init_module_text()
    for family in _STATE_DEPENDENT:
        for suffix, _math in family.required:
            name = f"{family.prefix}_{suffix}"
            assert re.search(rf"\b{re.escape(name)}\b", text), f"{name} (family {family.suffix}) is never read in s_initialize_eos_module"


def test_every_coefficients_fn_resolves():
    for family in _STATE_DEPENDENT:
        fn = getattr(eos_module, family.coefficients_fn, None)
        assert callable(fn), f"eos.py has no callable {family.coefficients_fn!r} (family {family.suffix})"


def _reference_curve_text():
    return re.search(r"subroutine s_reference_curve\(.*?\)(.*?)end subroutine s_reference_curve", _M_EOS, re.S).group(1)


def test_every_family_has_a_reference_curve_case():
    text = _reference_curve_text()
    for family in _STATE_DEPENDENT:
        assert f"case (eos_{family.suffix})" in text, f"no case (eos_{family.suffix}) in s_reference_curve"


# Verified by hand against the registry (docs/superpowers/plans/2026-09-12-eos-family-registry.md,
# R1 review). Hard-coded so an edit that garbles a LaTeX symbol while keeping the parameter name
# is caught here rather than slipping through unreviewed.
_EXPECTED_MATH = {
    ("mie_gruneisen", "rho0"): r"\f$\rho_{0,k}\f$",
    ("mie_gruneisen", "c0"): r"\f$c_{0,k}\f$",
    ("mie_gruneisen", "s"): r"\f$s_k\f$",
    ("mie_gruneisen", "gruneisen"): r"\f$\Gamma_{G,k}\f$",
    ("mie_gruneisen", "gruneisen_a"): r"\f$a_k\f$",
    ("mie_gruneisen", "t0"): r"\f$T_{0,k}\f$",
    ("mie_gruneisen", "s2"): r"\f$s_{2,k}\f$",
    ("mie_gruneisen", "s3"): r"\f$s_{3,k}\f$",
    ("jwl", "a"): r"\f$A_k\f$",
    ("jwl", "b"): r"\f$B_k\f$",
    ("jwl", "r1"): r"\f$R_{1,k}\f$",
    ("jwl", "r2"): r"\f$R_{2,k}\f$",
    ("jwl", "omega"): r"\f$\omega_k\f$",
    ("jwl", "rho0"): r"\f$\rho_{0,k}\f$",
    ("jwl", "t0"): r"\f$T_{0,k}\f$",
    ("vinet", "k0"): r"\f$K_{0,k}\f$",
    ("vinet", "k0p"): r"\f$K'_{0,k}\f$",
    ("vinet", "rho0"): r"\f$\rho_{0,k}\f$",
    ("vinet", "gruneisen"): r"\f$\Gamma_{G,k}\f$",
    ("vinet", "gruneisen_a"): r"\f$a_k\f$",
    ("vinet", "t0"): r"\f$T_{0,k}\f$",
}


def test_doxygen_math_symbols_match_expected():
    actual = {(family.suffix, suffix): math for family in EOS_FAMILIES for suffix, math in family.required + family.optional}
    assert actual == _EXPECTED_MATH


def _case_dispatched_fields():
    """eos_coeffs fields written by more than one family — the same rule generate_eos_fpp's
    _eos_case_fields uses, reimplemented here so this test does not depend on that function."""
    writers = {}
    for family in EOS_FAMILIES:
        for field_name in family.eos_coeffs:
            writers[field_name] = writers.get(field_name, 0) + 1
    return {name for name, count in writers.items() if count > 1}


def test_every_state_dependent_family_assigns_every_case_dispatched_field():
    case_fields = _case_dispatched_fields()
    assert case_fields  # sanity: the set this test exists to check is non-empty
    for family in _STATE_DEPENDENT:
        missing = case_fields - set(family.eos_coeffs)
        assert not missing, f"family {family.suffix} does not assign case-dispatched field(s) {sorted(missing)}"


def test_eos_coeff_defaults_keys_are_case_dispatched_fields():
    case_fields = _case_dispatched_fields()
    unknown = set(EOS_COEFF_DEFAULTS) - case_fields
    assert not unknown, f"EOS_COEFF_DEFAULTS key(s) {sorted(unknown)} are not case-dispatched fields"
