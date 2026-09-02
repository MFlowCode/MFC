"""Guards the per-fluid equation-of-state selector.

The eos_* constants are hand-written in m_constants.fpp because generate_constants_fpp skips
compound registry keys (fluid_pp(1)%eos is not a valid Fortran identifier), so nothing but this
test keeps the Fortran and Python enumerations in step.
"""

from ..case_validator import CaseValidator
from ..params.definitions import CONSTRAINTS
from ..params.namelist_parser import get_fortran_constants

EOS = CONSTRAINTS["fluid_pp(1)%eos"]["names"]


def _errors(**overrides):
    params = {"num_fluids": 1, "fluid_pp(1)%gamma": 2.5, "fluid_pp(1)%pi_inf": 0.0}
    params.update(overrides)
    validator = CaseValidator(params)
    validator.check_eos_selector()
    return validator.errors


def test_fortran_and_python_enums_agree():
    fortran = get_fortran_constants()
    for name, value in EOS.items():
        assert fortran[f"eos_{name}"] == value, f"eos_{name} disagrees between m_constants.fpp and definitions.py"


def test_stiffened_gas_accepts_stiffness():
    assert _errors(**{"fluid_pp(1)%eos": EOS["stiffened_gas"], "fluid_pp(1)%pi_inf": 1.0e5}) == []


def test_ideal_gas_accepts_zero_stiffness():
    assert _errors(**{"fluid_pp(1)%eos": EOS["ideal_gas"], "fluid_pp(1)%pi_inf": 0.0}) == []


def test_ideal_gas_rejects_nonzero_stiffness():
    errors = _errors(**{"fluid_pp(1)%eos": EOS["ideal_gas"], "fluid_pp(1)%pi_inf": 1.0e5})
    assert len(errors) == 1 and "pi_inf" in errors[0]


def test_value_outside_the_enumeration_is_rejected():
    # choices are enforced by validate_constraints, a different layer, so the membership
    # check here is not redundant.
    assert len(_errors(**{"fluid_pp(1)%eos": 99})) == 1


def test_unset_selector_is_skipped():
    assert _errors() == []
