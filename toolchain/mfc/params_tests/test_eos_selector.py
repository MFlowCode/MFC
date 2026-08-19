"""
Tests for the per-fluid equation-of-state selector, fluid_pp(i)%eos.

Covers the enum itself (Fortran/Python agreement), the readable-name to integer
resolution done by Case, and the check_eos constraints in case_validator.
"""

import unittest

from ..case import Case
from ..case_validator import CaseConstraintError, CaseValidator
from ..common import MFCException
from ..params.definitions import _EOS_NAMES
from ..params.namelist_parser import get_fortran_constants

# A minimal known-good pre_process case. Previously imported from
# params_tests.negative_tests, which was removed in #1717; inlined here because
# num_fluids = 1 and the absence of any chemistry parameter are what give the
# unused-slot and ideal_gas_mixture tests below their meaning.
BASE_CASE = {
    "m": 50,
    "n": 0,
    "p": 0,
    "model_eqns": 2,
    "num_fluids": 1,
    "num_patches": 1,
    "t_step_start": 0,
    "t_step_stop": 100,
    "t_step_save": 10,
    "dt": 1e-6,
    "weno_order": 5,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": 1.0,
    "patch_icpp(1)%alpha_rho(1)": 1.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 0.4,
    "fluid_pp(1)%pi_inf": 0.0,
}


def _eos_errors(overrides):
    """Validate BASE_CASE plus overrides, returning only the eos-related messages."""
    params = dict(BASE_CASE)
    params.update(overrides)
    try:
        CaseValidator(Case(params).params).validate("pre_process")
    except CaseConstraintError as exc:
        return [line for line in str(exc).splitlines() if "%eos" in line]
    return []


class TestEosEnum(unittest.TestCase):
    """The enum is hand-written in m_constants.fpp and restated in definitions.py."""

    def test_fortran_and_python_enums_agree(self):
        """_EOS_NAMES must match the eos_* parameters in m_constants.fpp.

        generate_constants_fpp skips compound registry keys, so these constants are
        hand-written on the Fortran side and nothing else forces the two to agree.
        """
        fortran = {name[len("eos_") :]: value for name, value in get_fortran_constants().items() if name.startswith("eos_")}
        self.assertEqual(fortran, _EOS_NAMES)

    def test_only_implemented_backends_are_exposed(self):
        """Reserved values must not appear until they have a backend and a check_eos branch."""
        self.assertEqual(set(_EOS_NAMES), {"stiffened_gas", "ideal_gas_mixture"})


class TestEosNameResolution(unittest.TestCase):
    """Case converts the readable name in a case file to the integer the namelist carries."""

    def test_name_resolves_to_integer(self):
        case = Case({"fluid_pp(1)%eos": "stiffened_gas"})
        self.assertEqual(case.params["fluid_pp(1)%eos"], _EOS_NAMES["stiffened_gas"])

    def test_integer_passes_through(self):
        case = Case({"fluid_pp(1)%eos": _EOS_NAMES["ideal_gas_mixture"]})
        self.assertEqual(case.params["fluid_pp(1)%eos"], _EOS_NAMES["ideal_gas_mixture"])

    def test_unknown_name_rejected(self):
        with self.assertRaises(MFCException) as ctx:
            Case({"fluid_pp(1)%eos": "jwl"})
        self.assertIn("stiffened_gas", str(ctx.exception))

    def test_resolution_applies_to_every_fluid_slot(self):
        """CONSTRAINTS is registered per slot, so slot 10 must resolve like slot 1."""
        case = Case({"fluid_pp(10)%eos": "stiffened_gas"})
        self.assertEqual(case.params["fluid_pp(10)%eos"], _EOS_NAMES["stiffened_gas"])


class TestCheckEos(unittest.TestCase):
    """check_eos constraints, on a non-chemistry build (BASE_CASE sets no chemistry)."""

    def test_base_case_has_no_eos_errors(self):
        self.assertEqual(_eos_errors({}), [])

    def test_stiffened_gas_accepted(self):
        self.assertEqual(_eos_errors({"fluid_pp(1)%eos": "stiffened_gas"}), [])

    def test_ideal_gas_mixture_requires_chemistry(self):
        errors = _eos_errors({"fluid_pp(1)%eos": "ideal_gas_mixture"})
        self.assertTrue(errors)
        self.assertIn("requires a chemistry build", " ".join(errors))

    def test_value_outside_enum_rejected(self):
        """The choices constraint covers integers the enum does not define."""
        errors = _eos_errors({"fluid_pp(1)%eos": 99})
        self.assertTrue(errors)

    def test_unused_slot_is_validated(self):
        """Slots above num_fluids are default-assigned and broadcast, so they are checked too.

        BASE_CASE sets num_fluids = 1; fluid_pp(3) is an unused slot whose value still
        reaches every rank through the fluid_pp member-loop broadcast.
        """
        errors = _eos_errors({"fluid_pp(3)%eos": "ideal_gas_mixture"})
        self.assertTrue(errors)
        self.assertIn("fluid_pp(3)%eos", " ".join(errors))


if __name__ == "__main__":
    unittest.main()
