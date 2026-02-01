"""
Unit tests for params/validate.py module.

Tests constraint validation, dependency checking, and error formatting.
"""

import unittest
from ..params.validate import (
    validate_constraints,
    check_dependencies,
    check_unknown_params,
    validate_case,
    format_validation_results,
)
from ..params.suggest import RAPIDFUZZ_AVAILABLE


class TestValidateConstraints(unittest.TestCase):
    """Tests for validate_constraints function."""

    def test_valid_params_no_errors(self):
        """Valid parameters should return empty error list."""
        params = {
            "m": 100,
            "n": 50,
            "dt": 1e-5,
        }
        errors = validate_constraints(params)
        self.assertEqual(errors, [])

    def test_unknown_params_ignored(self):
        """Unknown parameters should be silently ignored."""
        params = {
            "unknown_param_xyz": 123,
            "another_unknown": "value",
        }
        errors = validate_constraints(params)
        self.assertEqual(errors, [])

    def test_analytic_string_skipped(self):
        """Analytic expressions (strings) should skip numeric validation."""
        params = {
            "dt": "1e-5 * m / 100",  # String expression for numeric param
        }
        errors = validate_constraints(params)
        self.assertEqual(errors, [])

    def test_constraint_validation_called(self):
        """Constraint validation should be invoked for registered params."""
        # This tests the integration with ParamDef.validate_value
        params = {
            "m": 100,  # Valid integer
        }
        errors = validate_constraints(params)
        self.assertEqual(errors, [])


class TestCheckUnknownParams(unittest.TestCase):
    """Tests for check_unknown_params with 'did you mean?' suggestions."""

    def test_known_params_no_errors(self):
        """Known parameters should not generate errors."""
        params = {
            "m": 100,
            "n": 50,
            "dt": 1e-5,
        }
        errors = check_unknown_params(params)
        self.assertEqual(errors, [])

    def test_unknown_param_returns_error(self):
        """Unknown parameter should return an error."""
        params = {
            "totally_unknown_xyz_123": 42,
        }
        errors = check_unknown_params(params)
        self.assertEqual(len(errors), 1)
        self.assertIn("Unknown parameter", errors[0])
        self.assertIn("totally_unknown_xyz_123", errors[0])

    @unittest.skipUnless(RAPIDFUZZ_AVAILABLE, "rapidfuzz not installed")
    def test_similar_param_suggests_correction(self):
        """Typo near valid param should suggest 'did you mean?'."""
        # "model_eqn" is close to "model_eqns"
        params = {
            "model_eqn": 2,  # Missing 's'
        }
        errors = check_unknown_params(params)
        self.assertEqual(len(errors), 1)
        self.assertIn("model_eqn", errors[0])
        # Should suggest the correct parameter
        self.assertIn("model_eqns", errors[0])
        self.assertIn("Did you mean", errors[0])

    @unittest.skipUnless(RAPIDFUZZ_AVAILABLE, "rapidfuzz not installed")
    def test_weno_typo_suggests_correction(self):
        """Typo in weno_order should suggest correction."""
        params = {
            "weno_ordr": 5,  # Typo for weno_order
        }
        errors = check_unknown_params(params)
        self.assertEqual(len(errors), 1)
        self.assertIn("weno_order", errors[0])


class TestCheckDependencies(unittest.TestCase):
    """Tests for check_dependencies function."""

    def test_no_dependencies_empty_result(self):
        """Params without dependencies return empty errors/warnings."""
        params = {
            "m": 100,
            "n": 50,
        }
        errors, warnings = check_dependencies(params)
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])

    def test_unknown_params_ignored(self):
        """Unknown parameters should be silently ignored."""
        params = {
            "unknown_param_xyz": "T",
        }
        errors, warnings = check_dependencies(params)
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])


class TestValidateCase(unittest.TestCase):
    """Tests for validate_case function."""

    def test_valid_case_no_errors(self):
        """Valid case should return empty errors."""
        params = {
            "m": 100,
            "n": 50,
            "p": 0,
            "dt": 1e-5,
        }
        errors, _ = validate_case(params)
        self.assertEqual(errors, [])

    def test_warn_false_skips_warnings(self):
        """warn=False should skip dependency warnings."""
        params = {"m": 100}
        _, warnings = validate_case(params, warn=False)
        self.assertEqual(warnings, [])

    def test_returns_tuple(self):
        """Should return tuple of (errors, warnings)."""
        params = {"m": 100}
        result = validate_case(params)
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        self.assertIsInstance(result[0], list)
        self.assertIsInstance(result[1], list)

    def test_check_unknown_true_catches_unknown(self):
        """check_unknown=True should catch unknown parameters."""
        params = {
            "m": 100,
            "unknwn_prm": 42,  # Unknown param
        }
        errors, _ = validate_case(params, check_unknown=True)
        self.assertGreater(len(errors), 0)
        self.assertIn("unknwn_prm", str(errors))

    def test_check_unknown_false_ignores_unknown(self):
        """check_unknown=False should ignore unknown parameters."""
        params = {
            "m": 100,
            "unknwn_prm": 42,  # Unknown param
        }
        errors, _ = validate_case(params, check_unknown=False)
        # No errors from unknown params (may have other errors)
        unknown_errors = [e for e in errors if "unknwn_prm" in e]
        self.assertEqual(unknown_errors, [])


class TestFormatValidationResults(unittest.TestCase):
    """Tests for format_validation_results function."""

    def test_empty_results(self):
        """Empty errors and warnings should return empty string."""
        result = format_validation_results([], [])
        self.assertEqual(result, "")

    def test_errors_only(self):
        """Format should include errors with red markers."""
        errors = ["Error 1", "Error 2"]
        result = format_validation_results(errors, [])
        self.assertIn("Validation Errors", result)
        self.assertIn("Error 1", result)
        self.assertIn("Error 2", result)

    def test_warnings_only(self):
        """Format should include warnings with yellow markers."""
        warnings = ["Warning 1"]
        result = format_validation_results([], warnings)
        self.assertIn("Warnings", result)
        self.assertIn("Warning 1", result)

    def test_both_errors_and_warnings(self):
        """Format should include both errors and warnings."""
        errors = ["Error 1"]
        warnings = ["Warning 1"]
        result = format_validation_results(errors, warnings)
        self.assertIn("Validation Errors", result)
        self.assertIn("Error 1", result)
        self.assertIn("Warnings", result)
        self.assertIn("Warning 1", result)

    def test_rich_markup_included(self):
        """Format should include Rich markup for colors."""
        errors = ["Error 1"]
        result = format_validation_results(errors, [])
        self.assertIn("[red]", result)


class TestSchemaValidation(unittest.TestCase):
    """Tests for ParamDef.validate_value through validate_constraints."""

    def test_numeric_type_check_for_min_max(self):
        """String values should not crash on numeric comparisons."""
        # This tests the fix for the type comparison bug
        params = {
            "m": "100 + 50",  # String that looks numeric but isn't
        }
        # Should not raise TypeError
        errors = validate_constraints(params)
        self.assertIsInstance(errors, list)


if __name__ == "__main__":
    unittest.main()
