"""
Unit tests for params/validate.py module.

Tests constraint validation, dependency checking, and error formatting.
"""

import unittest
from ..params.validate import (
    validate_constraints,
    check_dependencies,
    validate_case,
    format_validation_results,
)


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
