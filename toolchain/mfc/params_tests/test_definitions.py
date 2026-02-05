"""
Unit tests for params/definitions.py module.

Tests parameter definitions, constraints, and dependencies.
"""

import unittest
from ..params import REGISTRY
from ..params.schema import ParamType
from ..params.definitions import (
    CONSTRAINTS,
    DEPENDENCIES,
    CASE_OPT_PARAMS,
    _validate_constraint,
    _validate_dependency,
)


class TestParameterDefinitions(unittest.TestCase):
    """Tests for parameter definitions."""

    def test_all_params_have_names(self):
        """Every parameter should have a non-empty name."""
        for name, param in REGISTRY.all_params.items():
            self.assertEqual(name, param.name)
            self.assertTrue(len(name) > 0)

    def test_all_params_have_valid_type(self):
        """Every parameter should have a valid ParamType."""
        for name, param in REGISTRY.all_params.items():
            self.assertIsInstance(
                param.param_type, ParamType,
                f"Parameter '{name}' has invalid type"
            )

    def test_core_params_exist(self):
        """Core parameters m, n, p, model_eqns should exist."""
        core_params = ["m", "n", "p", "model_eqns", "num_fluids"]
        for param_name in core_params:
            self.assertIn(
                param_name, REGISTRY.all_params,
                f"Core parameter '{param_name}' not found"
            )

    def test_domain_params_exist(self):
        """Domain parameters should exist for all directions."""
        for d in ["x", "y", "z"]:
            self.assertIn(f"{d}_domain%beg", REGISTRY.all_params)
            self.assertIn(f"{d}_domain%end", REGISTRY.all_params)
            self.assertIn(f"bc_{d}%beg", REGISTRY.all_params)
            self.assertIn(f"bc_{d}%end", REGISTRY.all_params)


class TestConstraintValidation(unittest.TestCase):
    """Tests for constraint schema validation."""

    def test_valid_choices_constraint(self):
        """Valid choices constraint should pass."""
        _validate_constraint("test", {"choices": [1, 2, 3]})

    def test_valid_range_constraint(self):
        """Valid range constraint should pass."""
        _validate_constraint("test", {"min": 0, "max": 100})

    def test_invalid_key_raises(self):
        """Invalid constraint key should raise ValueError."""
        with self.assertRaises(ValueError) as ctx:
            _validate_constraint("test", {"chioces": [1, 2]})  # Typo for choices

        self.assertIn("chioces", str(ctx.exception))
        # Verify "did you mean?" suggestion is provided
        self.assertIn("choices", str(ctx.exception))

    def test_choices_must_be_list(self):
        """choices value must be a list."""
        with self.assertRaises(ValueError):
            _validate_constraint("test", {"choices": "not a list"})

    def test_min_must_be_number(self):
        """min value must be a number."""
        with self.assertRaises(ValueError):
            _validate_constraint("test", {"min": "not a number"})

    def test_all_defined_constraints_are_valid(self):
        """All constraints in CONSTRAINTS dict should be valid."""
        for param_name, constraint in CONSTRAINTS.items():
            try:
                _validate_constraint(param_name, constraint)
            except ValueError as e:
                self.fail(f"Invalid constraint for '{param_name}': {e}")


class TestDependencyValidation(unittest.TestCase):
    """Tests for dependency schema validation."""

    def test_valid_when_true_dependency(self):
        """Valid when_true dependency should pass."""
        _validate_dependency("test", {
            "when_true": {
                "requires": ["other_param"],
                "recommends": ["another_param"],
            }
        })

    def test_invalid_top_level_key_raises(self):
        """Invalid top-level dependency key should raise."""
        with self.assertRaises(ValueError) as ctx:
            _validate_dependency("test", {"when_tru": {"requires": []}})  # Typo

        self.assertIn("when_tru", str(ctx.exception))

    def test_invalid_condition_key_raises(self):
        """Invalid condition key should raise."""
        with self.assertRaises(ValueError) as ctx:
            _validate_dependency("test", {
                "when_true": {"reqires": ["foo"]}  # Typo for requires
            })

        self.assertIn("reqires", str(ctx.exception))
        # Verify "did you mean?" suggestion is provided
        self.assertIn("requires", str(ctx.exception))

    def test_requires_must_be_list(self):
        """requires value must be a list."""
        with self.assertRaises(ValueError):
            _validate_dependency("test", {
                "when_true": {"requires": "not a list"}
            })

    def test_all_defined_dependencies_are_valid(self):
        """All dependencies in DEPENDENCIES dict should be valid."""
        for param_name, dependency in DEPENDENCIES.items():
            try:
                _validate_dependency(param_name, dependency)
            except ValueError as e:
                self.fail(f"Invalid dependency for '{param_name}': {e}")


class TestCaseOptimizationParams(unittest.TestCase):
    """Tests for case optimization parameter set."""

    def test_case_opt_params_exist_in_registry(self):
        """All CASE_OPT_PARAMS should exist in registry."""
        for param_name in CASE_OPT_PARAMS:
            self.assertIn(
                param_name, REGISTRY.all_params,
                f"Case opt param '{param_name}' not in registry"
            )

    def test_case_opt_params_have_flag_set(self):
        """Params in CASE_OPT_PARAMS should have case_optimization=True."""
        for param_name in CASE_OPT_PARAMS:
            param = REGISTRY.all_params[param_name]
            self.assertTrue(
                param.case_optimization,
                f"Parameter '{param_name}' should have case_optimization=True"
            )


class TestParameterCounts(unittest.TestCase):
    """Tests for expected parameter counts."""

    def test_total_param_count(self):
        """Total parameter count should be around 3400."""
        count = len(REGISTRY.all_params)
        self.assertGreater(count, 3000, "Too few parameters")
        self.assertLess(count, 4000, "Too many parameters")

    def test_log_params_count(self):
        """Should have many LOG type parameters."""
        log_count = sum(
            1 for p in REGISTRY.all_params.values()
            if p.param_type == ParamType.LOG
        )
        self.assertGreater(log_count, 300, "Too few LOG parameters")

    def test_tagged_params_exist(self):
        """Should have params with feature tags."""
        mhd_params = REGISTRY.get_params_by_tag("mhd")
        self.assertGreater(len(mhd_params), 5, "Too few MHD parameters")

        bubbles_params = REGISTRY.get_params_by_tag("bubbles")
        self.assertGreater(len(bubbles_params), 10, "Too few bubbles parameters")


if __name__ == "__main__":
    unittest.main()
