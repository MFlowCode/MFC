"""
Integration tests for params module with case_dicts.

Tests that the parameter registry integrates correctly with case_dicts.py.
"""
# pylint: disable=import-outside-toplevel

import unittest
from ..params import REGISTRY
from ..params.schema import Stage, ParamType
from ..params.generators.case_dicts_gen import (
    get_stage_dict,
    verify_type_systems,
)


class TestTypeSystems(unittest.TestCase):
    """Tests for type system synchronization."""

    def test_type_systems_are_synchronized(self):
        """verify_type_systems should return no errors."""
        errors = verify_type_systems()
        self.assertEqual(errors, [], f"Type system errors: {errors}")

    def test_all_registry_types_have_mapping(self):
        """Every ParamType in registry should have a case_dicts mapping."""
        from ..params.generators.case_dicts_gen import _get_type_map

        type_map = _get_type_map()

        for param_type in ParamType:
            self.assertIn(
                param_type, type_map,
                f"ParamType.{param_type.name} has no case_dicts mapping"
            )


class TestStageGeneration(unittest.TestCase):
    """Tests for stage dictionary generation."""

    def test_generate_common_stage(self):
        """Should generate COMMON stage dict."""
        result = get_stage_dict(Stage.COMMON, include_common=False)

        self.assertIsInstance(result, dict)
        self.assertGreater(len(result), 0)

    def test_generate_pre_process_stage(self):
        """Should generate PRE_PROCESS stage dict with COMMON params."""
        result = get_stage_dict(Stage.PRE_PROCESS, include_common=True)

        self.assertIsInstance(result, dict)
        self.assertGreater(len(result), 100)

    def test_generate_simulation_stage(self):
        """Should generate SIMULATION stage dict."""
        result = get_stage_dict(Stage.SIMULATION, include_common=True)

        self.assertIsInstance(result, dict)
        self.assertGreater(len(result), 100)

    def test_generate_post_process_stage(self):
        """Should generate POST_PROCESS stage dict."""
        result = get_stage_dict(Stage.POST_PROCESS, include_common=True)

        self.assertIsInstance(result, dict)
        self.assertGreater(len(result), 100)

    def test_include_common_adds_common_params(self):
        """include_common=True should add COMMON params to other stages."""
        without_common = get_stage_dict(Stage.SIMULATION, include_common=False)
        with_common = get_stage_dict(Stage.SIMULATION, include_common=True)

        self.assertGreater(len(with_common), len(without_common))

    def test_core_params_in_common(self):
        """Core params like m, n, p should be in COMMON stage."""
        result = get_stage_dict(Stage.COMMON, include_common=False)

        self.assertIn("m", result)
        self.assertIn("n", result)
        self.assertIn("p", result)


class TestCaseDictsIntegration(unittest.TestCase):
    """Tests for integration with actual case_dicts module."""

    def test_case_dicts_loads_from_registry(self):
        """case_dicts module should load successfully from registry."""
        from ..run import case_dicts

        # These should all be populated
        self.assertIsNotNone(case_dicts.COMMON)
        self.assertIsNotNone(case_dicts.PRE_PROCESS)
        self.assertIsNotNone(case_dicts.SIMULATION)
        self.assertIsNotNone(case_dicts.POST_PROCESS)
        self.assertIsNotNone(case_dicts.ALL)

    def test_case_dicts_all_contains_all_params(self):
        """case_dicts.ALL should contain all registry params."""
        from ..run import case_dicts

        # ALL should have approximately the same params as registry
        self.assertGreater(len(case_dicts.ALL), 3000)

    def test_case_optimization_params_from_registry(self):
        """CASE_OPTIMIZATION should be populated from registry."""
        from ..run import case_dicts

        self.assertIsInstance(case_dicts.CASE_OPTIMIZATION, list)
        self.assertGreater(len(case_dicts.CASE_OPTIMIZATION), 10)

    def test_json_schema_valid(self):
        """Generated JSON schema should be valid."""
        from ..run import case_dicts

        schema = case_dicts.SCHEMA
        self.assertIn("type", schema)
        self.assertEqual(schema["type"], "object")
        self.assertIn("properties", schema)
        self.assertGreater(len(schema["properties"]), 3000)


class TestValidatorIntegration(unittest.TestCase):
    """Tests for integration with case_validator."""

    def test_validator_gets_log_params_from_registry(self):
        """Validator should discover LOG params from registry."""
        from ..case_validator import _get_logical_params_from_registry

        log_params = _get_logical_params_from_registry()

        self.assertIsInstance(log_params, set)
        self.assertGreater(len(log_params), 300)

    def test_validator_log_params_match_registry(self):
        """Validator LOG params should match registry LOG params."""
        from ..case_validator import _get_logical_params_from_registry

        validator_log_params = _get_logical_params_from_registry()

        registry_log_params = {
            name for name, p in REGISTRY.all_params.items()
            if p.param_type == ParamType.LOG
        }

        self.assertEqual(validator_log_params, registry_log_params)


if __name__ == "__main__":
    unittest.main()
