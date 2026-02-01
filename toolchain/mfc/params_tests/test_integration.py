"""
Integration tests for params module with case_dicts.

Tests that the parameter registry integrates correctly with case_dicts.py
and provides correct JSON schema generation.
"""
# pylint: disable=import-outside-toplevel

import unittest
from ..params import REGISTRY
from ..params.schema import Stage, ParamType


class TestParamTypeJsonSchema(unittest.TestCase):
    """Tests for ParamType.json_schema property."""

    def test_all_types_have_json_schema(self):
        """Every ParamType should have a json_schema property."""
        for param_type in ParamType:
            schema = param_type.json_schema
            self.assertIsInstance(schema, dict)
            # Schema must have either "type" or "enum" key
            self.assertTrue(
                "type" in schema or "enum" in schema,
                f"{param_type.name} schema has neither 'type' nor 'enum'"
            )

    def test_int_schema(self):
        """INT should map to integer JSON schema."""
        schema = ParamType.INT.json_schema
        self.assertEqual(schema, {"type": "integer"})

    def test_real_schema(self):
        """REAL should map to number JSON schema."""
        schema = ParamType.REAL.json_schema
        self.assertEqual(schema, {"type": "number"})

    def test_log_schema(self):
        """LOG should map to enum with T/F."""
        schema = ParamType.LOG.json_schema
        self.assertEqual(schema, {"enum": ["T", "F"]})

    def test_str_schema(self):
        """STR should map to string JSON schema."""
        schema = ParamType.STR.json_schema
        self.assertEqual(schema, {"type": "string"})

    def test_analytic_int_schema(self):
        """ANALYTIC_INT should accept integer or string."""
        schema = ParamType.ANALYTIC_INT.json_schema
        self.assertEqual(schema, {"type": ["integer", "string"]})

    def test_analytic_real_schema(self):
        """ANALYTIC_REAL should accept number or string."""
        schema = ParamType.ANALYTIC_REAL.json_schema
        self.assertEqual(schema, {"type": ["number", "string"]})


class TestRegistryJsonSchema(unittest.TestCase):
    """Tests for registry JSON schema generation."""

    def test_get_json_schema_returns_valid_schema(self):
        """get_json_schema should return valid JSON schema structure."""
        schema = REGISTRY.get_json_schema()

        self.assertIn("type", schema)
        self.assertEqual(schema["type"], "object")
        self.assertIn("properties", schema)
        self.assertIn("additionalProperties", schema)
        self.assertFalse(schema["additionalProperties"])

    def test_get_json_schema_has_all_params(self):
        """Schema properties should include all registry params."""
        schema = REGISTRY.get_json_schema()
        properties = schema["properties"]

        self.assertEqual(len(properties), len(REGISTRY.all_params))

    def test_get_json_schema_for_stage(self):
        """get_json_schema(stage) should only include stage params."""
        full_schema = REGISTRY.get_json_schema()
        sim_schema = REGISTRY.get_json_schema(stage=Stage.SIMULATION)

        # Stage-specific should be subset of full
        self.assertLessEqual(
            len(sim_schema["properties"]),
            len(full_schema["properties"])
        )

    def test_core_params_in_schema(self):
        """Core params should be in JSON schema."""
        schema = REGISTRY.get_json_schema()
        props = schema["properties"]

        self.assertIn("m", props)
        self.assertIn("n", props)
        self.assertIn("p", props)
        self.assertIn("model_eqns", props)


class TestRegistryStageQueries(unittest.TestCase):
    """Tests for stage-based parameter queries."""

    def test_common_stage_params(self):
        """Should get COMMON stage params."""
        params = REGISTRY.get_params_by_stage(Stage.COMMON)

        self.assertIsInstance(params, dict)
        self.assertGreater(len(params), 0)
        self.assertIn("m", params)

    def test_simulation_stage_includes_common(self):
        """SIMULATION stage should include COMMON params."""
        sim_params = REGISTRY.get_params_by_stage(Stage.SIMULATION)

        # m is a COMMON param, should be included
        self.assertIn("m", sim_params)

    def test_pre_process_stage(self):
        """Should get PRE_PROCESS stage params."""
        params = REGISTRY.get_params_by_stage(Stage.PRE_PROCESS)

        self.assertIsInstance(params, dict)
        self.assertGreater(len(params), 100)

    def test_post_process_stage(self):
        """Should get POST_PROCESS stage params."""
        params = REGISTRY.get_params_by_stage(Stage.POST_PROCESS)

        self.assertIsInstance(params, dict)
        self.assertGreater(len(params), 100)


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
        self.assertEqual(len(case_dicts.ALL), len(REGISTRY.all_params))

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
        self.assertEqual(len(schema["properties"]), len(REGISTRY.all_params))

    def test_json_schema_matches_registry(self):
        """case_dicts.SCHEMA should match REGISTRY.get_json_schema()."""
        from ..run import case_dicts

        registry_schema = REGISTRY.get_json_schema()
        case_dicts_schema = case_dicts.SCHEMA

        self.assertEqual(registry_schema, case_dicts_schema)

    def test_get_validator_works(self):
        """get_validator should return a callable validator."""
        from ..run import case_dicts

        validator = case_dicts.get_validator()
        self.assertTrue(callable(validator))


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
