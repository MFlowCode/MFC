"""
Integration tests for params module with case_dicts.

Tests that the parameter registry integrates correctly with case_dicts.py
and provides correct JSON schema generation.
"""
# pylint: disable=import-outside-toplevel

import unittest
from ..params import REGISTRY
from ..params.schema import ParamType


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

    def test_core_params_in_schema(self):
        """Core params should be in JSON schema."""
        schema = REGISTRY.get_json_schema()
        props = schema["properties"]

        self.assertIn("m", props)
        self.assertIn("n", props)
        self.assertIn("p", props)
        self.assertIn("model_eqns", props)


class TestRegistryTagQueries(unittest.TestCase):
    """Tests for tag-based parameter queries."""

    def test_get_params_by_tag(self):
        """Should get params by feature tag."""
        mhd_params = REGISTRY.get_params_by_tag("mhd")

        self.assertIsInstance(mhd_params, dict)
        self.assertGreater(len(mhd_params), 5)
        self.assertIn("mhd", mhd_params)

    def test_get_all_tags(self):
        """Should get all registered tags."""
        tags = REGISTRY.get_all_tags()

        self.assertIsInstance(tags, set)
        self.assertIn("mhd", tags)
        self.assertIn("bubbles", tags)
        self.assertIn("weno", tags)

    def test_params_have_correct_tags(self):
        """Parameters should have their expected tags."""
        mhd_param = REGISTRY.all_params.get("mhd")
        self.assertIsNotNone(mhd_param)
        self.assertIn("mhd", mhd_param.tags)

        bubbles_param = REGISTRY.all_params.get("bubbles_euler")
        self.assertIsNotNone(bubbles_param)
        self.assertIn("bubbles", bubbles_param.tags)


class TestCaseDictsIntegration(unittest.TestCase):
    """Tests for integration with actual case_dicts module."""

    def test_case_dicts_loads_from_registry(self):
        """case_dicts module should load successfully from registry."""
        from ..run import case_dicts

        # ALL should be populated
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

    def test_get_input_dict_keys(self):
        """get_input_dict_keys should return target-specific params."""
        from ..run import case_dicts

        # Each target gets a filtered subset of params based on Fortran namelists
        pre_keys = case_dicts.get_input_dict_keys("pre_process")
        sim_keys = case_dicts.get_input_dict_keys("simulation")
        post_keys = case_dicts.get_input_dict_keys("post_process")

        # pre_process has most params (includes patch_icpp, patch_bc)
        self.assertGreater(len(pre_keys), 2500)
        # simulation and post_process have fewer (no patch_icpp, etc.)
        self.assertGreater(len(sim_keys), 500)
        self.assertGreater(len(post_keys), 400)

        # Verify target-specific filtering based on Fortran namelists
        self.assertIn("num_patches", pre_keys)
        self.assertNotIn("num_patches", sim_keys)
        self.assertNotIn("num_patches", post_keys)

        self.assertNotIn("run_time_info", pre_keys)
        self.assertIn("run_time_info", sim_keys)
        self.assertNotIn("run_time_info", post_keys)

        # Verify indexed params are filtered correctly
        patch_icpp_pre = [k for k in pre_keys if k.startswith("patch_icpp")]
        patch_icpp_sim = [k for k in sim_keys if k.startswith("patch_icpp")]
        self.assertGreater(len(patch_icpp_pre), 1000)  # Many patch_icpp params
        self.assertEqual(len(patch_icpp_sim), 0)  # None in simulation

        # Verify shared params are in all targets
        self.assertIn("m", pre_keys)
        self.assertIn("m", sim_keys)
        self.assertIn("m", post_keys)


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
