"""
Unit tests for params/registry.py module.

Tests registry functionality, freezing, and tag queries.
"""
# pylint: disable=import-outside-toplevel

import unittest
from ..params.registry import ParamRegistry, RegistryFrozenError
from ..params.schema import ParamDef, ParamType


class TestParamRegistry(unittest.TestCase):
    """Tests for ParamRegistry class."""

    def test_register_new_param(self):
        """Registering a new param should add it to the registry."""
        reg = ParamRegistry()
        param = ParamDef(name="test", param_type=ParamType.INT)
        reg.register(param)

        self.assertIn("test", reg.all_params)
        self.assertEqual(reg.all_params["test"].param_type, ParamType.INT)

    def test_register_merge_tags(self):
        """Registering same param twice should merge tags."""
        reg = ParamRegistry()
        reg.register(ParamDef(name="test", param_type=ParamType.INT, tags={"mhd"}))
        reg.register(ParamDef(name="test", param_type=ParamType.INT, tags={"physics"}))

        self.assertEqual(
            reg.all_params["test"].tags,
            {"mhd", "physics"}
        )

    def test_register_type_mismatch_raises(self):
        """Registering same param with different type should raise."""
        reg = ParamRegistry()
        reg.register(ParamDef(name="test", param_type=ParamType.INT))

        with self.assertRaises(ValueError) as ctx:
            reg.register(ParamDef(name="test", param_type=ParamType.REAL))

        self.assertIn("Type mismatch", str(ctx.exception))

    def test_get_params_by_tag(self):
        """get_params_by_tag should return params with that tag."""
        reg = ParamRegistry()
        reg.register(ParamDef(name="mhd_param", param_type=ParamType.INT, tags={"mhd"}))
        reg.register(ParamDef(name="other_param", param_type=ParamType.INT, tags={"other"}))

        mhd_params = reg.get_params_by_tag("mhd")

        self.assertIn("mhd_param", mhd_params)
        self.assertNotIn("other_param", mhd_params)

    def test_get_all_tags(self):
        """get_all_tags should return all registered tags."""
        reg = ParamRegistry()
        reg.register(ParamDef(name="p1", param_type=ParamType.INT, tags={"mhd", "physics"}))
        reg.register(ParamDef(name="p2", param_type=ParamType.INT, tags={"bubbles"}))

        tags = reg.get_all_tags()

        self.assertEqual(tags, {"mhd", "physics", "bubbles"})


class TestRegistryFreezing(unittest.TestCase):
    """Tests for registry freeze functionality."""

    def test_freeze_sets_flag(self):
        """freeze() should set is_frozen to True."""
        reg = ParamRegistry()
        self.assertFalse(reg.is_frozen)

        reg.freeze()

        self.assertTrue(reg.is_frozen)

    def test_freeze_is_idempotent(self):
        """Calling freeze() multiple times should be safe."""
        reg = ParamRegistry()
        reg.freeze()
        reg.freeze()  # Should not raise

        self.assertTrue(reg.is_frozen)

    def test_register_after_freeze_raises(self):
        """Registering after freeze should raise RegistryFrozenError."""
        reg = ParamRegistry()
        reg.freeze()

        with self.assertRaises(RegistryFrozenError):
            reg.register(ParamDef(name="test", param_type=ParamType.INT))

    def test_all_params_readonly_after_freeze(self):
        """all_params should be read-only after freeze."""
        reg = ParamRegistry()
        reg.register(ParamDef(name="test", param_type=ParamType.INT))
        reg.freeze()

        params = reg.all_params

        with self.assertRaises(TypeError):
            params["hacked"] = "value"


class TestGlobalRegistry(unittest.TestCase):
    """Tests for the global REGISTRY instance."""

    def test_global_registry_is_frozen(self):
        """Global REGISTRY should be frozen after import."""
        from ..params import REGISTRY
        self.assertTrue(REGISTRY.is_frozen)

    def test_global_registry_has_params(self):
        """Global REGISTRY should have parameters loaded."""
        from ..params import REGISTRY
        self.assertGreater(len(REGISTRY.all_params), 3000)

    def test_global_registry_cannot_be_modified(self):
        """Global REGISTRY should reject new registrations."""
        from ..params import REGISTRY

        with self.assertRaises(RegistryFrozenError):
            REGISTRY.register(
                ParamDef(name="injected", param_type=ParamType.INT)
            )


if __name__ == "__main__":
    unittest.main()
