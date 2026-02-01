"""
Unit tests for params/registry.py module.

Tests registry functionality, freezing, and stage queries.
"""

import unittest
from ..params.registry import ParamRegistry, RegistryFrozenError
from ..params.schema import ParamDef, ParamType, Stage


class TestParamRegistry(unittest.TestCase):
    """Tests for ParamRegistry class."""

    def test_register_new_param(self):
        """Registering a new param should add it to the registry."""
        reg = ParamRegistry()
        param = ParamDef(name="test", param_type=ParamType.INT, stages={Stage.COMMON})
        reg.register(param)

        self.assertIn("test", reg.all_params)
        self.assertEqual(reg.all_params["test"].param_type, ParamType.INT)

    def test_register_merge_stages(self):
        """Registering same param twice should merge stages."""
        reg = ParamRegistry()
        reg.register(ParamDef(name="test", param_type=ParamType.INT, stages={Stage.PRE_PROCESS}))
        reg.register(ParamDef(name="test", param_type=ParamType.INT, stages={Stage.SIMULATION}))

        self.assertEqual(
            reg.all_params["test"].stages,
            {Stage.PRE_PROCESS, Stage.SIMULATION}
        )

    def test_register_type_mismatch_raises(self):
        """Registering same param with different type should raise."""
        reg = ParamRegistry()
        reg.register(ParamDef(name="test", param_type=ParamType.INT, stages={Stage.COMMON}))

        with self.assertRaises(ValueError) as ctx:
            reg.register(ParamDef(name="test", param_type=ParamType.REAL, stages={Stage.COMMON}))

        self.assertIn("Type mismatch", str(ctx.exception))

    def test_get_params_by_stage_common(self):
        """get_params_by_stage should return COMMON params for any stage."""
        reg = ParamRegistry()
        reg.register(ParamDef(name="common_param", param_type=ParamType.INT, stages={Stage.COMMON}))
        reg.register(ParamDef(name="sim_param", param_type=ParamType.INT, stages={Stage.SIMULATION}))

        sim_params = reg.get_params_by_stage(Stage.SIMULATION)

        self.assertIn("common_param", sim_params)
        self.assertIn("sim_param", sim_params)

    def test_get_params_by_stage_excludes_other_stages(self):
        """get_params_by_stage should not include params from other stages."""
        reg = ParamRegistry()
        reg.register(ParamDef(name="pre_param", param_type=ParamType.INT, stages={Stage.PRE_PROCESS}))
        reg.register(ParamDef(name="sim_param", param_type=ParamType.INT, stages={Stage.SIMULATION}))

        sim_params = reg.get_params_by_stage(Stage.SIMULATION)

        self.assertNotIn("pre_param", sim_params)
        self.assertIn("sim_param", sim_params)


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
            reg.register(ParamDef(name="test", param_type=ParamType.INT, stages={Stage.COMMON}))

    def test_all_params_readonly_after_freeze(self):
        """all_params should be read-only after freeze."""
        reg = ParamRegistry()
        reg.register(ParamDef(name="test", param_type=ParamType.INT, stages={Stage.COMMON}))
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
        from ..params import REGISTRY, RegistryFrozenError

        with self.assertRaises(RegistryFrozenError):
            REGISTRY.register(
                ParamDef(name="injected", param_type=ParamType.INT, stages={Stage.COMMON})
            )


if __name__ == "__main__":
    unittest.main()
