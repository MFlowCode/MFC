"""
Unit tests for params/registry.py module.

Tests registry functionality, freezing, tag queries, and indexed families.
"""

import unittest

from ..params.registry import IndexedFamily, ParamRegistry, RegistryFrozenError
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

        self.assertEqual(reg.all_params["test"].tags, {"mhd", "physics"})

    def test_register_type_mismatch_raises(self):
        """Registering same param with different type should raise."""
        reg = ParamRegistry()
        reg.register(ParamDef(name="test", param_type=ParamType.INT))

        with self.assertRaises(ValueError) as ctx:
            reg.register(ParamDef(name="test", param_type=ParamType.REAL))

        self.assertIn("existing type is", str(ctx.exception))
        self.assertIn("new type is", str(ctx.exception))

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


class TestIndexedFamily(unittest.TestCase):
    """Tests for indexed family registration and resolution."""

    def _make_registry_with_family(self, max_index=None):
        """Helper: registry with one indexed family."""
        reg = ParamRegistry()
        reg.register_family(
            IndexedFamily(
                base_name="thing",
                attrs={
                    "geom": (ParamType.INT, {"tag1"}),
                    "vel(1)": (ParamType.REAL, {"tag1"}),
                },
                tags={"tag1"},
                max_index=max_index,
            )
        )
        reg.freeze()
        return reg

    def test_valid_family_param_is_known(self):
        """Valid indexed family param should be recognized."""
        reg = self._make_registry_with_family()
        self.assertTrue(reg.is_known_param("thing(1)%geom"))
        self.assertTrue(reg.is_known_param("thing(999)%vel(1)"))

    def test_bogus_attr_rejected(self):
        """Unknown attribute on a valid family should be rejected."""
        reg = self._make_registry_with_family()
        self.assertFalse(reg.is_known_param("thing(1)%nonexistent"))

    def test_zero_index_rejected(self):
        """Index 0 should be rejected (1-based indexing)."""
        reg = self._make_registry_with_family()
        self.assertFalse(reg.is_known_param("thing(0)%geom"))

    def test_unknown_family_rejected(self):
        """Unknown family base name should be rejected."""
        reg = self._make_registry_with_family()
        self.assertFalse(reg.is_known_param("unknown(1)%geom"))

    def test_max_index_enforced(self):
        """Indices beyond max_index should be rejected."""
        reg = self._make_registry_with_family(max_index=5)
        self.assertTrue(reg.is_known_param("thing(5)%geom"))
        self.assertFalse(reg.is_known_param("thing(6)%geom"))

    def test_unlimited_index_when_no_max(self):
        """Without max_index, arbitrarily large indices are valid."""
        reg = self._make_registry_with_family(max_index=None)
        self.assertTrue(reg.is_known_param("thing(999999)%geom"))

    def test_get_param_def_returns_correct_type(self):
        """get_param_def should return correct ParamType for family params."""
        reg = self._make_registry_with_family()
        pdef = reg.get_param_def("thing(3)%geom")
        self.assertIsNotNone(pdef)
        self.assertEqual(pdef.param_type, ParamType.INT)
        self.assertEqual(pdef.name, "thing(3)%geom")

    def test_get_param_def_none_for_invalid(self):
        """get_param_def should return None for invalid family params."""
        reg = self._make_registry_with_family()
        self.assertIsNone(reg.get_param_def("thing(1)%bogus"))
        self.assertIsNone(reg.get_param_def("unknown(1)%geom"))

    def test_all_params_iteration_bounded(self):
        """Iterating all_params should yield scalars + one example per family attr."""
        reg = ParamRegistry()
        reg.register(ParamDef(name="scalar1", param_type=ParamType.INT))
        reg.register(ParamDef(name="scalar2", param_type=ParamType.REAL))
        reg.register_family(
            IndexedFamily(
                base_name="thing",
                attrs={
                    "geom": (ParamType.INT, {"tag1"}),
                    "vel(1)": (ParamType.REAL, {"tag1"}),
                },
                tags={"tag1"},
            )
        )
        reg.freeze()

        keys = list(reg.all_params)
        # 2 scalars + 2 family attrs (one example each at index=1)
        self.assertEqual(len(reg.all_params), 4)
        self.assertEqual(len(keys), 4)
        self.assertIn("scalar1", keys)
        self.assertIn("scalar2", keys)
        # Family examples use index=1
        self.assertIn("thing(1)%geom", keys)
        self.assertIn("thing(1)%vel(1)", keys)
        # Arbitrary indices should NOT appear in iteration
        self.assertNotIn("thing(42)%geom", keys)

    def test_all_params_contains_family(self):
        """all_params mapping should resolve family params via __contains__."""
        reg = self._make_registry_with_family()
        self.assertIn("thing(42)%geom", reg.all_params)
        self.assertNotIn("thing(1)%bogus", reg.all_params)

    def test_all_params_getitem_family(self):
        """all_params[family_param] should return a ParamDef."""
        reg = self._make_registry_with_family()
        pdef = reg.all_params["thing(7)%vel(1)"]
        self.assertEqual(pdef.param_type, ParamType.REAL)

    def test_all_params_get_family(self):
        """all_params.get() should resolve family params and respect defaults."""
        reg = self._make_registry_with_family()
        pdef = reg.all_params.get("thing(5)%geom")
        self.assertIsNotNone(pdef)
        self.assertEqual(pdef.param_type, ParamType.INT)
        self.assertEqual(reg.all_params.get("thing(1)%nonexistent", "default"), "default")

    def test_all_params_getitem_raises_for_invalid(self):
        """all_params[invalid] should raise KeyError."""
        reg = self._make_registry_with_family()
        with self.assertRaises(KeyError):
            _ = reg.all_params["thing(1)%bogus"]

    def test_register_family_after_freeze_raises(self):
        """Registering a family after freeze should raise."""
        reg = ParamRegistry()
        reg.freeze()
        with self.assertRaises(RegistryFrozenError):
            reg.register_family(IndexedFamily(base_name="late", attrs={}, tags=set()))

    def test_invalid_base_name_raises(self):
        """IndexedFamily with empty or invalid base_name should raise."""
        with self.assertRaises(ValueError):
            IndexedFamily(base_name="")
        with self.assertRaises(ValueError):
            IndexedFamily(base_name="123invalid")

    def test_invalid_max_index_raises(self):
        """IndexedFamily with max_index < 1 should raise."""
        with self.assertRaises(ValueError):
            IndexedFamily(base_name="test", max_index=0)
        with self.assertRaises(ValueError):
            IndexedFamily(base_name="test", max_index=-5)

    def test_all_params_before_freeze_with_families_raises(self):
        """Accessing all_params before freeze with families should raise."""
        reg = ParamRegistry()
        reg.register_family(
            IndexedFamily(
                base_name="fam",
                attrs={"x": (ParamType.INT, {"t"})},
                tags={"t"},
            )
        )
        with self.assertRaises(RuntimeError):
            _ = reg.all_params

    def test_family_params_by_tag(self):
        """get_params_by_tag should include family example params."""
        reg = self._make_registry_with_family()
        tagged = reg.get_params_by_tag("tag1")
        self.assertGreater(len(tagged), 0)


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
            REGISTRY.register(ParamDef(name="injected", param_type=ParamType.INT))


if __name__ == "__main__":
    unittest.main()
