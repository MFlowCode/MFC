"""
Tests for CaseValidator.flag() and ast_analyzer._build_local_param_map fix.
"""

import ast
import unittest

from ..case_validator import CaseValidator
from ..params.ast_analyzer import CaseValidatorAnalyzer


def _param_map(source: str) -> dict:
    """Parse a single method body and return its local-param map."""
    wrapped = "class _V:\n" + "\n".join("    " + l for l in source.splitlines())
    func = ast.parse(wrapped).body[0].body[0]
    a = CaseValidatorAnalyzer.__new__(CaseValidatorAnalyzer)
    return a._build_local_param_map(func)


class TestFlag(unittest.TestCase):

    def test_true_when_set_to_T(self):
        self.assertTrue(CaseValidator({"chemistry": "T"}).flag("chemistry"))

    def test_false_when_set_to_F(self):
        self.assertFalse(CaseValidator({"chemistry": "F"}).flag("chemistry"))

    def test_false_when_absent(self):
        self.assertFalse(CaseValidator({}).flag("chemistry"))

    def test_only_uppercase_T_is_truthy(self):
        """Lowercase, integers, and booleans must not count as set."""
        for val in ("t", "True", 1, True, None):
            self.assertFalse(CaseValidator({"x": val}).flag("x"), f"Expected False for {val!r}")

    def test_returns_bool(self):
        self.assertIsInstance(CaseValidator({"x": "T"}).flag("x"), bool)


class TestASTCoupling(unittest.TestCase):

    def test_flag_call_recorded_by_analyzer(self):
        """self.flag('x') must register 'x' in the param map, same as self.get('x','F')=='T'."""
        m = _param_map("def check(self):\n    chemistry = self.flag('chemistry')")
        self.assertEqual(m.get("chemistry"), "chemistry")

    def test_get_call_still_recorded(self):
        """Existing self.get() pattern must not regress."""
        m = _param_map("def check(self):\n    d = self.get('chem_params%diffusion', 'F') == 'T'")
        self.assertEqual(m.get("d"), "chem_params%diffusion")

    def test_old_matcher_would_have_missed_flag(self):
        """Demonstrate the pre-fix bug: attr=='get' silently dropped flag() calls."""
        source = "def check(self):\n    chemistry = self.flag('chemistry')"
        wrapped = "class _V:\n" + "\n".join("    " + l for l in source.splitlines())
        func = ast.parse(wrapped).body[0].body[0]

        broken = {}
        for node in ast.walk(func):
            if isinstance(node, ast.Assign):
                v = node.value.left if isinstance(node.value, ast.Compare) else node.value
                if isinstance(v, ast.Call) and isinstance(v.func, ast.Attribute):
                    if v.func.value.id == "self" and v.func.attr == "get":  # old check
                        for t in node.targets:
                            if isinstance(t, ast.Name):
                                broken[t.id] = v.args[0].value

        self.assertNotIn("chemistry", broken, "Old matcher should miss flag() — that was the bug")

        a = CaseValidatorAnalyzer.__new__(CaseValidatorAnalyzer)
        fixed = a._build_local_param_map(func)
        self.assertIn("chemistry", fixed, "Fixed matcher must capture flag() calls")


if __name__ == "__main__":
    unittest.main()