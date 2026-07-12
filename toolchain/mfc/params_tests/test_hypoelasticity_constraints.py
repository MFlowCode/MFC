"""Regression tests for hypoelasticity feature-compatibility constraints."""

import unittest

from ..case_validator import CaseValidator
from ..params import REGISTRY


class TestHypoelasticityConstraints(unittest.TestCase):
    BASE = {
        "hypoelasticity": "T",
        "model_eqns": 2,
        "riemann_solver": 1,
        "num_fluids": 2,
        "n": 1,
        "p": 0,
    }

    def errors_for(self, **updates):
        validator = CaseValidator({**self.BASE, **updates})
        validator.check_hypoelasticity()
        return validator.errors

    def assert_rejected(self, message, **updates):
        self.assertIn(message, self.errors_for(**updates))

    def test_removed_energy_guard_is_not_registered(self):
        self.assertIsNone(REGISTRY.get_param_def("hypo_energy_guard"))

    def test_mhd_is_prohibited_for_every_hypo_solver(self):
        for solver in (1, 2, 4):
            with self.subTest(solver=solver):
                self.assert_rejected("MHD and hypoelasticity cannot be enabled together", riemann_solver=solver, mhd="T")

    def test_hyperelasticity_is_prohibited(self):
        self.assert_rejected("Hypoelasticity and hyperelasticity cannot be enabled together", hyperelasticity="T")

    def test_euler_bubbles_are_prohibited_for_every_hypo_solver(self):
        for solver in (1, 2, 4):
            with self.subTest(solver=solver):
                self.assert_rejected("Hypoelasticity does not support Euler-Euler bubbles", riemann_solver=solver, bubbles_euler="T")

    def test_continuum_damage_with_alt_soundspeed_is_prohibited(self):
        validator = CaseValidator(
            {
                **self.BASE,
                "hypoelasticity": "F",
                "cont_damage": "T",
                "alt_soundspeed": "T",
                "tau_star": 0.0,
                "cont_damage_s": 2.0,
                "alpha_bar": 1e-4,
            }
        )
        validator.check_continuum_damage()
        self.assertIn("Continuum damage does not support alt_soundspeed", validator.errors)

    def test_hll_and_hllc_allow_more_than_two_fluids_without_alt_soundspeed(self):
        for solver in (1, 2):
            with self.subTest(solver=solver):
                self.assertEqual([], self.errors_for(riemann_solver=solver, num_fluids=3, alt_soundspeed="F"))

    def test_more_than_two_fluids_remain_prohibited_with_alt_soundspeed(self):
        self.assert_rejected("hypoelastic alt_soundspeed requires exactly 2 fluid components", num_fluids=3, alt_soundspeed="T")

    def test_hlld_still_requires_exactly_two_fluids(self):
        self.assert_rejected("HLLD hypoelasticity requires exactly 2 fluid components", riemann_solver=4, num_fluids=3)


if __name__ == "__main__":
    unittest.main()
