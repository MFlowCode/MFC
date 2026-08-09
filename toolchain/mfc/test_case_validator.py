"""
Unit tests for case_validator.py constraint checks.

These cover constraints that are enforced only in Python: the Fortran
m_checker*.fpp counterparts were removed, so a check that silently stops
firing would otherwise go unnoticed (validating the example cases only
exercises configurations that are meant to pass).
"""

import unittest

from .case_validator import CaseValidator

# A minimal 1D case that passes simulation validation.
BASE = {
    "m": 50,
    "n": 0,
    "p": 0,
    "model_eqns": 2,
    "num_fluids": 1,
    "num_patches": 1,
    "t_step_start": 0,
    "t_step_stop": 100,
    "t_step_save": 10,
    "dt": 1e-6,
    "weno_order": 5,
    "weno_eps": 1e-6,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": 1.0,
    "patch_icpp(1)%alpha_rho(1)": 1.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 0.4,
    "fluid_pp(1)%pi_inf": 0.0,
}

BASE_2D = {
    **BASE,
    "n": 50,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "y_domain%beg": 0.0,
    "y_domain%end": 1.0,
    "patch_icpp(1)%y_centroid": 0.5,
    "patch_icpp(1)%length_y": 1.0,
    "patch_icpp(1)%vel(2)": 0.0,
}

# A reactive-burn case satisfying every rburn constraint.
REACTIVE_BURN = {
    **BASE,
    "num_fluids": 2,
    "reactive_burn": "T",
    "rburn%k": 1.0e6,
    "rburn%pign": 1.0e8,
    "rburn%pref": 1.0e8,
    "rburn%n": 1.0,
    "rburn%ta": 0.0,
    "fluid_pp(1)%qv": 1.0e6,
    "fluid_pp(2)%gamma": 0.4,
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(2)%qv": 0.0,
    "patch_icpp(1)%alpha_rho(2)": 0.0,
    "patch_icpp(1)%alpha(2)": 0.0,
}

CHEMISTRY = {**BASE, "chemistry": "T", "cantera_file": "h2o2.yaml"}


class ConstraintTestCase(unittest.TestCase):
    """Base class providing assertions over simulation-stage validation."""

    def errors_for(self, params) -> str:
        """Return all simulation-stage validation errors for params, joined."""
        validator = CaseValidator(dict(params))
        try:
            validator.validate("simulation")
        except Exception as exc:  # CaseConstraintError
            return str(exc)
        return ""

    def assertRejects(self, params, expected: str):
        """params must fail validation with expected in the message."""
        errors = self.errors_for(params)
        self.assertIn(expected, errors)

    def assertAccepts(self, params, unexpected: str):
        """params must not trip the check identified by unexpected."""
        errors = self.errors_for(params)
        self.assertNotIn(unexpected, errors)


class TestImmersedBoundaryFlags(ConstraintTestCase):
    MSG = "many_ib_patch_parallelism requires ib"

    def test_requires_ib(self):
        self.assertRejects({**BASE, "many_ib_patch_parallelism": "T"}, self.MSG)

    def test_not_tripped_when_disabled(self):
        self.assertAccepts(BASE, self.MSG)


class TestBodyForceSpatialSupport(ConstraintTestCase):
    MSG = "bf_spatial_support is implemented for 2D only"

    def test_rejects_1d(self):
        self.assertRejects({**BASE, "bf_spatial_support": "T"}, self.MSG)

    def test_rejects_3d(self):
        self.assertRejects({**BASE_2D, "p": 50, "bf_spatial_support": "T"}, self.MSG)

    def test_accepts_2d(self):
        self.assertAccepts({**BASE_2D, "bf_spatial_support": "T"}, self.MSG)


class TestChemistrySubstepping(ConstraintTestCase):
    def test_rejects_negative_substeps(self):
        self.assertRejects({**CHEMISTRY, "chem_params%reaction_substeps": -1}, "reaction_substeps must be >= 0")

    def test_rejects_substepping_with_igr(self):
        self.assertRejects({**CHEMISTRY, "igr": "T", "chem_params%reaction_substeps": 2}, "not supported with igr")

    def test_rejects_adaptive_without_substeps(self):
        """adap_substeps with reaction_substeps unset: the Fortran default is 0, below the floor of 1."""
        self.assertRejects({**CHEMISTRY, "chem_params%adap_substeps": "T"}, "requires reaction_substeps >= 1")

    def test_rejects_adaptive_with_zero_substeps(self):
        self.assertRejects({**CHEMISTRY, "chem_params%adap_substeps": "T", "chem_params%reaction_substeps": 0}, "requires reaction_substeps >= 1")

    def test_rejects_max_below_floor(self):
        params = {**CHEMISTRY, "chem_params%adap_substeps": "T", "chem_params%reaction_substeps": 5, "chem_params%reaction_substeps_max": 2}
        self.assertRejects(params, "reaction_substeps_max must be >=")

    def test_rejects_max_unset_below_floor(self):
        """reaction_substeps_max unset defaults to 0 in the Fortran, below a floor of 5."""
        self.assertRejects({**CHEMISTRY, "chem_params%adap_substeps": "T", "chem_params%reaction_substeps": 5}, "reaction_substeps_max must be >=")

    def test_accepts_no_substepping(self):
        self.assertAccepts(CHEMISTRY, "reaction_substeps")

    def test_accepts_valid_adaptive_substepping(self):
        params = {**CHEMISTRY, "chem_params%adap_substeps": "T", "chem_params%reaction_substeps": 2, "chem_params%reaction_substeps_max": 8}
        self.assertAccepts(params, "reaction_substeps")

    def test_accepts_igr_without_substepping(self):
        self.assertAccepts({**CHEMISTRY, "igr": "T", "chem_params%reaction_substeps": 0}, "not supported with igr")


class TestReactiveBurnFluidPairing(ConstraintTestCase):
    def test_rejects_wrong_num_fluids(self):
        self.assertRejects({**REACTIVE_BURN, "num_fluids": 3}, "reactive_burn requires num_fluids = 2")

    def test_rejects_gamma_mismatch(self):
        self.assertRejects({**REACTIVE_BURN, "fluid_pp(2)%gamma": 0.5}, "fluid_pp(1)%gamma == fluid_pp(2)%gamma")

    def test_rejects_pi_inf_mismatch(self):
        self.assertRejects({**REACTIVE_BURN, "fluid_pp(2)%pi_inf": 1.0e5}, "fluid_pp(1)%pi_inf == fluid_pp(2)%pi_inf")

    def test_rejects_equal_qv(self):
        self.assertRejects({**REACTIVE_BURN, "fluid_pp(1)%qv": 0.0}, "fluid_pp(1)%qv > fluid_pp(2)%qv")

    def test_rejects_inverted_qv(self):
        self.assertRejects({**REACTIVE_BURN, "fluid_pp(1)%qv": 0.0, "fluid_pp(2)%qv": 1.0e6}, "fluid_pp(1)%qv > fluid_pp(2)%qv")

    def test_rejects_unset_qv(self):
        """qv defaults to 0 in the Fortran, so leaving both unset means no energy release."""
        params = {k: v for k, v in REACTIVE_BURN.items() if not k.endswith("%qv")}
        self.assertRejects(params, "fluid_pp(1)%qv > fluid_pp(2)%qv")

    def test_accepts_valid_configuration(self):
        self.assertEqual(self.errors_for(REACTIVE_BURN), "")


class TestSyntheticTurbulence(ConstraintTestCase):
    """A 2D case with one fully specified forcing zone."""

    ENABLED = {
        **BASE_2D,
        "synthetic_turbulence": "T",
        "num_turbulent_sources": 1,
        "turb_pos(1,1)": 0.5,
        "turb_pos(1,2)": 0.5,
        "synth_L(1,1)": 1.0,
        "synth_L(1,2)": 1.0,
    }

    def test_rejects_zero_sources(self):
        self.assertRejects({**self.ENABLED, "num_turbulent_sources": 0}, "num_turbulent_sources must be > 0")

    def test_rejects_unset_sources(self):
        params = {k: v for k, v in self.ENABLED.items() if k != "num_turbulent_sources"}
        self.assertRejects(params, "num_turbulent_sources must be > 0")

    def test_rejects_missing_position(self):
        params = {k: v for k, v in self.ENABLED.items() if k != "turb_pos(1,2)"}
        self.assertRejects(params, "turb_pos(1,2) must be specified")

    def test_rejects_missing_extent(self):
        params = {k: v for k, v in self.ENABLED.items() if k != "synth_L(1,2)"}
        self.assertRejects(params, "synth_L(1,2) must be positive")

    def test_rejects_nonpositive_extent(self):
        self.assertRejects({**self.ENABLED, "synth_L(1,2)": 0.0}, "synth_L(1,2) must be positive")

    def test_accepts_fully_specified_zone(self):
        self.assertEqual(self.errors_for(self.ENABLED), "")

    def test_third_dimension_not_required_in_2d(self):
        """The Fortran loops d = 1, num_dims, so a 2D case needs no z components."""
        self.assertAccepts(self.ENABLED, "turb_pos(1,3)")
        self.assertAccepts(self.ENABLED, "synth_L(1,3)")

    def test_not_checked_when_disabled(self):
        self.assertAccepts(BASE_2D, "num_turbulent_sources")


class TestTimeStepPositivity(ConstraintTestCase):
    MSG = "dt must be positive"

    def test_rejects_negative_dt(self):
        self.assertRejects({**BASE, "dt": -1.0}, self.MSG)

    def test_rejects_zero_dt(self):
        self.assertRejects({**BASE, "dt": 0.0}, self.MSG)

    def test_accepts_positive_dt(self):
        self.assertAccepts(BASE, self.MSG)


if __name__ == "__main__":
    unittest.main()
