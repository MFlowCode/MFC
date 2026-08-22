"""
Unit tests for case_validator.py constraint checks.

These cover constraints that are enforced only in Python: the Fortran
m_checker*.fpp counterparts were removed, so a check that silently stops
firing would otherwise go unnoticed (validating the example cases only
exercises configurations that are meant to pass).
"""

import unittest

from .case_validator import CaseConstraintError, CaseValidator

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

# Two-fluid variants, which alt_soundspeed requires (the Kapila K coefficient is a
# two-fluid closure).
TWO_FLUID = {
    **BASE,
    "num_fluids": 2,
    "fluid_pp(2)%gamma": 0.4,
    "fluid_pp(2)%pi_inf": 0.0,
    "patch_icpp(1)%alpha_rho(2)": 0.0,
    "patch_icpp(1)%alpha(2)": 0.0,
}

# 2D axisymmetric (p = 0) and 3D cylindrical (p > 0, odd) two-fluid cases.
CYL_2D = {
    **BASE_2D,
    **{k: v for k, v in TWO_FLUID.items() if k not in BASE},
    "num_fluids": 2,
    "cyl_coord": "T",
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "y_domain%beg": 0.0,
}

CYL_3D = {
    **CYL_2D,
    "p": 49,
    "bc_y%beg": -14,
    "bc_z%beg": -1,
    "bc_z%end": -1,
    "z_domain%beg": 0.0,
    "z_domain%end": 6.28,
    "patch_icpp(1)%z_centroid": 3.14,
    "patch_icpp(1)%length_z": 6.28,
    "patch_icpp(1)%vel(3)": 0.0,
}


class ConstraintTestCase(unittest.TestCase):
    """Base class providing assertions over simulation-stage validation."""

    def errors_for(self, params) -> str:
        """Return the constraint violations for params, or "" if it validates.

        Only CaseConstraintError is caught. Anything else -- a KeyError or
        TypeError from a regression in the validator -- propagates and fails
        the test, rather than being stringified into something an assertion
        might accept.
        """
        validator = CaseValidator(dict(params))
        try:
            validator.validate("simulation")
        except CaseConstraintError as exc:
            return str(exc)
        return ""

    def assertRejects(self, params, expected: str):
        """params must fail validation, and expected must name the reason."""
        errors = self.errors_for(params)
        self.assertNotEqual(errors, "", f"expected validation to fail with {expected!r}, but the case was accepted")
        self.assertIn(expected, errors)

    def assertAccepts(self, params):
        """params must validate cleanly -- no violations at all.

        Asserting full validity rather than the absence of one message means a
        fixture that breaks for an unrelated reason fails here instead of
        silently satisfying a narrower check.
        """
        self.assertEqual(self.errors_for(params), "")


class TestImmersedBoundaryFlags(ConstraintTestCase):
    MSG = "many_ib_patch_parallelism requires ib"

    def test_requires_ib(self):
        self.assertRejects({**BASE, "many_ib_patch_parallelism": "T"}, self.MSG)

    def test_not_tripped_when_disabled(self):
        self.assertAccepts(BASE)


class TestBodyForceSpatialSupport(ConstraintTestCase):
    MSG = "bf_spatial_support is implemented for 2D only"

    def test_rejects_1d(self):
        self.assertRejects({**BASE, "bf_spatial_support": "T"}, self.MSG)

    def test_rejects_3d(self):
        self.assertRejects({**BASE_2D, "p": 50, "bf_spatial_support": "T"}, self.MSG)

    def test_accepts_2d(self):
        self.assertAccepts({**BASE_2D, "bf_spatial_support": "T"})


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
        self.assertAccepts(CHEMISTRY)

    def test_accepts_valid_adaptive_substepping(self):
        params = {**CHEMISTRY, "chem_params%adap_substeps": "T", "chem_params%reaction_substeps": 2, "chem_params%reaction_substeps_max": 8}
        self.assertAccepts(params)

    def test_accepts_igr_without_substepping(self):
        self.assertAccepts({**CHEMISTRY, "igr": "T", "chem_params%reaction_substeps": 0})


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

    def test_rejects_unset_fluid2_eos(self):
        """Unset gamma/pi_inf defaults to dflt_real in the solver, so fluid 2 would
        carry a negative stiffened-gas EOS. The Fortran caught this by comparing
        against the sentinel; an `is not None` guard would silently pass it."""
        for prop in ("gamma", "pi_inf"):
            params = {k: v for k, v in REACTIVE_BURN.items() if k != f"fluid_pp(2)%{prop}"}
            self.assertRejects(params, f"both fluid_pp(1)%{prop} and fluid_pp(2)%{prop} to be set")

    def test_rejects_unset_num_fluids(self):
        params = {k: v for k, v in REACTIVE_BURN.items() if k != "num_fluids"}
        self.assertRejects(params, "reactive_burn requires num_fluids = 2")

    def test_rejects_unset_model_eqns(self):
        params = {k: v for k, v in REACTIVE_BURN.items() if k != "model_eqns"}
        self.assertRejects(params, "reactive_burn requires model_eqns = 2 or 3")

    def test_accepts_valid_configuration(self):
        self.assertAccepts(REACTIVE_BURN)


class TestPhaseChangeFluidPairing(ConstraintTestCase):
    MSG = "phase change requires num_fluids >= 2 (liquid = 1, vapor = 2)"

    def test_rejects_single_fluid(self):
        self.assertRejects({**BASE, "relax": "T", "relax_model": 5}, self.MSG)

    def test_rejects_unset_num_fluids(self):
        params = {k: v for k, v in TWO_FLUID.items() if k != "num_fluids"}
        self.assertRejects({**params, "relax": "T", "relax_model": 5}, self.MSG)

    def test_accepts_two_fluids(self):
        self.assertAccepts({**TWO_FLUID, "relax": "T", "relax_model": 5})

    def test_accepts_three_fluids(self):
        params = {
            **TWO_FLUID,
            "num_fluids": 3,
            "fluid_pp(3)%gamma": 0.4,
            "fluid_pp(3)%pi_inf": 0.0,
            "patch_icpp(1)%alpha_rho(3)": 0.0,
            "patch_icpp(1)%alpha(3)": 0.0,
            "relax": "T",
            "relax_model": 5,
        }
        self.assertAccepts(params)

    def test_not_checked_when_disabled(self):
        self.assertAccepts({**BASE, "relax": "F"})


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
        """ENABLED sets d = 1, 2 only. The Fortran loops d = 1, num_dims, so a 2D
        case needs no z components -- accepting it proves turb_pos(1,3) and
        synth_L(1,3) are not required."""
        self.assertAccepts(self.ENABLED)

    def test_not_checked_when_disabled(self):
        self.assertAccepts(BASE_2D)


class TestTimeStep(ConstraintTestCase):
    """The Fortran checked `if (.not. cfl_dt) dt <= 0`, which fired both on an
    explicitly bad dt and on the dflt_real sentinel left by an unset one."""

    def test_rejects_negative_dt(self):
        self.assertRejects({**BASE, "dt": -1.0}, "dt must be positive")

    def test_rejects_zero_dt(self):
        self.assertRejects({**BASE, "dt": 0.0}, "dt must be positive")

    def test_rejects_unset_dt_under_fixed_stepping(self):
        params = {k: v for k, v in BASE.items() if k != "dt"}
        self.assertRejects(params, "dt must be set when using fixed time stepping")

    def test_rejects_unset_dt_with_adap_dt(self):
        """adap_dt is not an exemption -- it uses dt as its initial value, and the
        Fortran aborted on the sentinel regardless of adap_dt."""
        params = {k: v for k, v in BASE.items() if k != "dt"}
        params.update({"adap_dt": "T", "bubbles_euler": "T", "polytropic": "T", "adv_n": "T", "nb": 1})
        self.assertRejects(params, "dt must be set when using fixed time stepping")

    def test_accepts_cfl_adap_dt_without_dt(self):
        """CFL-driven stepping genuinely needs no dt (e.g. 2D_lagrange_rising_bubble)."""
        params = {k: v for k, v in BASE.items() if k not in ("dt", "t_step_start", "t_step_stop", "t_step_save")}
        params.update({"cfl_adap_dt": "T", "cfl_target": 0.5, "t_stop": 1.0, "t_save": 0.1, "n_start": 0})
        self.assertAccepts(params)

    def test_accepts_positive_dt(self):
        self.assertAccepts(BASE)


class TestHllMethodTwoGeometry(ConstraintTestCase):
    MSG = "HLL Method 2 is not supported for 3D cylindrical geometry"

    def test_rejects_3d_cylindrical(self):
        self.assertRejects({**CYL_3D, "riemann_solver": 1, "hll_u_interface": "T"}, self.MSG)

    def test_accepts_2d_axisymmetric(self):
        self.assertAccepts({**CYL_2D, "riemann_solver": 1, "hll_u_interface": "T"})

    def test_accepts_3d_cylindrical_without_method_two(self):
        self.assertAccepts({**CYL_3D, "riemann_solver": 1})


class TestAltSoundspeedGeometry(ConstraintTestCase):
    """alt_soundspeed with HLL has no cylindrical geometric-source treatment."""

    AXISYM_MSG = "alt_soundspeed with HLL Method 1 is not supported for 2D axisymmetric geometry"
    CYL_3D_MSG = "alt_soundspeed with HLL is not currently supported for 3D cylindrical geometry"

    def test_rejects_hll_method_one_2d_axisymmetric(self):
        self.assertRejects({**CYL_2D, "riemann_solver": 1, "alt_soundspeed": "T"}, self.AXISYM_MSG)

    def test_accepts_hll_method_two_2d_axisymmetric(self):
        """Method 2 carries the shared interface velocity the source term needs."""
        self.assertAccepts({**CYL_2D, "riemann_solver": 1, "hll_u_interface": "T", "alt_soundspeed": "T"})

    def test_rejects_hll_3d_cylindrical(self):
        self.assertRejects({**CYL_3D, "riemann_solver": 1, "hll_u_interface": "T", "alt_soundspeed": "T"}, self.CYL_3D_MSG)

    def test_accepts_hllc_cartesian(self):
        self.assertAccepts({**TWO_FLUID, "riemann_solver": 2, "alt_soundspeed": "T"})


class TestAltSoundspeedHlld(ConstraintTestCase):
    MSG = "alt_soundspeed with HLLD requires hypoelasticity = T"

    def test_rejects_hlld_without_hypoelasticity(self):
        self.assertRejects({**TWO_FLUID, "riemann_solver": 4, "alt_soundspeed": "T"}, self.MSG)

    def test_not_tripped_without_alt_soundspeed(self):
        self.assertNotIn(self.MSG, self.errors_for({**TWO_FLUID, "riemann_solver": 4}))


if __name__ == "__main__":
    unittest.main()
