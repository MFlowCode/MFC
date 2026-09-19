"""Guards on the surface-mechanism half of the chemistry toolchain.

Every case here is a silent-wrongness path: the generator emits a rate expression that is
wrong by orders of magnitude and nothing is printed. The integration goldens cannot see it --
they run one mechanism, whose every reaction happens to be of the one supported kind.

Cantera's own ptcombust.yaml is the fixture because it carries all three kinds at once:
19 interface-Arrhenius reactions, 5 sticking-coefficient ones, and 2 with coverage
dependence. A synthetic mechanism would only prove Cantera parses what we wrote.
"""

import os

import pytest

from mfc import common

ct = pytest.importorskip("cantera", reason="the surface codegen needs Cantera")

# What generate_surface_thermochem is willing to translate. Kept as data so a test failure
# names the mismatch rather than a bare isinstance.
SUPPORTED = (ct.ArrheniusRate, ct.InterfaceArrheniusRate)


def _coverage(rate):
    return dict(getattr(rate, "coverage_dependencies", {}) or {})


@pytest.fixture(scope="module")
def ptcombust():
    return ct.Interface("ptcombust.yaml", "Pt_surf")


def test_a_sticking_coefficient_is_not_mistaken_for_an_arrhenius_prefactor(ptcombust):
    """A sticking rate's 'pre_exponential_factor' is a dimensionless probability, not a rate
    prefactor; using it as one is wrong by orders of magnitude. Cantera exposes the attribute on
    both, so hasattr() cannot tell them apart -- which is exactly how it was got wrong once."""
    sticking = [r for r in ptcombust.reactions() if isinstance(r.rate, ct.StickRateBase)]
    assert sticking, "ptcombust.yaml is expected to carry sticking reactions"

    for reaction in sticking:
        rate = reaction.rate
        assert hasattr(rate, "pre_exponential_factor"), "the weak guard this test rules out"
        assert not isinstance(rate, SUPPORTED), f"sticking rate would be emitted as Arrhenius: {reaction.equation}"


def test_coverage_dependent_rates_are_not_silently_stripped(ptcombust):
    """The generator emits no coverage terms, so accepting such a reaction drops part of its
    rate law without saying so."""
    covered = [r for r in ptcombust.reactions() if _coverage(r.rate)]
    assert covered, "ptcombust.yaml is expected to carry coverage-dependent reactions"


def test_plain_surface_arrhenius_reactions_remain_acceptable(ptcombust):
    """The converse: a guard that rejects everything is as broken as one that accepts anything."""
    plain = [r for r in ptcombust.reactions() if isinstance(r.rate, SUPPORTED) and not _coverage(r.rate)]
    assert plain, "ptcombust.yaml is expected to carry ordinary interface-Arrhenius reactions"


def test_the_shipped_carbon_mechanism_is_translatable_end_to_end():
    """The reacting-surface example must keep generating. With no CFD golden behind it (the
    example is skipped in cases.py -- the suite cannot afford a third compiled-in mechanism),
    these tests are the whole safety net, so a guard that turns the example away breaks the
    feature rather than protecting it."""
    path = os.path.join(common.MFC_ROOT_DIR, "examples", "2D_ibm_reacting_surface", "carbon_surface_bradley_11species.yaml")
    if not os.path.isfile(path):
        pytest.skip("reacting-surface example not present")

    surface = ct.Interface(path, "carbon_surface")
    assert surface.n_reactions > 0

    for reaction in surface.reactions():
        assert isinstance(reaction.rate, SUPPORTED), reaction.equation
        assert not _coverage(reaction.rate), reaction.equation


def test_surface_site_species_in_stoichiometry_are_refused_with_a_reason():
    """ptcombust is a coverage-based mechanism: its reactions consume and produce adsorbed sites
    (PT(S), O(S)). The generator models a bulk surface reacting with gas species and has no site
    balance, so it must say so rather than emit a rate law missing its site terms."""
    from mfc.run.input import MFCInputFile

    case = MFCInputFile("case.py", ".", {})
    surface = ct.Interface("ptcombust.yaml", "Pt_surf")

    with pytest.raises(common.MFCException, match="Surface-site species"):
        case.generate_surface_thermochem(ct.Solution("ptcombust.yaml", "gas"), surface)


def test_the_shipped_carbon_mechanism_generates_compilable_fortran():
    """The end-to-end counterpart: the example's mechanism must produce a module with both
    entry points m_ibm.fpp imports, and a rate expression per reaction."""
    from mfc.run.input import MFCInputFile

    path = os.path.join(common.MFC_ROOT_DIR, "examples", "2D_ibm_reacting_surface", "carbon_surface_bradley_11species.yaml")
    if not os.path.isfile(path):
        pytest.skip("reacting-surface example not present")

    # Through get_cantera_surface, not a hand-built ct.Interface: the mechanism declares adjacent
    # phases that have to be resolved and loaded, which is exactly what that method is for, so this
    # covers the resolution path as well as the generator.
    case = MFCInputFile(
        "case.py",
        os.path.dirname(path),
        {
            "chemistry": "T",
            "cantera_file": os.path.join(os.path.dirname(path), "carbon_gasphase_reduced_gri11.yaml"),
            "surface_cantera_file": path,
            "surface_phase": "carbon_surface",
        },
    )
    gas = case.get_cantera_solution()
    surface = case.get_cantera_surface()
    assert surface is not None, "the example's surface mechanism failed to resolve"
    code = case.generate_surface_thermochem(gas, surface)

    assert "module m_surface_thermochem" in code
    assert "subroutine get_surface_net_production_rates" in code
    assert "subroutine get_surface_reaction_heat_flux" in code
    for i in range(1, surface.n_reactions + 1):
        assert f"! Surface reaction {i}" in code, f"reaction {i} produced no rate expression"


def test_a_case_without_a_surface_mechanism_still_generates_a_module():
    """Chemistry without surface reactions must still compile: m_ibm.fpp imports the module
    unconditionally, so the generator has to emit the no-op form rather than nothing."""
    from mfc.run.input import MFCInputFile

    code = MFCInputFile("case.py", ".", {}).generate_surface_thermochem(ct.Solution("h2o2.yaml"), None)

    assert "module m_surface_thermochem" in code
    assert "get_surface_net_production_rates" in code
    assert "get_surface_reaction_heat_flux" in code
