"""Tests for the species-bounds assertion the test harness applies to chemistry cases (#2154)."""

from types import SimpleNamespace

import pytest

from mfc.common import MFCException
from mfc.packer.pack import Pack, PackEntry
from mfc.test.test import _assert_species_bounds

N_SPECIES = 2


def _case(chemistry="T"):
    solution = SimpleNamespace(n_species=N_SPECIES)
    return SimpleNamespace(params={"chemistry": chemistry}, to_input_file=lambda: SimpleNamespace(get_cantera_solution=lambda: solution))


def _pack(rho_y, rho=(1.0, 1.0)):
    """Equations 1 (density), 2 (energy) and the species; one cell per entry of rho_y[k]."""
    entries = {1: list(rho), 2: [5.0, 5.0], 3: rho_y[0], 4: rho_y[1]}
    return Pack([PackEntry(f"D/cons.{i}.00.000050.dat", doubles) for i, doubles in entries.items()])


CLEAN = "Species cleanup at step 7: 2 cell cleanings, sum |d rho| =  4.045E-12\n"


def test_mass_fractions_on_the_simplex_pass():
    _assert_species_bounds(_case(), _pack([[0.25, 1.0], [0.75, 0.0]]), "")


def test_a_negative_mass_fraction_fails():
    with pytest.raises(MFCException, match="leaves"):
        _assert_species_bounds(_case(), _pack([[-1.0e-6, 1.0], [1.0, 0.0]]), "")


def test_a_mass_fraction_above_one_fails():
    with pytest.raises(MFCException, match="leaves"):
        _assert_species_bounds(_case(), _pack([[1.0 + 1.0e-6, 1.0], [0.0, 0.0]]), "")


def test_roundoff_scale_cleanup_is_allowed():
    _assert_species_bounds(_case(), _pack([[0.25, 1.0], [0.75, 0.0]]), CLEAN)


def test_a_cleanup_that_moves_real_mass_fails():
    with pytest.raises(MFCException, match="clean-up at step 9"):
        _assert_species_bounds(_case(), _pack([[0.25, 1.0], [0.75, 0.0]]), CLEAN + "Species cleanup at step 9: 40 cell cleanings, sum |d rho| =  3.2E-03\n")


def test_missing_output_is_an_error_not_a_pass():
    with pytest.raises(MFCException, match="no D/cons"):
        _assert_species_bounds(_case(), Pack(), "")


def test_cases_without_chemistry_are_skipped():
    _assert_species_bounds(_case(chemistry="F"), Pack(), "Species cleanup at step 9: 40 cell cleanings, sum |d rho| =  3.2E-03\n")


def test_mass_fractions_are_relative_to_the_density():
    # rho = 2: partial densities up to 2 are fine, and 1.5 + 0.5 sums to it
    _assert_species_bounds(_case(), _pack([[1.5, 0.0], [0.5, 2.0]], rho=(2.0, 2.0)), "")
    with pytest.raises(MFCException, match="leaves"):
        _assert_species_bounds(_case(), _pack([[2.5, 0.0], [0.0, 0.0]], rho=(2.0, 2.0)), "")


def test_a_cell_left_without_mass_fails_even_though_no_density_moved():
    with pytest.raises(MFCException, match="clean-up at step 3"):
        _assert_species_bounds(_case(), _pack([[0.25, 1.0], [0.75, 0.0]]), "Species cleanup at step 3: 1 cell cleanings, sum |d rho| =  0.000E+00\n")


def test_a_nan_cleanup_report_fails():
    with pytest.raises(MFCException, match="clean-up"):
        _assert_species_bounds(_case(), _pack([[0.25, 1.0], [0.75, 0.0]]), "Species cleanup at step 3: 1 cell cleanings, sum |d rho| =       NaN\n")


def test_a_three_digit_exponent_is_read():
    # Fortran ES10.3 prints 1.000+100 for 1e100 and 1.000-100 for 1e-100
    with pytest.raises(MFCException, match="clean-up"):
        _assert_species_bounds(_case(), _pack([[0.25, 1.0], [0.75, 0.0]]), "Species cleanup at step 3: 9 cell cleanings, sum |d rho| =  1.000+100\n")
    _assert_species_bounds(_case(), _pack([[0.25, 1.0], [0.75, 0.0]]), "Species cleanup at step 3: 1 cell cleanings, sum |d rho| =  1.000-100\n")


def test_single_precision_uses_the_looser_limits():
    report = "Species cleanup at step 3: 1 cell cleanings, sum |d rho| =  1.000E-06\n"
    with pytest.raises(MFCException, match="clean-up"):
        _assert_species_bounds(_case(), _pack([[0.25, 1.0], [0.75, 0.0]]), report)
    _assert_species_bounds(_case(), _pack([[0.25, 1.0], [0.75, 0.0]]), report, single=True)
    _assert_species_bounds(_case(), _pack([[0.25, 1.0], [0.75 + 1.0e-6, 0.0]]), "", single=True)
