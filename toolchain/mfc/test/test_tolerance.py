"""override_tol and override_tol_single: which tolerance a golden comparison uses in double and in single precision."""

import pytest

from mfc import state
from mfc.test.case import TestCase


@pytest.fixture
def precision(monkeypatch):
    def set_single(single):
        monkeypatch.setattr(state, "gARG", {**state.gARG, "single": single})

    return set_single


def test_double_uses_override_tol(precision):
    precision(False)
    assert TestCase("t", {}, override_tol=1e-10, override_tol_single=1e-3).compute_tolerance() == 1e-10


def test_single_uses_override_tol_single(precision):
    precision(True)
    assert TestCase("t", {}, override_tol=1e-10, override_tol_single=1e-3).compute_tolerance() == 1e-3


def test_single_floor_without_override_tol_single(precision):
    precision(True)
    assert TestCase("t", {}, override_tol=1e-10).compute_tolerance() == pytest.approx(1e-4)


def test_single_floor_bounds_a_tight_override_tol_single(precision):
    precision(True)
    assert TestCase("t", {}, override_tol=1e-10, override_tol_single=1e-6).compute_tolerance() == pytest.approx(1e-4)
