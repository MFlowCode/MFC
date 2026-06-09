"""Tests for named parameter values (the "names" key in CONSTRAINTS)."""

import re


def test_names_are_wellformed_and_cover_choices():
    from mfc.params.definitions import CONSTRAINTS

    name_re = re.compile(r"^[a-z0-9][a-z0-9_]*$")
    for param, c in CONSTRAINTS.items():
        names = c.get("names")
        if names is None:
            continue
        assert isinstance(names, dict), param
        for name, value in names.items():
            assert name_re.match(name), f"{param}: bad name {name!r}"
            assert isinstance(value, int), f"{param}: {name} -> {value!r} is not int"
        assert len(set(names.values())) == len(names), f"{param}: duplicate values in names"
        if "choices" in c:
            assert set(names.values()) == set(c["choices"]), f"{param}: names do not cover choices"
    assert any("names" in c for c in CONSTRAINTS.values()), "no CONSTRAINTS entry defines names"


def test_riemann_solver_names():
    from mfc.params.definitions import CONSTRAINTS

    assert CONSTRAINTS["riemann_solver"]["names"] == {"hll": 1, "hllc": 2, "hlld": 4, "lax_friedrichs": 5}


def test_invalid_names_rejected():
    import pytest

    from mfc.params.definitions import _validate_constraint

    with pytest.raises(ValueError):
        _validate_constraint("x", {"choices": [1, 2], "names": {"BadName": 1, "ok": 2}})
    with pytest.raises(ValueError):
        _validate_constraint("x", {"choices": [1, 2], "names": {"a": 1, "b": 3}})
    with pytest.raises(ValueError):
        _validate_constraint("x", {"choices": [1, 2], "names": ["a", "b"]})
    with pytest.raises(ValueError):
        _validate_constraint("x", {"choices": [1, 2], "names": {1: "a", 2: "b"}})
