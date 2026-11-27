"""bench_diff compares simulation grind time and nothing else.

Execution times and pre_process/post_process timings bury the one number a
reviewer needs, so they are neither shown, warned about, nor required.
"""

import pytest
import yaml

from mfc import state
from mfc.bench import diff
from mfc.common import MFCException


def _results(tmp_path, name, cases):
    path = tmp_path / f"{name}.yaml"
    path.write_text(yaml.safe_dump({"metadata": {"invocation": ["bench"], "lock": {}}, "cases": {slug: {"output_summary": s} for slug, s in cases.items()}}))
    return str(path)


def _diff(tmp_path, monkeypatch, lhs, rhs):
    summary = tmp_path / "summary.md"
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(summary))
    monkeypatch.setitem(state.gARG, "lhs", _results(tmp_path, "lhs", lhs))
    monkeypatch.setitem(state.gARG, "rhs", _results(tmp_path, "rhs", rhs))
    diff()
    return summary.read_text()


def _all_targets(grind, exec_time):
    return {"pre_process": {"exec": exec_time}, "simulation": {"exec": exec_time, "grind": grind}, "post_process": {"exec": exec_time}}


def test_reports_only_grind_speedup(tmp_path, monkeypatch, capsys):
    text = _diff(tmp_path, monkeypatch, {"ibm": _all_targets(3.0, 1.0)}, {"ibm": _all_targets(2.0, 1.0)})

    assert "| Case | Grind speedup |" in text
    assert "| `ibm` | 1.50 |" in text
    printed = capsys.readouterr().out
    for unwanted in ("Exec", "Pre Process", "Post Process"):
        assert unwanted not in text
        assert unwanted not in printed


def test_slower_exec_time_alone_does_not_warn(tmp_path, monkeypatch):
    text = _diff(tmp_path, monkeypatch, {"ibm": _all_targets(3.0, 1.0)}, {"ibm": _all_targets(3.0, 9.0)})

    assert "Below threshold" not in text


def test_slower_grind_time_warns(tmp_path, monkeypatch):
    text = _diff(tmp_path, monkeypatch, {"ibm": _all_targets(3.0, 1.0)}, {"ibm": _all_targets(4.0, 1.0)})

    assert "grind speedup 0.75 < 0.95 in `ibm`" in text


def test_needs_only_the_simulation_summary(tmp_path, monkeypatch):
    simulation_only = {"simulation": {"grind": 3.0}}

    text = _diff(tmp_path, monkeypatch, {"ibm": simulation_only}, {"ibm": simulation_only})

    assert "| `ibm` | 1.00 |" in text


@pytest.mark.parametrize("broken", [{}, {"simulation": {}}, {"simulation": {"grind": float("nan")}}])
def test_missing_or_non_finite_grind_time_fails(tmp_path, monkeypatch, broken):
    with pytest.raises(MFCException):
        _diff(tmp_path, monkeypatch, {"ibm": _all_targets(3.0, 1.0)}, {"ibm": broken})
