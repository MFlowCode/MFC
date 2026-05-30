import tempfile
import types as _types
from pathlib import Path
from unittest.mock import patch

from mfc.test.coverage import format_summary, get_changed_files, is_always_run_all, load_map, map_health, param_hash, save_map, select_tests


def test_param_hash_is_order_independent():
    a = param_hash({"m": 100, "weno_order": 5, "bubbles_euler": "T"})
    b = param_hash({"weno_order": 5, "bubbles_euler": "T", "m": 100})
    assert a == b


def test_param_hash_changes_with_value():
    a = param_hash({"weno_order": 5})
    b = param_hash({"weno_order": 3})
    assert a != b


def test_param_hash_is_hex_and_short():
    h = param_hash({"m": 1})
    assert len(h) == 16 and all(c in "0123456789abcdef" for c in h)


def test_param_hash_nested_order_independent():
    a = param_hash({"patch": {"x": 1, "y": 2}})
    b = param_hash({"patch": {"y": 2, "x": 1}})
    assert a == b


def test_save_then_load_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "m.json.gz"
        save_map(p, {"abc": ["src/simulation/m_rhs.fpp"]}, n_tests=1, git_sha="deadbee", gfortran_version="13")
        entries, meta = load_map(p)
        assert entries == {"abc": ["src/simulation/m_rhs.fpp"]}
        assert meta["n_tests"] == 1 and meta["git_sha"] == "deadbee"
        assert "built_at" in meta


def test_load_missing_returns_none():
    assert load_map(Path("/nonexistent/m.json.gz")) == (None, None)


def test_load_corrupt_returns_none():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "m.json.gz"
        p.write_bytes(b"not gzip")
        assert load_map(p) == (None, None)


def test_macro_file_forces_all():
    assert is_always_run_all({"src/common/include/parallel_macros.fpp"})


def test_cmake_forces_all():
    assert is_always_run_all({"CMakeLists.txt"})
    assert is_always_run_all({"toolchain/cmake/foo.cmake"})


def test_param_codegen_forces_all():
    assert is_always_run_all({"toolchain/mfc/params/definitions.py"})


def test_ordinary_common_module_does_not_force_all():
    assert not is_always_run_all({"src/common/m_helper.fpp"})


def test_ordinary_sim_module_does_not_force_all():
    assert not is_always_run_all({"src/simulation/m_rhs.fpp"})


class _Case:
    def __init__(self, ph, params=None):
        self._ph = ph
        self.params = params or {}

    def coverage_key(self):
        return self._ph


def _cases(*phs):
    return [_Case(p) for p in phs]


def test_rung1_no_changed_files_runs_all():
    cases = _cases("a", "b")
    run, skip, reason = select_tests(cases, {"a": ["src/x.fpp"]}, None)
    assert len(run) == 2 and skip == [] and reason.startswith("rung1")


def test_rung2_always_run_all():
    cases = _cases("a", "b")
    run, skip, reason = select_tests(cases, {"a": [], "b": []}, {"CMakeLists.txt"})
    assert len(run) == 2 and reason.startswith("rung2")


def test_rung3_f90_change_runs_all():
    cases = _cases("a")
    run, skip, reason = select_tests(cases, {"a": []}, {"src/common/m_precision_select.f90"})
    assert len(run) == 1 and reason.startswith("rung3")


def test_rung4_changed_fpp_with_zero_coverage_runs_all():
    cases = _cases("a")
    # m_gpu_only.fpp is covered by no test in the map
    run, skip, reason = select_tests(cases, {"a": ["src/simulation/m_rhs.fpp"]}, {"src/simulation/m_gpu_only.fpp"})
    assert len(run) == 1 and reason.startswith("rung4")


def test_rung5_unmapped_test_is_included():
    cases = _cases("a", "new")  # 'new' not in map
    run, skip, _ = select_tests(cases, {"a": ["src/simulation/m_rhs.fpp"]}, {"src/simulation/m_rhs.fpp"})
    assert {c.coverage_key() for c in run} == {"a", "new"}


def test_rung6_and_7_overlap_selects_subset():
    cases = _cases("hit", "miss")
    cov = {"hit": ["src/simulation/m_bubbles_EE.fpp"], "miss": ["src/simulation/m_rhs.fpp"]}
    run, skip, _ = select_tests(cases, cov, {"src/simulation/m_bubbles_EE.fpp"})
    assert [c.coverage_key() for c in run] == ["hit"]
    assert [c.coverage_key() for c in skip] == ["miss"]


def test_case_coverage_key_uses_full_params():
    from mfc.test.case import TestCase

    tc = TestCase("1D -> Foo", {"m": 100, "weno_order": 5})
    assert tc.coverage_key() == param_hash(tc.params)


def test_case_coverage_key_changes_with_params():
    from mfc.test.case import TestCase

    a = TestCase("1D -> Foo", {"weno_order": 5})
    b = TestCase("1D -> Foo", {"weno_order": 3})
    assert a.coverage_key() != b.coverage_key()


def test_case_coverage_key_ignores_trace():
    from mfc.test.case import TestCase

    a = TestCase("1D -> Foo", {"m": 100})
    b = TestCase("totally -> different -> trace", {"m": 100})
    assert a.coverage_key() == b.coverage_key()


def test_changed_files_prefers_explicit_list():
    files = get_changed_files("/repo", "master", explicit="src/a.fpp\nsrc/b.fpp\n")
    assert files == {"src/a.fpp", "src/b.fpp"}


def test_changed_files_deepens_then_recovers():
    state = {"deepened": False}

    def fake_run(cmd, **kw):
        sub = cmd[1] if len(cmd) > 1 else ""
        if sub == "fetch":
            state["deepened"] = True
            return _types.SimpleNamespace(returncode=0, stdout="", stderr="")
        if sub == "merge-base":
            return _types.SimpleNamespace(returncode=0 if state["deepened"] else 1, stdout="base\n", stderr="")
        if sub == "diff":
            return _types.SimpleNamespace(returncode=0, stdout="src/x.fpp\n", stderr="")
        return _types.SimpleNamespace(returncode=0, stdout="", stderr="")

    with patch("subprocess.run", fake_run):
        assert get_changed_files("/repo", "master") == {"src/x.fpp"}


def test_changed_files_returns_none_when_unrecoverable():
    def fake_run(cmd, **kw):
        rc = 1 if (len(cmd) > 1 and cmd[1] == "merge-base") else 0
        return _types.SimpleNamespace(returncode=rc, stdout="", stderr="boom")

    with patch("subprocess.run", fake_run):
        assert get_changed_files("/repo", "master") is None


def test_summary_mentions_counts_age_reason():
    s = format_summary(
        ran=47,
        total=610,
        reason="selected 47/610 by coverage overlap",
        meta={"built_at": "2026-05-20T00:00:00+00:00"},
        now="2026-05-29T00:00:00+00:00",
    )
    assert "47/610" in s and "9d" in s and "coverage overlap" in s


def test_summary_handles_missing_meta():
    s = format_summary(
        ran=610,
        total=610,
        reason="rung1: changed-file list unavailable",
        meta=None,
        now="2026-05-29T00:00:00+00:00",
    )
    assert "610/610" in s and "map age unknown" in s


def test_health_ok():
    ok, msg = map_health(
        meta={"built_at": "2026-05-28T00:00:00+00:00", "n_tests": 600},
        current_keys=set(str(i) for i in range(600)),
        mapped_keys=set(str(i) for i in range(580)),
        now="2026-05-29T00:00:00+00:00",
        max_age_days=10,
        min_fraction=0.8,
    )
    assert ok, msg


def test_health_stale_fails():
    ok, msg = map_health(
        meta={"built_at": "2026-05-01T00:00:00+00:00", "n_tests": 600},
        current_keys=set(["a"]),
        mapped_keys=set(["a"]),
        now="2026-05-29T00:00:00+00:00",
        max_age_days=10,
        min_fraction=0.8,
    )
    assert not ok and "stale" in msg.lower()


def test_health_undercoverage_fails():
    ok, msg = map_health(
        meta={"built_at": "2026-05-28T00:00:00+00:00", "n_tests": 10},
        current_keys=set(str(i) for i in range(100)),
        mapped_keys=set(str(i) for i in range(50)),
        now="2026-05-29T00:00:00+00:00",
        max_age_days=10,
        min_fraction=0.8,
    )
    assert not ok and "coverage" in msg.lower()


def test_builder_has_coverage_key_matching_case():
    from mfc.test.case import TestCaseBuilder

    b = TestCaseBuilder(trace="1D -> Foo", mods={"m": 100, "weno_order": 5}, path="", args=[], ppn=1, functor=None)
    assert b.coverage_key() == b.to_case().coverage_key()


def test_rung5_empty_coverage_is_included():
    # a test whose map entry is [] (uncertain) must be RUN, not skipped, on a .fpp change.
    # "anchor" covers the changed .fpp so rung4 passes and we reach the per-test rungs.
    cases = _cases("hasempty", "anchor")
    cov_map = {
        "hasempty": [],
        "anchor": ["src/simulation/m_rhs.fpp"],
    }
    run, skip, _ = select_tests(cases, cov_map, {"src/simulation/m_rhs.fpp"})
    run_keys = {c.coverage_key() for c in run}
    assert "hasempty" in run_keys and skip == []


def test_changed_files_explicit_space_and_comma_separated():
    from mfc.test.coverage import get_changed_files

    assert get_changed_files("/r", "master", explicit="src/a.fpp src/b.fpp") == {"src/a.fpp", "src/b.fpp"}
    assert get_changed_files("/r", "master", explicit="src/a.fpp,src/b.fpp") == {"src/a.fpp", "src/b.fpp"}


def test_sim_include_fpp_forces_all():
    # gcov can't reliably attribute #:include'd files; any src include change runs all.
    assert is_always_run_all({"src/simulation/include/inline_riemann.fpp"})
    assert is_always_run_all({"src/pre_process/include/2dHardcodedIC.fpp"})


def test_empty_explicit_is_uncertainty_not_skipall():
    # An empty/whitespace --changed-files must NOT become an empty set (skip-all under
    # enforce). It falls through to git detection -> None when that fails -> run all.
    def fail_git(cmd, **kw):
        rc = 1 if (len(cmd) > 1 and cmd[1] == "merge-base") else 0
        return _types.SimpleNamespace(returncode=rc, stdout="", stderr="x")

    with patch("subprocess.run", fail_git):
        assert get_changed_files("/r", "master", explicit="") is None
        assert get_changed_files("/r", "master", explicit="  ,  ") is None


def test_run_and_test_infra_force_all():
    assert is_always_run_all({"toolchain/mfc/run/input.py"})
    assert is_always_run_all({"toolchain/mfc/test/case.py"})
    assert is_always_run_all({"toolchain/mfc/test/test.py"})


def test_cases_py_is_not_always_run():
    assert not is_always_run_all({"toolchain/mfc/test/cases.py"})


def test_cases_py_change_runs_new_tests_not_skipall():
    # cases.py-only change must run the NEW/modified tests (rung 5), not skip everything.
    cases = _cases("mapped", "newtest")  # "newtest" absent from map
    run, skip, _ = select_tests(cases, {"mapped": ["src/simulation/m_rhs.fpp"]}, {"toolchain/mfc/test/cases.py"})
    assert [c.coverage_key() for c in run] == ["newtest"]
    assert [c.coverage_key() for c in skip] == ["mapped"]


def test_docs_only_still_skips_all():
    cases = _cases("a")
    run, skip, reason = select_tests(cases, {"a": ["src/x.fpp"]}, {"README.md"})
    assert run == [] and len(skip) == 1 and "rung7" in reason


def test_uppercase_fortran_extension_forces_all():
    cases = _cases("a")
    run, skip, reason = select_tests(cases, {"a": []}, {"src/common/m_x.F90"})
    assert len(run) == 1 and reason.startswith("rung3")


def test_toolchain_py_change_forces_all_except_cases():
    assert is_always_run_all({"toolchain/mfc/case.py"})
    assert is_always_run_all({"toolchain/mfc/build.py"})
    assert is_always_run_all({"toolchain/mfc/common.py"})
    assert not is_always_run_all({"toolchain/mfc/test/cases.py"})


def test_empty_map_with_fpp_change_runs_all_rung4():
    cases = _cases("a", "b")
    run, skip, reason = select_tests(cases, {}, {"src/simulation/m_rhs.fpp"})
    assert len(run) == 2 and reason.startswith("rung4")
